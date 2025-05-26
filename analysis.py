import tkinter as tk
from tkinter import ttk
import numpy as np
import pandas as pd

from simulation import (
    simulate_deterministic,
    simulate_stochastic,
    simulate_dde,
    simulate_stochastic_dde,
)


class AnalysisTools:
    def __init__(self, app):
        self.app = app

    def run(self, tool_name: str):
        if tool_name == "Sensitivity analysis":
            self._open_sensitivity_window()
        elif tool_name == "Bifurcation analysis":
            self._open_bifurcation_window()
        elif tool_name == "Nullcline analysis":
            self._open_nullcline_window()

    def _open_sensitivity_window(self):
        win = tk.Toplevel(self.app)
        win.title("Sensitivity Analysis")
        frm = tk.Frame(win, padx=10, pady=10)
        frm.pack(fill="both", expand=True)

        tk.Label(frm, text="Parameter:").grid(row=0, column=0, sticky="w")
        param = tk.StringVar(value="alpha")
        ttk.OptionMenu(frm, param, "alpha", "alpha", "n", "sigma").grid(row=0, column=1, sticky="ew")

        tk.Label(frm, text="Min value:").grid(row=1, column=0, sticky="w")
        vmin = tk.DoubleVar(value=self.app.param_alpha.get())
        ttk.Entry(frm, textvariable=vmin).grid(row=1, column=1, sticky="ew")

        tk.Label(frm, text="Max value:").grid(row=2, column=0, sticky="w")
        vmax = tk.DoubleVar(value=self.app.param_alpha.get())
        ttk.Entry(frm, textvariable=vmax).grid(row=2, column=1, sticky="ew")

        tk.Label(frm, text="Steps:").grid(row=3, column=0, sticky="w")
        steps = tk.IntVar(value=10)
        ttk.Entry(frm, textvariable=steps).grid(row=3, column=1, sticky="ew")

        def _go():
            win.destroy()
            self._run_sensitivity(param.get(), vmin.get(), vmax.get(), steps.get())

        ttk.Button(frm, text="Run", command=_go).grid(row=4, column=0, columnspan=2, pady=8)
        frm.columnconfigure(1, weight=1)

    def _run_sensitivity(self, pname: str, vmin: float, vmax: float, steps: int):
        app = self.app
        vals = np.linspace(vmin, vmax, steps)
        metric = []

        gvar = getattr(app, f"param_{pname}")
        g_orig = gvar.get()
        per_orig = {nid: p[pname].get() for nid, p in app.node_params.items() if pname in p}

        for v in vals:
            gvar.set(v)
            for nid in per_orig:
                app.node_params[nid][pname].set(v)
            has_delay = any(getattr(e, "delay", None) and e.delay.get() > 0 for e in app.edges)
            stype = app.param_sim_type.get()
            p = {
                "alpha": app.param_alpha.get(),
                "n": app.param_n.get(),
                "basal": app.param_basal.get(),
                "t0": app.param_t0.get(),
                "t1": app.param_t1.get(),
                "steps": app.param_steps.get(),
                "sigma": app.param_sigma.get(),
                "simulate_all": app.param_sim_all.get(),
                "sim_type": stype,
            }
            if has_delay:
                if stype == "stochastic":
                    sol, *_ = simulate_stochastic_dde(app.nodes, app.edges, app.node_params, p)
                else:
                    sol, *_ = simulate_dde(app.nodes, app.edges, app.node_params, p)
            else:
                if stype == "stochastic":
                    sol, *_ = simulate_stochastic(app.nodes, app.edges, app.node_params, p)
                else:
                    sol, *_ = simulate_deterministic(app.nodes, app.edges, app.node_params, p)
            metric.append(sol.y[0].max() if sol.y.size else np.nan)

        gvar.set(g_orig)
        for nid, val in per_orig.items():
            app.node_params[nid][pname].set(val)

        ax = app.analysis_ax
        ax.clear()
        ax.plot(vals, metric, marker="o")
        ax.set_xlabel(pname)
        ax.set_ylabel("Max expression (Gene 1)")
        ax.set_title(f"Sensitivity: {pname}")
        app.analysis_canvas.draw()
        app.notebook.select(app.analysis_tab)

        pd.DataFrame({pname: vals, "max": metric}).to_csv(f"sensitivity_{pname}.csv", index=False)
        tk.messagebox.showinfo("Done", "Sensitivity sweep finished.")

    def _open_bifurcation_window(self):
        win = tk.Toplevel(self.app)
        win.title("Bifurcation Analysis")
        frm = tk.Frame(win, padx=10, pady=10)
        frm.pack(fill="both", expand=True)

        tk.Label(frm, text="Parameter:").grid(row=0, column=0, sticky="w")
        param = tk.StringVar(value="alpha")
        ttk.OptionMenu(frm, param, "alpha", "alpha", "n", "sigma").grid(row=0, column=1, sticky="ew")

        tk.Label(frm, text="Start:").grid(row=1, column=0, sticky="w")
        v0 = tk.DoubleVar(value=self.app.param_alpha.get())
        ttk.Entry(frm, textvariable=v0).grid(row=1, column=1, sticky="ew")

        tk.Label(frm, text="End:").grid(row=2, column=0, sticky="w")
        v1 = tk.DoubleVar(value=self.app.param_alpha.get())
        ttk.Entry(frm, textvariable=v1).grid(row=2, column=1, sticky="ew")

        tk.Label(frm, text="Steps:").grid(row=3, column=0, sticky="w")
        nsteps = tk.IntVar(value=50)
        ttk.Entry(frm, textvariable=nsteps).grid(row=3, column=1, sticky="ew")

        tk.Label(frm, text="Gene:").grid(row=4, column=0, sticky="w")
        genes = [f"Gene {n.id}" for n in self.app.nodes if n.type == "gene"]
        genes = ["All"] + genes
        gsel = tk.StringVar(value=genes[0])
        ttk.OptionMenu(frm, gsel, *genes).grid(row=4, column=1, sticky="ew")

        def _go():
            win.destroy()
            self._run_bifurcation(param.get(), v0.get(), v1.get(), nsteps.get(), gsel.get())

        ttk.Button(frm, text="Run", command=_go).grid(row=5, column=0, columnspan=2, pady=8)
        frm.columnconfigure(1, weight=1)

    def _run_bifurcation(self, pname, v0, v1, steps, gsel):
        app = self.app
        vals = np.linspace(v0, v1, steps)
        genes = [n for n in app.nodes if n.type == "gene"]
        if gsel == "All":
            idxs = range(len(genes))
        else:
            gid = int(gsel.split()[-1])
            idxs = [next(i for i, n in enumerate(genes) if n.id == gid)]
        curves = {i: [] for i in idxs}

        gvar = getattr(app, f"param_{pname}")
        g_orig = gvar.get()
        per_orig = {nid: p[pname].get() for nid, p in app.node_params.items() if pname in p}

        for v in vals:
            gvar.set(v)
            for nid in per_orig: app.node_params[nid][pname].set(v)
            p = {
                "alpha": app.param_alpha.get(),
                "n": app.param_n.get(),
                "basal": app.param_basal.get(),
                "t0": app.param_t0.get(),
                "t1": app.param_t1.get(),
                "steps": app.param_steps.get(),
                "sigma": app.param_sigma.get(),
                "simulate_all": app.param_sim_all.get(),
                "sim_type": "deterministic",
            }
            sol, *_ = simulate_deterministic(app.nodes, app.edges, app.node_params, p)
            for i in idxs:
                curves[i].append(sol.y[i, -1])

        gvar.set(g_orig)
        for nid, val in per_orig.items(): app.node_params[nid][pname].set(val)

        ax = app.analysis_ax
        ax.clear()
        for i in idxs:
            ax.plot(vals, curves[i], ".", label=f"Gene {genes[i].id}")
        ax.set_xlabel(pname)
        ax.set_ylabel("Steady state")
        ax.set_title(f"Bifurcation: {pname}")
        ax.legend()
        app.analysis_canvas.draw()
        app.notebook.select(app.analysis_tab)

        df = pd.DataFrame({pname: vals, **{f"g{genes[i].id}": curves[i] for i in idxs}})
        df.to_csv(f"bifurcation_{pname}.csv", index=False)
        tk.messagebox.showinfo("Done", "Bifurcation finished")

    def _open_nullcline_window(self):
        win = tk.Toplevel(self.app)
        win.title("Nullcline Analysis")
        frm = tk.Frame(win, padx=10, pady=10)
        frm.pack(fill="both", expand=True)

        tk.Label(frm, text="Gene:").grid(row=0, column=0, sticky="w")
        genes = [f"Gene {n.id}" for n in self.app.nodes if n.type == "gene"] or ["None"]
        gsel = tk.StringVar(value=genes[0])
        ttk.OptionMenu(frm, gsel, *genes).grid(row=0, column=1, sticky="ew")

        tk.Label(frm, text="X min:").grid(row=1, column=0, sticky="w")
        xmin = tk.DoubleVar(value=0.0)
        tk.Entry(frm, textvariable=xmin).grid(row=1, column=1, sticky="ew")

        tk.Label(frm, text="X max:").grid(row=2, column=0, sticky="w")
        xmax = tk.DoubleVar(value=10.0)
        tk.Entry(frm, textvariable=xmax).grid(row=2, column=1, sticky="ew")

        tk.Label(frm, text="Steps:").grid(row=3, column=0, sticky="w")
        steps = tk.IntVar(value=200)
        tk.Entry(frm, textvariable=steps).grid(row=3, column=1, sticky="ew")

        def _go():
            win.destroy()
            self._run_nullcline(gsel.get(), xmin.get(), xmax.get(), steps.get())

        ttk.Button(frm, text="Run", command=_go) \
            .grid(row=4, column=0, columnspan=2, pady=8)

        frm.columnconfigure(1, weight=1)

    def _run_nullcline(self, gene_str, xmin, xmax, steps):
        import numpy as np

        app = self.app
        if not app.nodes:
            tk.messagebox.showwarning("No network", "Create a network first.")
            return

        gid = int(gene_str.split()[-1])
        params = app.node_params[gid]
        alpha = params["alpha"].get()
        hill = params["n"].get()
        deg = params["deg"].get()
        basal = app.param_basal.get()

        act_w, rep_w = [], []
        for e in app.edges:
            if e.source.type == "gene" and e.target.id == gid and e.target.type == "gene":
                w = e.weight.get() if hasattr(e, "weight") else 1.0
                (act_w if e.type == "activator" else rep_w).append(w)

        xs = np.linspace(xmin, xmax, steps)
        prod = np.array([
            basal + alpha *
            (np.prod([(x ** hill / (1 + x ** hill)) ** w for w in act_w]) if act_w else 1.0) *
            (np.prod([(1 / (1 + x ** hill)) ** w for w in rep_w]) if rep_w else 1.0)
            for x in xs
        ])
        deg_line = deg * xs

        diff = prod - deg_line
        roots = []
        for i in range(len(diff) - 1):
            if diff[i] == 0 or diff[i] * diff[i + 1] < 0:
                x0, x1 = xs[i], xs[i + 1]
                y0, y1 = diff[i], diff[i + 1]
                roots.append(x0 - y0 * (x1 - x0) / (y1 - y0))

        ax = app.analysis_ax
        ax.clear()
        ax.plot(xs, prod, label="Production")
        ax.plot(xs, deg_line, label="Degradation")
        if roots:
            y_roots = [
                basal + alpha *
                (np.prod([(r ** hill / (1 + r ** hill)) ** w for w in act_w]) if act_w else 1.0) *
                (np.prod([(1 / (1 + r ** hill)) ** w for w in rep_w]) if rep_w else 1.0)
                for r in roots
            ]
            ax.scatter(roots, y_roots, color="k", zorder=5, label="Fixed points")
        ax.set_xlabel("Expression level x")
        ax.set_ylabel("Rate")
        ax.set_title(f"Nullcline for Gene {gid}")
        ax.legend()
        app.analysis_canvas.draw()
        app.notebook.select(app.analysis_tab)



