import tkinter as tk
from tkinter import ttk, BooleanVar
import re
import math
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from nodes import Node
from edges import Edge
from simulation import simulate_network
from analysis_tools import AnalysisTools
from live_simulation import LiveSimulationPanel


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Gene Regulatory Network Editor")

        self.param_alpha = tk.DoubleVar(value=10.0)
        self.param_n = tk.DoubleVar(value=2.0)
        self.param_basal = tk.DoubleVar(value=0.1)
        self.param_t0 = tk.DoubleVar(value=0.0)
        self.param_t1 = tk.DoubleVar(value=50.0)
        self.param_steps = tk.IntVar(value=500)
        self.param_sigma = tk.DoubleVar(value=0.1)
        self.param_sim_type = tk.StringVar(value='deterministic')
        self.param_sim_all = BooleanVar(value=False)


        self.node_type = tk.StringVar(value='gene')
        self.edge_type = tk.StringVar(value='activator')
        self.analysis_tool = tk.StringVar(value='Sensitivity analysis')

        self._create_ui()

        self.reporter_color_map = {
            'GFP': 'green',
            'RFP': 'red',
            'BFP': 'blue',
            'YFP': 'yellow',
            'CFP': 'cyan',
            'mCherry': 'magenta'
        }

        outer = tk.PanedWindow(
            self,
            orient=tk.HORIZONTAL,
            sashwidth=8,  # thickness of the divider
            sashrelief=tk.RAISED,
            sashpad=4  # extra grab padding
        )
        outer.pack(fill='both', expand=True)

        # LEFT pane
        left = tk.Frame(outer)
        outer.add(left, minsize=200)
        self.tree = ttk.Treeview(left, show='tree')
        self.tree.pack(fill='both', expand=True)
        for iid, lbl in (
                ('genes_root', 'Genes'),
                ('inducers_root', 'Inducers'),
                ('reporters_root', 'Reporters'),
                ('edges_root', 'Edges')  
        ):
            self.tree.insert('', 'end', iid=iid, text=lbl, open=True)
        self.tree.bind('<<TreeviewSelect>>', self._on_tree_select)

        self.node_settings_frame = tk.LabelFrame(left,
                                                 text="Node Settings", padx=5, pady=5)
        self.node_settings_frame.pack(fill='x', padx=5, pady=5)

        right = tk.PanedWindow(
            outer,
            orient=tk.VERTICAL,
            sashwidth=8,
            sashrelief=tk.RAISED,
            sashpad=4
        )
        outer.add(right, stretch='always')

        editor_frame = tk.Frame(right)
        right.add(editor_frame, minsize=200)
        self.canvas = tk.Canvas(editor_frame, bg='white')
        self.canvas.pack(fill='both', expand=True)

        sep = ttk.Separator(editor_frame, orient='horizontal')
        sep.pack(fill='x', pady=2)

        plot_frame = tk.Frame(right)
        right.add(plot_frame, minsize=150)
        self.plot_frame = plot_frame

        self.notebook = ttk.Notebook(self.plot_frame)
        self.notebook.pack(fill='both', expand=True)

        self.sim_tab = tk.Frame(self.notebook)
        self.notebook.add(self.sim_tab, text='Simulation')

        self.analysis_tab = tk.Frame(self.notebook)
        self.notebook.add(self.analysis_tab, text='Analysis')

        self.sim_fig = Figure(dpi=100)
        self.sim_fig.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.10)
        self.sim_ax = self.sim_fig.add_subplot(111)
        self.sim_canvas = FigureCanvasTkAgg(self.sim_fig, master=self.sim_tab)
        self.sim_canvas.draw()
        self.sim_canvas.get_tk_widget().pack(fill='both', expand=True)
        sim_toolbar = NavigationToolbar2Tk(self.sim_canvas, self.sim_tab)
        sim_toolbar.update()
        self.sim_canvas._tkcanvas.pack(fill='x')

        self.sim_tab.bind("<Configure>", self._on_sim_tab_resize)

        self.analysis_fig = Figure(dpi=100)
        self.analysis_fig.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.10)
        self.analysis_ax = self.analysis_fig.add_subplot(111)

        self.analysis_canvas = FigureCanvasTkAgg(self.analysis_fig, master=self.analysis_tab)
        self.analysis_canvas.draw()
        self.analysis_canvas.get_tk_widget().pack(fill='both', expand=True)
        ana_toolbar = NavigationToolbar2Tk(self.analysis_canvas, self.analysis_tab)
        ana_toolbar.update()

        self.analysis_tab.bind("<Configure>", self._on_analysis_tab_resize)

        self.nodes = []
        self.node_params = {}
        self.edges = []
        self.node_count = 0
        self.edge_count = 0

        self.selected_node = None

        self.canvas.bind('<Button-1>', self._on_click)
        self.canvas.bind('<Button-2>', self._on_right_click)
        self.canvas.bind('<Button-3>', self._on_right_click)

        self.analysis = AnalysisTools(self)

    def _create_ui(self):
        frame = tk.Frame(self)
        frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        tk.Label(frame, text="Node Type:").pack(side=tk.LEFT)
        tk.OptionMenu(frame, self.node_type, "gene", "inducer", "reporter").pack(side=tk.LEFT)
        tk.Label(frame, text="Edge Type:").pack(side=tk.LEFT, padx=(10, 0))
        tk.OptionMenu(frame, self.edge_type, "activator", "repressor").pack(side=tk.LEFT)

        tk.Button(frame, text="Declutter", command=self._declutter).pack(side=tk.LEFT, padx=5)

        self.sim_mode = tk.StringVar(value="Static simulation")
        tk.Label(frame, text="Simulation:").pack(side=tk.LEFT, padx=(20, 0))
        tk.OptionMenu(frame, self.sim_mode, "Static simulation", "Live simulation").pack(side=tk.LEFT)
        tk.Button(frame, text="Run", command=self._open_simulation_window).pack(side=tk.LEFT, padx=4)

        tk.Label(frame, text="Analysis tool:").pack(side=tk.LEFT, padx=(20, 0))
        tk.OptionMenu(frame, self.analysis_tool, "Sensitivity analysis", "Bifurcation analysis", "Nullcline analysis").pack(side=tk.LEFT)
        tk.Button(frame, text="Run", command=self._run_analysis).pack(side=tk.LEFT, padx=4)

    def _open_simulation_window(self):
        if self.sim_mode.get() == "Static simulation":
            self._open_static_window()
        else:
            self._open_live_window()

    def _open_static_window(self):
        win = tk.Toplevel(self);
        win.title("Static-simulation options")
        frm = ttk.Frame(win, padding=10);
        frm.pack(fill="both", expand=True)

        rows = [
            ("t₀:", self.param_t0),
            ("t₁:", self.param_t1),
            ("Steps:", self.param_steps),
        ]
        for r, (lbl, var) in enumerate(rows):
            ttk.Label(frm, text=lbl).grid(row=r, column=0, sticky="w")
            ttk.Entry(frm, textvariable=var, width=8).grid(row=r, column=1, sticky="ew")

        ttk.Label(frm, text="Sim-type:").grid(row=3, column=0, sticky="w")
        ttk.OptionMenu(frm, self.param_sim_type,
                       self.param_sim_type.get(), "deterministic", "stochastic") \
            .grid(row=3, column=1, sticky="ew")
        ttk.Label(frm, text="σ (noise):").grid(row=4, column=0, sticky="w")
        sigma = ttk.Entry(frm, textvariable=self.param_sigma, width=8)
        sigma.grid(row=4, column=1, sticky="ew")

        # toggle σ widget
        sigma.state(["!disabled"] if self.param_sim_type.get() == "stochastic"
                    else ["disabled"])
        self.param_sim_type.trace_add("write", lambda *_: sigma.state(["!disabled"] if self.param_sim_type.get() == "stochastic" else ["disabled"]))

        ttk.Checkbutton(frm, text="Simulate ALL nodes",
                        variable=self.param_sim_all).grid(row=5, column=0,
                                                          columnspan=2, sticky="w")
        ttk.Button(frm, text="Run",
                   command=lambda: (win.destroy(), self._simulate())) \
            .grid(row=6, column=0, columnspan=2, pady=8)
        frm.columnconfigure(1, weight=1)

    def _open_live_window(self):
        win = tk.Toplevel(self);
        win.title("Live-simulation options")
        frm = ttk.Frame(win, padding=10);
        frm.pack(fill="both", expand=True)

        win_len = tk.DoubleVar(value=self.param_steps.get()) 
        dt_var = tk.DoubleVar(value=0.5)
        ref_var = tk.IntVar(value=100)

        rows = [("Visible window (t):", win_len),
                ("Integrator dt:", dt_var),
                ("Refresh (ms):", ref_var)]
        for r, (lbl, var) in enumerate(rows):
            ttk.Label(frm, text=lbl).grid(row=r, column=0, sticky="w")
            ttk.Entry(frm, textvariable=var, width=8).grid(row=r, column=1, sticky="ew")

        live_ctl = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm, text="Show inducer control panel",
                        variable=live_ctl).grid(row=3, column=0,
                                                columnspan=2, sticky="w")

        def _run():
            win.destroy()
            for w in self.sim_tab.winfo_children(): w.destroy()
            self.live_panel = LiveSimulationPanel(
                self.sim_tab, self,
                window=win_len.get(),
                dt=dt_var.get(),
                refresh_ms=ref_var.get(),
                live_control=live_ctl.get())
            self.live_panel.pack(fill="both", expand=True)

        ttk.Button(frm, text="Run", command=_run) \
            .grid(row=4, column=0, columnspan=2, pady=8)
        frm.columnconfigure(1, weight=1)

    def _open_settings_window(self):
        win = tk.Toplevel(self)
        win.title("Global Node Parameters")
        frm = ttk.Frame(win, padding=10)
        frm.pack(fill="both", expand=True)

        row_vars = [("α (production rate):", self.param_alpha),
                    ("n (Hill coefficient):", self.param_n),
                    ("basal (leak):", self.param_basal)]
        for i, (lbl, var) in enumerate(row_vars):
            ttk.Label(frm, text=lbl).grid(row=i, column=0, sticky="w", pady=2)
            ttk.Entry(frm, textvariable=var, width=10).grid(row=i, column=1, sticky="ew")

        def _apply_all():
            for nid, p in self.node_params.items():
                if "alpha" in p:
                    p["alpha"].set(self.param_alpha.get())
                    p["n"].set(self.param_n.get())

        ttk.Button(frm, text="Apply to ALL genes", command=_apply_all) \
            .grid(row=len(row_vars), column=0, columnspan=2, pady=(10, 0))

        frm.columnconfigure(1, weight=1)

    def _add_node(self, x, y):
        self.node_count += 1
        ntype = self.node_type.get()  
        node = Node(self.canvas, x, y, self.node_count, node_type=ntype)
        self.nodes.append(node)

        if ntype == 'gene':
            self.node_params[node.id] = {
                'ic': tk.DoubleVar(value=0.0),
                'alpha': tk.DoubleVar(value=self.param_alpha.get()),
                'n': tk.DoubleVar(value=self.param_n.get()),
                'deg': tk.DoubleVar(value=1.0),
            }
            parent, text = 'genes_root', f"Gene {node.id}"
        elif ntype == 'inducer':
            self.node_params[node.id] = {
                'on_time': tk.DoubleVar(value=1.0)
            }
            parent, text = 'inducers_root', f"Inducer {node.id}"
        else:  # reporter
            self.node_params[node.id] = {
                'reporter': tk.StringVar(value='GFP')
            }
            parent, text = 'reporters_root', f"Reporter {node.id}"
            self.canvas.itemconfig(node.oval,
                                   fill=self.reporter_color_map['GFP'])

        self.tree.insert(parent, 'end', iid=str(node.id), text=text)

    def _on_click(self, event):
        x, y = event.x, event.y
        clicked = self.canvas.find_overlapping(x, y, x, y)
        node = None
        for item in clicked:
            for tag in self.canvas.gettags(item):
                m = re.fullmatch(r'node(\d+)', tag)
                if m:
                    nid = int(m.group(1))
                    node = next(n for n in self.nodes if n.id == nid)
                    break
            if node:
                break
        if node:
            self._handle_node_click(node)
        else:
            self._add_node(x, y)

    def _on_right_click(self, event):
        x, y = event.x, event.y
        clicked_items = self.canvas.find_overlapping(x, y, x, y)
        node = None
        for item in clicked_items:
            for tag in self.canvas.gettags(item):
                m = re.fullmatch(r'node(\d+)', tag)
                if m:
                    nid = int(m.group(1))
                    node = next(n for n in self.nodes if n.id == nid)
                    break
            if node:
                break

        if node:
            removed = [e for e in self.edges
                       if e.source is node or e.target is node]

            self.canvas.delete(node.oval)
            self.canvas.delete(node.label)

            self.nodes.remove(node)
            self.node_params.pop(node.id, None)

            self.tree.delete(str(node.id))

            for e in removed:

                self.tree.delete(f"edge{e.id}")
            self.edges = [e for e in self.edges if e not in removed]

            for w in self.node_settings_frame.winfo_children():
                w.destroy()

            self.canvas.delete('edge')
            for e in self.edges:
                e.draw()

    def _on_tree_select(self, event):
        sel = self.tree.selection()
        if not sel: return
        iid = sel[0]
        if iid in ('genes_root', 'inducers_root', 'reporters_root', 'edges_root'):
            for w in self.node_settings_frame.winfo_children():
                w.destroy()
            return
        if iid.startswith('edge'):
            edge_id = int(iid.replace('edge', ''))
            self._build_edge_settings(edge_id)
        else:
            node_id = int(iid)
            self._build_node_settings(node_id)

    def _build_node_settings(self, node_id):
        for w in self.node_settings_frame.winfo_children():
            w.destroy()

        node = next(n for n in self.nodes if n.id == node_id)
        params = self.node_params[node.id]

        tk.Label(self.node_settings_frame,
                 text=f"{node.type.capitalize()} {node.id} settings",
                 font=("TkDefaultFont", 10, "bold")) \
            .pack(anchor="w", pady=(0, 4))

        if node.type == "gene":
            tk.Label(self.node_settings_frame, text="Initial condition:").pack(anchor="w")
            tk.Entry(self.node_settings_frame, textvariable=params["ic"],
                     width=10).pack(anchor="w", pady=(0, 4))

            for lbl, key in (("α  (production rate):", "alpha"),
                             ("n  (Hill coeff):", "n"),
                             ("deg (degradation):", "deg")):
                tk.Label(self.node_settings_frame, text=lbl).pack(anchor="w")
                tk.Entry(self.node_settings_frame, textvariable=params[key],
                         width=10).pack(anchor="w", pady=(0, 4))

            def _apply_all_genes():
                for p in self.node_params.values():
                    if "alpha" in p: p["alpha"].set(params["alpha"].get())
                    if "n" in p: p["n"].set(params["n"].get())
                    if "deg" in p: p["deg"].set(params["deg"].get())

            tk.Button(self.node_settings_frame, text="Apply to ALL genes",
                      command=_apply_all_genes) \
                .pack(anchor="w", pady=(6, 0))

        elif node.type == "inducer":
            if "level" not in params: params["level"] = tk.DoubleVar(value=1.0)
            if "on_time" not in params: params["on_time"] = tk.DoubleVar(value=0.0)

            tk.Label(self.node_settings_frame, text="Level (0 – 1):").pack(anchor="w")
            tk.Scale(self.node_settings_frame, from_=0, to=1,
                     orient="horizontal", variable=params["level"],
                     length=120).pack(anchor="w", pady=(0, 6))

            tk.Label(self.node_settings_frame, text="On-time (t₀ units):").pack(anchor="w")
            tk.Entry(self.node_settings_frame, textvariable=params["on_time"],
                     width=10).pack(anchor="w")

        else:  
            tk.Label(self.node_settings_frame, text="Reporter type:").pack(anchor="w")
            var = params["reporter"]
            ttk.OptionMenu(self.node_settings_frame, var, var.get(),
                           *self.reporter_color_map.keys(),
                           command=lambda v, nid=node.id:
                           self._set_reporter_color(nid, v)) \
                .pack(anchor="w")

    def _handle_node_click(self, node):
        if not self.selected_node:
            self.selected_node = node
            self.canvas.itemconfig(node.oval, outline='red')
        else:
            et = self.edge_type.get()
            edge = Edge(self.canvas, self.selected_node, node, et)
            self.edge_count += 1
            edge.id = self.edge_count
            edge.weight = tk.DoubleVar(value=1.0)
            edge.delay = tk.DoubleVar(value=0.0)
            self.edges.append(edge)
            txt = f"{edge.source.type.capitalize()}{edge.source.id}→{edge.target.type.capitalize()}{edge.target.id}"
            self.tree.insert('edges_root', 'end',
                             iid=f"edge{edge.id}", text=txt)
            self.canvas.itemconfig(self.selected_node.oval, outline='black')
            self.selected_node = None

    def _reset_selection(self):
        if self.selected_node:
            self.canvas.itemconfig(self.selected_node.oval, outline='black')
            self.selected_node = None

    def _declutter(self):
        self.canvas.update_idletasks()
        width = self.canvas.winfo_width()
        height = self.canvas.winfo_height()
        n = len(self.nodes)
        if n < 2:
            return

        area = width * height
        k = math.sqrt(area / n)

        pos = {node: [node.x, node.y] for node in self.nodes}
        margin = max(node.radius for node in self.nodes) + 10

        temp = min(width, height) / 10.0
        for _ in range(50):
            disp = {node: [0.0, 0.0] for node in self.nodes}

            for i, v in enumerate(self.nodes):
                for u in self.nodes[i + 1:]:
                    dx = pos[v][0] - pos[u][0]
                    dy = pos[v][1] - pos[u][1]
                    dist = math.hypot(dx, dy) or 0.01
                    force = k * k / dist
                    disp[v][0] += dx / dist * force
                    disp[v][1] += dy / dist * force
                    disp[u][0] -= dx / dist * force
                    disp[u][1] -= dy / dist * force

            for e in self.edges:
                if e.source not in self.nodes or e.target not in self.nodes:
                    continue
                sx, sy = pos[e.source]
                tx, ty = pos[e.target]
                dx = sx - tx
                dy = sy - ty
                dist = math.hypot(dx, dy) or 0.01
                force = dist * dist / k
                disp[e.source][0] -= dx / dist * force
                disp[e.source][1] -= dy / dist * force
                disp[e.target][0] += dx / dist * force
                disp[e.target][1] += dy / dist * force

            width = self.canvas.winfo_width()
            height = self.canvas.winfo_height()
            for node in self.nodes:
                dx, dy = disp[node]
                d = math.hypot(dx, dy)
                if d:
                    step = min(d, temp)
                    pos[node][0] += dx / d * step
                    pos[node][1] += dy / d * step
                pos[node][0] = min(max(pos[node][0], margin), width - margin)
                pos[node][1] = min(max(pos[node][1], margin), height - margin)

            temp *= 0.95

        for node in self.nodes:
            node.x, node.y = pos[node]
            self.canvas.coords(
                node.oval,
                node.x - node.radius, node.y - node.radius,
                node.x + node.radius, node.y + node.radius
            )
            self.canvas.coords(node.label, node.x, node.y)

        self.canvas.delete('edge')
        for e in self.edges:
            e.draw()

    def _ensure_static_canvas(self):
        try:
            if self.sim_canvas.get_tk_widget().winfo_exists():
                return  
        except Exception:
            pass  

        self.sim_fig = Figure(dpi=100)
        self.sim_ax = self.sim_fig.add_subplot(111)
        self.sim_canvas = FigureCanvasTkAgg(self.sim_fig, master=self.sim_tab)
        self.sim_canvas.draw()
        self.sim_canvas.get_tk_widget().pack(fill="both", expand=True)

        tb = NavigationToolbar2Tk(self.sim_canvas, self.sim_tab)
        tb.update()
        self.sim_canvas._tkcanvas.pack(fill="x")

    def _simulate(self):
        self._ensure_static_canvas()
        """Run a static simulation and show it in the Simulation tab."""
        if hasattr(self, "live_panel"):
            self.live_panel.stop()
            self.live_panel.destroy()
            del self.live_panel
            self.sim_canvas.get_tk_widget().pack(fill="both", expand=True)

        if not hasattr(self, "_sim_toolbar_packed"):
            self._sim_toolbar_packed = True
        else:
            for child in self.sim_tab.winfo_children():
                if isinstance(child, NavigationToolbar2Tk):
                    child.pack(fill="x")
                    break
        params = {
            "alpha": self.param_alpha.get(),
            "n": self.param_n.get(),
            "basal": self.param_basal.get(),
            "t0": self.param_t0.get(),
            "t1": self.param_t1.get(),
            "steps": self.param_steps.get(),
            "sigma": self.param_sigma.get(),
            "simulate_all": self.param_sim_all.get(),
            "sim_type": self.param_sim_type.get()
        }

        self.sim_ax.clear()
        simulate_network(self.nodes, self.edges,
                         self.node_params, params, ax=self.sim_ax)
        self.sim_canvas.draw()
        self.notebook.select(self.sim_tab)

    def _on_plot_frame_resize(self, event):
        w, h = event.width, event.height
        dpi = self.fig.get_dpi()
        self.fig.set_size_inches(w / dpi, h / dpi)
        self.canvas_fig.draw()

    def _on_plot_frame_resize(self, event):

        w, h = event.width, event.height

        dpi = self.fig.get_dpi()

        self.fig.set_size_inches(w / dpi, h / dpi)

        self.canvas_fig.draw()

    def _build_edge_settings(self, edge_id):
        """Rebuild the settings panel for a selected edge."""
        for w in self.node_settings_frame.winfo_children():
            w.destroy()

        edge = next(e for e in self.edges if e.id == edge_id)

        tk.Label(self.node_settings_frame, text=f"Edge {edge_id}",
                 font=("TkDefaultFont", 10, "bold")) \
            .pack(anchor="w", pady=(0, 4))

        for lbl, var in (("Weight:", edge.weight),
                         ("Delay:", edge.delay)):
            tk.Label(self.node_settings_frame, text=lbl).pack(anchor="w")
            tk.Entry(self.node_settings_frame, textvariable=var,
                     width=10).pack(anchor="w", pady=(0, 4))

        def _reset():
            edge.weight.set(1.0)
            edge.delay.set(0.0)

        def _apply_all():
            for e in self.edges:
                e.weight.set(edge.weight.get())
                e.delay.set(edge.delay.get())

        row = tk.Frame(self.node_settings_frame)
        row.pack(anchor="w", pady=(6, 0))
        tk.Button(row, text="Reset", command=_reset).pack(side="left", padx=(0, 6))
        tk.Button(row, text="Apply to ALL edges", command=_apply_all).pack(side="left")

    def _reset_node_params(self, node_id):
        p = self.node_params[node_id]
        p['alpha'].set(self.param_alpha.get())
        p['n'].set(self.param_n.get())
        p['deg'].set(1.0)

    def _reset_edge_params(self, edge):
        edge.weight.set(1.0)
        edge.delay.set(0.0)

    def _reset_all_params(self):
        for nid, p in self.node_params.items():
            if 'alpha' in p:
                self._reset_node_params(nid)
        for e in self.edges:
            self._reset_edge_params(e)

    def _set_reporter_color(self, node_id, reporter_type):
        node = next(node for node in self.nodes if node.id == node_id)
        color = self.reporter_color_map.get(reporter_type, 'white')
        self.canvas.itemconfig(node.oval, fill=color)

    def _run_analysis(self):
        self.analysis.run(self.analysis_tool.get())

    def _on_sim_tab_resize(self, event):
        w, h = event.width, event.height
        dpi = self.sim_fig.get_dpi()
        self.sim_fig.set_size_inches(w / dpi, h / dpi)
        self.sim_canvas.draw()

    def _on_analysis_tab_resize(self, event):
        w, h = event.width, event.height
        dpi = self.analysis_fig.get_dpi()
        self.analysis_fig.set_size_inches(w / dpi, h / dpi)
        self.analysis_canvas.draw()


if __name__=='__main__':
    app=App()
    app.mainloop()
