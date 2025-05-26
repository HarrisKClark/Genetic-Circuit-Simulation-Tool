import time, threading, tkinter as tk
from typing import Optional
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import threading
import time

from simulation import simulate_step
from control_panel import ControlPanel


class LiveSimulationPanel(tk.Frame):
    def __init__(self, master, app,
                 t_start=0.0, t_end=None, window=100.0,
                 live_control=False, dt=0.5, refresh_ms=100):

        super().__init__(master)
        self.app = app

        self.dt = float(dt)
        self.refresh_ms = max(80, int(refresh_ms))
        self._steps_per_refresh = max(1, int(round((self.refresh_ms / 1000) / self.dt)))
        self.t_end = t_end
        self.window = float(window)
        self.buffer_max = 2_000

        self.genes = [n for n in app.nodes if n.type == "gene"]
        self.gene_index = {n.id: i for i, n in enumerate(self.genes)} 

        if not app.param_sim_all.get():
            reps = [n for n in app.nodes if n.type == "reporter"]
            self.plot_nodes = reps if reps else self.genes  
        else:
            self.plot_nodes = self.genes + [n for n in app.nodes if n.type == "reporter"]

        self.expr_idx = []
        for n in self.plot_nodes:
            if n.type == "gene":
                self.expr_idx.append(self.gene_index[n.id])
            else:
                src_gene = next((e.source for e in app.edges
                                 if e.target is n and e.source.type == "gene"), None)
                self.expr_idx.append(self.gene_index[src_gene.id] if src_gene else 0)

        self.fig = Figure(dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("t");
        self.ax.set_ylabel("Expression")

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self._background = None
        self.canvas.mpl_connect(
            "draw_event",
            lambda ev: setattr(self, "_background",
                               self.canvas.copy_from_bbox(self.ax.bbox)))

        self._t = [t_start]
        self._ys = [np.zeros(len(self.genes))]
        self.lines = []
        for n in self.plot_nodes:
            col = (app.reporter_color_map.get(app.node_params[n.id]
                                              .get("reporter", "").get(), "black")
                   if n.type == "reporter" else None)
            ln, = self.ax.plot([], [], label=f"{n.type.capitalize()} {n.id}",
                               color=col)
            self.lines.append(ln)
        if len(self.lines) > 1:
            self.ax.legend(loc="upper right")

        self._live_vars = {}
        self._control_win = None
        self._spawn_or_refresh_panel()

        self._running = True
        self._paused = False
        threading.Thread(target=self._loop, daemon=True).start()

    def _cache_background(self):
        self._background = self.canvas.copy_from_bbox(self.ax.bbox)

    def _initial_state(self):
        return np.array([self.app.node_params[n.id]["ic"].get()
                         for n in self.genes], dtype=float)

    def _loop(self):
        t, y = self._t[0], self._ys[0]

        while self._running and (self.t_end is None or t < self.t_end):
            if getattr(self, "_paused", False):
                time.sleep(0.05)
                continue

            for ind_id, var in self._live_vars.items():
                p = self.app.node_params.get(ind_id)
                if p and "level" in p:
                    p["level"].set(var.get())

            for _ in range(self._steps_per_refresh):
                y = simulate_step(
                    self.app.nodes, self.app.edges, self.app.node_params,
                    y, t, t + self.dt,
                    dt=self.dt,
                    sigma=self.app.param_sigma.get(),
                    sim_type=self.app.param_sim_type.get(),
                    alpha_global=self.app.param_alpha.get(),
                    n_global=self.app.param_n.get(),
                    basal_global=self.app.param_basal.get())

                t += self.dt
                self._t.append(t)
                self._ys.append(y.copy())

            self._t.append(t)
            self._ys.append(y.copy())
            if len(self._t) > self.buffer_max:
                trim = len(self._t) - self.buffer_max
                self._t = self._t[trim:]
                self._ys = self._ys[trim:]

            self.after(0, self._redraw)
            time.sleep(self.refresh_ms / 1000.0)

    def _redraw(self):
        import numpy as np

        t_arr = np.asarray(self._t)
        y_arr = np.vstack(self._ys).T

        for ln, idx in zip(self.lines, self.expr_idx):
            ln.set_data(t_arr, y_arr[idx])

        t_end = t_arr[-1]
        t_start = t_end - self.window
        if t_start < t_arr[0]:
            t_start = t_arr[0]
        self.ax.set_xlim(t_start, t_end)

        y_min = np.nanmin(y_arr)
        y_max = np.nanmax(y_arr)
        if np.isfinite(y_min) and np.isfinite(y_max):
            margin = (y_max - y_min) * 0.1
            if margin == 0:
                margin = max(1.0, abs(y_min) * 0.1)
            self.ax.set_ylim(y_min - margin, y_max + margin)

        self.canvas.draw_idle()

    def _build_control_panel(self):
        if not hasattr(self, "_control_win") or not self._control_win.winfo_exists():
            self._control_win = ControlPanel(self.app, self._live_vars)
        else:
            self._control_win.refresh(self._live_vars)
        return self._control_win

    def refresh_controls(self):
        self._build_control_panel(None)

    def _spawn_or_refresh_panel(self):
        self._live_vars.clear()
        for n in self.app.nodes:
            if n.type == "inducer":
                var = self.app.node_params[n.id].setdefault(
                    "level", tk.DoubleVar(value=1.0)
                )
                self._live_vars[n.id] = var

        has_window = hasattr(self, "_control_win") and \
                     self._control_win is not None and \
                     self._control_win.winfo_exists()

        if not has_window:
            self._control_win = ControlPanel(self, self.app, self._live_vars)
            self._control_win.attributes("-topmost", True)
            self._control_win.lift()
        else:
            self._control_win.refresh(self._live_vars)

    def stop(self):
        self._running = False
        if hasattr(self, "_panel") and self._panel.winfo_exists():
            self._panel.destroy()

    def pause(self):   self._paused = True
    def resume(self):  self._paused = False

    def stop(self):
        self._running = False
        try:
            self.ctrl.destroy()
        except Exception:
            pass
