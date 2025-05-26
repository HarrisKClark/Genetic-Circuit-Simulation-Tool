import tkinter as tk
from tkinter import ttk

class ControlPanel(tk.Toplevel):

    def __init__(self, live_panel, app, live_vars):
        super().__init__(app)
        self.app = app
        self.live_panel = live_panel
        self.live_vars = live_vars

        self.transient(app)
        self.attributes('-topmost', True)
        self.lift()

        self.title("Control Panel")
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        frm = ttk.Frame(self, padding=10)
        frm.pack(fill='both', expand=True)

        ttk.Label(frm, text="Inducer Levels:", font=(None, 10, 'bold'))\
            .pack(anchor='w', pady=(0,5))
        self.slider_frame = ttk.Frame(frm)
        self.slider_frame.pack(fill='both', expand=True)
        self._build_sliders()

        btn_frame = ttk.Frame(frm)
        btn_frame.pack(fill='x', pady=(10,0))
        self.pause_btn = ttk.Button(btn_frame, text="Pause",
                                    command=self._toggle_pause)
        self.pause_btn.pack(side='left', padx=(0,5))
        ttk.Button(btn_frame, text="Stop", command=self._stop_live)\
            .pack(side='left')

    def _build_sliders(self):
        for w in self.slider_frame.winfo_children():
            w.destroy()
        for ind_id, var in self.live_vars.items():
            row = ttk.Frame(self.slider_frame)
            row.pack(fill='x', pady=2)
            ttk.Label(row, text=f"Inducer {ind_id}:").pack(side='left')
            ttk.Scale(row, from_=0.0, to=1.0, orient='horizontal',
                      variable=var).pack(side='left', fill='x', expand=True, padx=(5,0))

    def refresh(self, live_vars):
        self.live_vars = live_vars
        self._build_sliders()

    def _toggle_pause(self):
        if getattr(self.live_panel, '_paused', False):
            self.live_panel.resume()
            self.pause_btn.config(text="Pause")
        else:
            self.live_panel.pause()
            self.pause_btn.config(text="Resume")

    def _stop_live(self):
        self.live_panel.stop()
        self.destroy()

    def _on_close(self):
        self.destroy()
