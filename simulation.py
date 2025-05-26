import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import tkinter as tk     



def _prepare_network(nodes, edges):
    gene_nodes = [n for n in nodes if n.type == "gene"]
    G = len(gene_nodes)
    idx = {n.id: i for i, n in enumerate(gene_nodes)}

    act_g = {i: [] for i in range(G)}
    rep_g = {i: [] for i in range(G)}
    act_i = {i: [] for i in range(G)}
    rep_i = {i: [] for i in range(G)}

    for e in edges:
        if e.target.type != "gene":
            continue
        gi = idx[e.target.id]
        if e.source.type == "gene":
            sj = idx[e.source.id]
            (act_g if e.type == "activator" else rep_g)[gi].append((sj, e.weight.get()))
        else:
            (act_i if e.type == "activator" else rep_i)[gi].append((e.source.id, e.weight.get()))
    return gene_nodes, idx, act_g, rep_g, act_i, rep_i


def simulate_deterministic(nodes, edges, node_params, params):
    gene_nodes     = [n for n in nodes if n.type == 'gene']
    inducer_nodes  = [n for n in nodes if n.type == 'inducer']
    reporter_nodes = [n for n in nodes if n.type == 'reporter']

    G = len(gene_nodes)
    gene_index = {n.id: i for i, n in enumerate(gene_nodes)}

    act_g = {i: [] for i in range(G)}
    rep_g = {i: [] for i in range(G)}
    act_i = {i: [] for i in range(G)}
    rep_i = {i: [] for i in range(G)}

    on_times = {
        nid: node_params[nid]['on_time'].get()
        for nid in node_params
        if nid in [u.id for u in inducer_nodes]
    }

    for e in edges:
        if e.target.type != 'gene':
            continue
        ti = gene_index[e.target.id]
        wvar = getattr(e, 'weight', None)
        w = wvar.get() if wvar is not None else 1.0
        dvar = getattr(e, 'delay', None)
        d = dvar.get() if dvar is not None else 0.0

        if e.source.type == 'gene':
            si = gene_index.get(e.source.id)
            if si is None:
                continue
            if e.type == 'activator':
                act_g[ti].append((si, w))
            else:
                rep_g[ti].append((si, w))
        else:
            if e.type == 'activator':
                act_i[ti].append((e.source.id, d, w))
            else:
                rep_i[ti].append((e.source.id, d, w))

    basal = params['basal']
    t0, t1, steps = params['t0'], params['t1'], params['steps']

    x0 = np.array([node_params[n.id]['ic'].get() for n in gene_nodes])

    def odes(t, x):
        dx = np.zeros(G)
        for i, node in enumerate(gene_nodes):
            nid = node.id
            p = node_params[nid]
            alpha_i = p['alpha'].get()
            hill_i  = p['n'].get()
            deg_i   = p['deg'].get()

            pa_g = np.prod([ (x[j]**hill_i/(1+x[j]**hill_i))**w for j, w in act_g[i] ]) if act_g[i] else 1.0
            pr_g = np.prod([ (1.0/(1+x[j]**hill_i))**w for j, w in rep_g[i] ]) if rep_g[i] else 1.0
            pa_i = np.prod([ 1.0 if t >= on_times[jid] + d else 0.0
                             for jid, d, w in act_i[i] ]) if act_i[i] else 1.0
            pr_i = np.prod([ 1.0/(1 + (1.0 if t >= on_times[jid] + d else 0.0)**hill_i)**w
                             for jid, d, w in rep_i[i] ]) if rep_i[i] else 1.0

            dx[i] = basal + alpha_i * pa_g * pa_i * pr_g * pr_i - deg_i * x[i]
        return dx

    t_eval = np.linspace(t0, t1, steps)
    sol = solve_ivp(odes, (t0, t1), x0, t_eval=t_eval)
    return sol, gene_nodes, reporter_nodes


def simulate_stochastic(nodes, edges, node_params, params):
    sol_det, gene_nodes, reporter_nodes = simulate_deterministic(nodes, edges, node_params, params)
    alpha, hill, basal = params['alpha'], params['n'], params['basal']
    t0, t1, steps     = params['t0'], params['t1'], params['steps']
    sigma = params.get('sigma', 0.1)
    dt    = (t1 - t0) / max(1, steps - 1)
    G     = len(gene_nodes)

    gene_index = {n.id: i for i,n in enumerate(gene_nodes)}
    act_g = {i: [] for i in range(G)}; rep_g = {i: [] for i in range(G)}
    act_i = {i: [] for i in range(G)}; rep_i = {i: [] for i in range(G)}
    on_times = {nid: node_params[nid]['on_time'].get() for nid in node_params if nid in [u.id for u in nodes if u.type=='inducer']}
    for e in edges:
        tid = e.target.id
        if e.target.type != 'gene': continue
        ti = gene_index[tid]
        w = getattr(e, 'weight', None)
        wval = w.get() if w is not None else 1.0
        d = getattr(e, 'delay', None)
        dval = d.get() if d is not None else 0.0
        if e.source.type == 'gene':
            si = gene_index.get(e.source.id)
            if si is None: continue
            if e.type=='activator': act_g[ti].append((si, wval))
            else:                 rep_g[ti].append((si, wval))
        elif e.source.type=='inducer':
            if e.type=='activator': act_i[ti].append((e.source.id, dval, wval))
            else:                    rep_i[ti].append((e.source.id, dval, wval))

    times = np.linspace(t0, t1, steps)
    x = np.zeros((G, steps))
    x[:,0] = sol_det.y[:,0] if sol_det.y.shape[1] else np.zeros(G)

    def drift(xv, t):
        dx = np.zeros(G)
        for i in range(G):
            pa_g = np.prod([ (xv[j]**hill/(1+xv[j]**hill))**w for j,w in act_g[i] ]) if act_g[i] else 1.0
            pr_g = np.prod([ (1.0/(1+xv[j]**hill))**w       for j,w in rep_g[i] ]) if rep_g[i] else 1.0
            pa_i = np.prod([1.0 if t>= on_times[jid]+d else 0.0 for jid,d,w in act_i[i]]) if act_i[i] else 1.0
            pr_i = np.prod([1.0/(1+(1.0 if t>= on_times[jid]+d else 0.0)**hill)**w for jid,d,w in rep_i[i]]) if rep_i[i] else 1.0
            dx[i] = basal + alpha*pa_g*pa_i*pr_g*pr_i - xv[i]
        return dx

    for k in range(1, steps):
        t = times[k-1]
        d = drift(x[:,k-1], t)
        noise = sigma * np.random.randn(G)
        x[:,k] = x[:,k-1] + d*dt + noise*np.sqrt(dt)

    class S:
        t = times
        y = x
    return S, gene_nodes, reporter_nodes


def simulate_dde(nodes, edges, node_params, params):
    gene_nodes     = [n for n in nodes if n.type=='gene']
    inducer_nodes  = [n for n in nodes if n.type=='inducer']
    reporter_nodes = [n for n in nodes if n.type=='reporter']

    G = len(gene_nodes)
    gene_index = {n.id:i for i,n in enumerate(gene_nodes)}

    act_g = {i:[] for i in range(G)}
    rep_g = {i:[] for i in range(G)}
    act_i = {i:[] for i in range(G)}
    rep_i = {i:[] for i in range(G)}
    on_times = {nid: node_params[nid]['on_time'].get()
                for nid in node_params if nid in [u.id for u in inducer_nodes]}

    for e in edges:
        if e.target.type!='gene': continue
        ti = gene_index[e.target.id]
        wvar = getattr(e,'weight',None)
        dvar = getattr(e,'delay',None)
        w = wvar.get() if wvar is not None else 1.0
        d = dvar.get() if dvar is not None else 0.0

        if e.source.type=='gene':
            si = gene_index.get(e.source.id)
            if si is None: continue
            if e.type=='activator': act_g[ti].append((si,d,w))
            else:                 rep_g[ti].append((si,d,w))
        else:
            if e.type=='activator': act_i[ti].append((e.source.id,d,w))
            else:                    rep_i[ti].append((e.source.id,d,w))

    basal = params['basal']
    t0, t1, steps = params['t0'], params['t1'], params['steps']
    times = np.linspace(t0, t1, steps)
    dt = (t1 - t0) / max(steps-1,1)

    x0 = np.array([node_params[n.id]['ic'].get() for n in gene_nodes])
    x = np.zeros((G, steps))
    x[:,0] = x0

    for k in range(1, steps):
        t = times[k-1]
        dx = np.zeros(G)
        for i,node in enumerate(gene_nodes):
            p = node_params[node.id]
            alpha_i = p['alpha'].get()
            hill_i  = p['n'].get()
            deg_i   = p['deg'].get()

            pa_g = 1.0
            for si,d,w in act_g[i]:
                td = t - d
                if td <= t0: xj = x0[si]
                else:
                    idx = (td - t0)/dt
                    lo = int(np.floor(idx))
                    hi = min(lo+1, k-1)
                    frac = idx - lo
                    xj = x[si,lo]*(1-frac) + x[si,hi]*frac
                pa_g *= (xj**hill_i/(1+xj**hill_i))**w
            pr_g = 1.0
            for si,d,w in rep_g[i]:
                td = t - d
                if td <= t0: xj = x0[si]
                else:
                    idx = (td - t0)/dt
                    lo = int(np.floor(idx))
                    hi = min(lo+1, k-1)
                    frac = idx - lo
                    xj = x[si,lo]*(1-frac) + x[si,hi]*frac
                pr_g *= (1.0/(1+xj**hill_i))**w

            pa_i = 1.0
            for jid,d,w in act_i[i]:
                if t >= on_times[jid] + d: pa_i *= 1.0**w
                else:                    pa_i *= 0.0
            pr_i = 1.0
            for jid,d,w in rep_i[i]:
                if t >= on_times[jid] + d: pr_i *= (1.0/(1+1.0**hill_i))**w
                else:                    pr_i *= (1.0/(1+0.0**hill_i))**w

            dx[i] = basal + alpha_i*pa_g*pa_i*pr_g*pr_i - deg_i*x[i,k-1]
        x[:,k] = x[:,k-1] + dx*dt

    class S: t = times; y = x
    return S, gene_nodes, reporter_nodes


def simulate_stochastic_dde(nodes, edges, node_params, params):
    S_det, gene_nodes, reporter_nodes = simulate_dde(nodes, edges, node_params, params)
    times = S_det.t
    x_det = S_det.y
    G = len(gene_nodes)
    sigma = params.get('sigma', 0.1)
    x = np.zeros_like(x_det)
    x[:,0] = x_det[:,0]
    dt = times[1] - times[0] if len(times)>1 else 0.0


    gene_index = {n.id:i for i,n in enumerate(gene_nodes)}
    act_g = {i:[] for i in range(G)}; rep_g = {i:[] for i in range(G)}
    act_i = {i:[] for i in range(G)}; rep_i = {i:[] for i in range(G)}
    on_times = {nid: node_params[nid]['on_time'].get() for nid in node_params if nid in [u.id for u in nodes if u.type=='inducer']}
    for e in edges:
        if e.target.type!='gene': continue
        ti = gene_index[e.target.id]
        w = getattr(e, 'weight', None).get() if getattr(e, 'weight', None) else 1.0
        d = getattr(e, 'delay', None).get() if getattr(e, 'delay', None) else 0.0
        if e.source.type=='gene':
            si = gene_index.get(e.source.id)
            if si is None: continue
            if e.type=='activator': act_g[ti].append((si,d,w))
            else:                 rep_g[ti].append((si,d,w))
        else:
            if e.type=='activator': act_i[ti].append((e.source.id,d,w))
            else:                    rep_i[ti].append((e.source.id,d,w))

    def drift_vector(xv, t):
        dx = np.zeros(G)
        for i,node in enumerate(gene_nodes):
            p = node_params[node.id]
            alpha_i = p['alpha'].get(); hill_i = p['n'].get(); deg_i = p['deg'].get()
            pa_g = 1.0
            for si,d,w in act_g[i]:
                td = t - d; idx = max(int((td - times[0])/dt),0)
                xj = xv[si] if td<=times[0] else x[si,idx]
                pa_g *= (xj**hill_i/(1+xj**hill_i))**w
            pr_g = 1.0
            for si,d,w in rep_g[i]:
                td = t - d; idx = max(int((td - times[0])/dt),0)
                xj = xv[si] if td<=times[0] else x[si,idx]
                pr_g *= (1.0/(1+xj**hill_i))**w
            pa_i = np.prod([1.0 if t>= on_times[jid]+d else 0.0 for (jid,d,w) in act_i[i]]) if act_i[i] else 1.0
            pr_i = np.prod([1.0/(1+(1.0 if t>= on_times[jid]+d else 0.0)**hill_i)**w for (jid,d,w) in rep_i[i]]) if rep_i[i] else 1.0
            dx[i] = params['basal'] + alpha_i*pa_g*pa_i*pr_g*pr_i - deg_i*xv[i]
        return dx

    for k in range(1, x.shape[1]):
        t_prev = times[k-1]
        xv_prev = x[:,k-1]
        drift_k = drift_vector(xv_prev, t_prev)
        noise = sigma * np.random.randn(G)
        x[:,k] = xv_prev + drift_k*dt + noise*np.sqrt(dt)

    class S2: t = times; y = x
    return S2, gene_nodes, reporter_nodes




def simulate_network(nodes, edges, node_params, params, ax=None):
    has_delay = any(getattr(e, 'delay', None) and e.delay.get() > 0 for e in edges)
    stype = params.get('sim_type', 'deterministic')

    if has_delay:
        if stype == 'stochastic':
            sol, gene_nodes, reporter_nodes = simulate_stochastic_dde(nodes, edges, node_params, params)
        else:
            sol, gene_nodes, reporter_nodes = simulate_dde(nodes, edges, node_params, params)
    else:
        if stype == 'stochastic':
            sol, gene_nodes, reporter_nodes = simulate_stochastic(nodes, edges, node_params, params)
        else:
            sol, gene_nodes, reporter_nodes = simulate_deterministic(nodes, edges, node_params, params)

    created = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6,4)); created = True

    color_map = {
        'GFP':'green','RFP':'red','BFP':'blue',
        'YFP':'yellow','CFP':'cyan','mCherry':'magenta'
    }
    use_rep = bool(reporter_nodes) and not params.get('simulate_all', False)
    if use_rep:
        idx = {n.id: i for i,n in enumerate(gene_nodes)}
        reps = {}
        for rep in reporter_nodes:
            typ = node_params[rep.id]['reporter'].get()
            reps.setdefault(typ, []).append(rep)
        for typ, rs in reps.items():
            sig = np.zeros_like(sol.t)
            for rep in rs:
                src = next((e.source.id for e in edges if e.target.id==rep.id and e.source.type=='gene'), None)
                if src in idx: sig += sol.y[idx[src]]
            ax.plot(sol.t, sig, label=typ, color=color_map.get(typ))
    else:
        for i,n in enumerate(gene_nodes): ax.plot(sol.t, sol.y[i], label=f'Gene {n.id}')

    ax.set_xlabel('Time'); ax.set_ylabel('Expression')
    title = ('Stochastic DDE' if has_delay and stype=='stochastic'
             else 'DDE' if has_delay else stype.capitalize())
    ax.set_title(f"Simulation ({title})")
    ax.legend(loc='upper right'); ax.grid(True)

    if created: plt.tight_layout(); plt.show()


def _rhs(t, x, nodes, edges, node_params, global_params, pre=None):
    if pre is None:
        pre = _prepare_network(nodes, edges)
    gene_nodes, idx, act_g, rep_g, act_i, rep_i = pre

    alpha = global_params["alpha"]
    n_h   = global_params["n"]
    basal = global_params["basal"]

    level = {
        nid: node_params.get(nid, {}).get("level", tk.DoubleVar(value=0.0)).get()
        for nid in list(node_params.keys())
    }

    dx = np.zeros_like(x)
    for i in range(len(gene_nodes)):
        pa_g = np.prod([(x[j] ** n_h / (1 + x[j] ** n_h)) ** w
                        for j, w in act_g[i]]) if act_g[i] else 1.0
        pr_g = np.prod([(1 / (1 + x[j] ** n_h)) ** w
                        for j, w in rep_g[i]]) if rep_g[i] else 1.0

        pa_i = np.prod([(level[nid]) ** w
                        for nid, w in act_i[i]]) if act_i[i] else 1.0
        pr_i = np.prod([(1 / (1 + level[nid] ** n_h)) ** w
                        for nid, w in rep_i[i]]) if rep_i[i] else 1.0

        dx[i] = basal + alpha * pa_g * pa_i * pr_g * pr_i - x[i]
    return dx


def simulate_step(nodes, edges, node_params, y0, t0, t1,
                  *, dt=1.0, sigma=0.0, sim_type="deterministic",
                  alpha_global=10.0, n_global=2.0, basal_global=0.1):

    import numpy as np

    global_p = dict(alpha=alpha_global,
                    n=n_global,
                    basal=basal_global)

    rhs = _rhs                      

    y = y0.copy()
    t = t0
    n_sub = max(1, int(round((t1 - t0) / dt)))
    h = (t1 - t0) / n_sub
    G = y.size

    if sim_type == "stochastic" and sigma > 0:
        sqrt2sig = (2 * sigma) ** 0.5
        rng = np.random.default_rng()

        for _ in range(n_sub):
            dy = rhs(t, y, nodes, edges, node_params, global_p)
            y += h * dy + sqrt2sig * np.sqrt(h) * rng.standard_normal(G)
            y = np.clip(y, 0.0, None)
            t += h
    else: 
        for _ in range(n_sub):
            dy = rhs(t, y, nodes, edges, node_params, global_p)
            y += h * dy
            y = np.clip(y, 0.0, None)
            t += h

    return y

