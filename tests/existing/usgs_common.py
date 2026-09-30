"""Shared set-up of the tests with USGS example models (existing-model mode, :OverlapWeights).

The example models are not georeferenced, so the Liard Raven model is laid over them: the coupled HRUs of subbasin 52 get
zones of the model's active top-layer cells (nearest of a few k-means seeds on the cell centres; weight 1.0 per cell) and
.rvh areas equal to their zones, so that Raven's recharge depths become the same depths on the model.
"""
import os, re, shutil
import numpy as np
import testconfig as tc


def active_cells(disfile):
    """(nrow, ncol, delr, delc, [(row, col) of active top-layer cells]) of a DIS file with CONSTANT delr/delc"""
    dis = open(disfile).read()
    nrow = int(re.search(r'NROW\s+(\d+)', dis, re.I).group(1)); ncol = int(re.search(r'NCOL\s+(\d+)', dis, re.I).group(1))
    delr = float(re.search(r'delr\s+CONSTANT\s+(\S+)', dis, re.I).group(1)); delc = float(re.search(r'delc\s+CONSTANT\s+(\S+)', dis, re.I).group(1))
    m = re.search(r'idomain\s*\n\s*(INTERNAL[^\n]*|CONSTANT\s+\S+)\n', dis, re.I)
    if m is None:
        idom = np.ones((nrow, ncol), int)
    elif m.group(1).upper().startswith('CONSTANT'):
        idom = np.full((nrow, ncol), int(float(m.group(1).split()[1])))
    else:
        vals = dis[m.end():].split()[:nrow * ncol]
        idom = np.array([int(float(v)) for v in vals]).reshape(nrow, ncol)
    return nrow, ncol, delr, delc, [(r, c) for r in range(nrow) for c in range(ncol) if idom[r, c] > 0]


def make_case(W, name, model_dir, disfile, rvg_extra, days, start='1985-10-01', rvc=None):
    """a Liard case coupled to the model in model_dir through zones of its active cells; returns the case folder"""
    nrow, ncol, delr, delc, act = active_cells(disfile)
    d = os.path.join(W, name)
    shutil.rmtree(d, ignore_errors=True)
    tc.copy_case(tc.LIARD, d)
    tc.sb52(d)                                   # groundwater under subbasin 52 only
    rvh = open(d + '/Liard.rvh').read()
    head, rest = rvh.split(':HRUs', 1); body, tail = rest.split(':EndHRUs', 1)
    lines = body.split('\n'); hdr = [l for l in lines if ':Attributes' in l][0].replace(',', ' ').split()[1:]
    ia, iar = hdr.index('AQUIFER_PROFILE') + 1, hdr.index('AREA') + 1
    coupled = [l.split(',')[0].strip() for l in lines if l.strip()[:1].isdigit() and l.split(',')[ia].strip() != '[NONE]']
    nz = min(len(coupled), len(act))
    xy = np.array([(c * delr + delr / 2, (nrow - r) * delc - delc / 2) for r, c in act])
    rng = np.random.default_rng(52); seeds = xy[rng.choice(len(xy), nz, replace=False)]
    for _ in range(10):
        lab = np.argmin(((xy[:, None, :] - seeds[None, :, :]) ** 2).sum(-1), axis=1)
        seeds = np.array([xy[lab == k].mean(0) if np.any(lab == k) else seeds[k] for k in range(nz)])
    area = {coupled[k]: (lab == k).sum() * delr * delc / 1e6 for k in range(nz)}
    out = []
    for l in lines:
        cc = [x.strip() for x in l.split(',')]
        if l.strip()[:1].isdigit() and cc[0] in coupled:
            if area.get(cc[0], 0) > 0:
                cc[iar] = '%.6f' % area[cc[0]]
            else:
                cc[ia] = '[NONE]'                # more HRUs than cells: the rest stay uncoupled
            l = ', '.join(cc)
        out.append(l)
    open(d + '/Liard.rvh', 'w').write(head + ':HRUs' + '\n'.join(out) + ':EndHRUs' + tail)
    w = ':OverlapWeights\n  :Attributes, HRU_ID, CELL, WEIGHT\n' + ''.join(
        '  %s %d 1.0\n' % (coupled[lab[i]], r * ncol + c + 1) for i, (r, c) in enumerate(act)) + ':EndOverlapWeights\n'
    open(d + '/Liard.rvg', 'w').write(':MF6Simulation  %s/mfsim.nam\n:MF6StartDate   1985-10-01\n%s\n%s' % (model_dir, rvg_extra, w))
    t = open(d + '/Liard.rvi').read()
    t = re.sub(r':Duration\s+\S+(\s+\S+)?', ':Duration %d' % days, t)
    t = re.sub(r':StartDate\s+\S+\s+\S+', ':StartDate %s 00:00:00' % start, t)
    open(d + '/Liard.rvi', 'w').write(t)
    if rvc:
        shutil.copy(rvc, d + '/Liard.rvc')
    return d


def audit(d, resid=1e-8):
    """the audit's identities for an existing-model run: [(name, value, limit)]"""
    import pandas as pd
    A = tc.audit.audit(d + '/out')
    b = pd.read_csv(d + '/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]'])
    chk = [('recharge identity', abs(A['recharge: (Raven - MF6 - dStore - toCHD)/Raven']), 1e-9),
           ('MODFLOW cumulative residual', abs(A['MF6 cumulative residual / gross flow']), resid),
           ('river bookkeeping', abs(A['river: applied - net MF6 - carry change (rel)']), 1e-8),
           ('unconverged steps', A['unconverged steps'], 0.5)]
    if 'seepage returned vs MF6 seepage (rel)' in A:
        chk.append(('seepage returned', abs(A['seepage returned vs MF6 seepage (rel)']), 1e-9))
    return A, b, chk
