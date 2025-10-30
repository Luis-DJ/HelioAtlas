#!/usr/bin/env python3
"""
HelioAtlas — 3I/ATLAS interactive plots (presentation-ready)

Features
- Robust ephemerides: HORIZONS (validated) with automatic MPC fallback (configurable)
- UTC timeline, nearest alignment (±12 h)
- Figure 1: top-down ecliptic view with monthly labels (Sep→Feb), optional weekly dots
  Sun (gold), Earth (blue), Mars (red), Comet (custom marker, purple)
  Compact legend (bottom-right), readout box (top-right): date, distances, λ/β, phase, elongation
- Figure 2: distances vs time with monthly/weekly vertical guides
  Compact legend (bottom-right) + mirrored readout box
- Clean timeline controls:
  Slider with NO label and NO numeric readout on the right
  Start/end dates (dd-mmm-yyyy) pinned under slider ends
  Centered “Days since …” counter under the slider (DISCOVERY date optional)
"""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.path import Path
from datetime import datetime, timezone
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB
from astropy.time import Time

# -------------------- CONFIG --------------------
FORCE_MPC = True        # False: try HORIZONS first (auto-fallback); True: force MPC comet track
START = '2025-07-01'
STOP  = '2026-03-31'
STEP  = '1d'
MARK  = '2025-12-19'
PERI  = '2025-10-29 00:00Z'

# Month labels to place along the comet trajectory (kept if within span)
EXPLICIT_MONTH_LABELS = ["2025-09-01","2025-10-01","2025-11-01","2025-12-01","2026-01-01","2026-02-01"]
SHOW_WEEKLY_DOTS = True   # set False to hide weekly dots on the orbit figure

# Presentation
DATE_FMT   = '%d-%b-%Y'   # dd-mmm-yyyy
DISCOVERY  = None         # e.g., '2025-01-09' (UTC). If None, counter uses START.

# Colors
SUN_COLOR   = 'gold'
EARTH_COLOR = 'tab:blue'
MARS_COLOR  = 'tab:red'
COMET_COLOR = 'purple'
MONTH_DOT   = 'gold'
WEEK_DOT    = 'lightgray'
GUIDE_COLOR = 'lightgray'

# MPC fallback orbital elements (update if newer)
Q_AU = 1.3568482
E    = 6.1423568
INC  = np.radians(175.11282)
OMEG = np.radians(322.15521)
ARGP = np.radians(128.00073)
T_PERI = datetime(2025, 10, 29, 11, 29, 0, tzinfo=timezone.utc)

CANDIDATE_KEYS = [
    "DES=I003","I003","3I/2025 A3 (ATLAS)","3I/2025 A3","3I/ATLAS","3I ATLAS","A/2025 A3","3I",
]

# -------------------- Utilities --------------------
def ensure_utc(obj):
    s = pd.to_datetime(obj, errors='raise')
    if isinstance(s, pd.DatetimeIndex):
        s = pd.Series(s)
    try:
        tzinfo = getattr(s.dt, "tz", None)
    except Exception:
        tzinfo = None
    return s.dt.tz_localize("UTC") if tzinfo is None else s.dt.tz_convert("UTC")

def parse_horizons_datetime_str(series):
    s = (series.astype(str)
         .str.replace('A.D. ', '', regex=False)
         .str.replace(' TDB', '', regex=False)
         .str.replace(' UT',  '', regex=False)
         .str.replace(' UTC', '', regex=False))
    try:
        return pd.to_datetime(s, format='%Y-%b-%d %H:%M:%S.%f', utc=True)
    except Exception:
        return pd.to_datetime(s, utc=True)

def horizons_times_to_utc(tbl):
    if 'datetime_str' in tbl.colnames:
        try:
            return parse_horizons_datetime_str(tbl['datetime_str'].astype(str))
        except Exception:
            pass
    for cand in ('JD','datetime_jd'):
        if cand in tbl.colnames:
            t = Time(np.array(tbl[cand], dtype=float), format='jd', scale='tdb')
            return pd.to_datetime(t.utc.datetime)
    raise RuntimeError("No time column in HORIZONS table")

def angle_deg(u, v):
    u = np.asarray(u, dtype=float); v = np.asarray(v, dtype=float)
    nu = np.linalg.norm(u); nv = np.linalg.norm(v)
    if nu == 0 or nv == 0: return np.nan
    cosang = np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0)
    return float(np.degrees(np.arccos(cosang)))

def lonlat_from_xyz(x, y, z):
    lam = math.degrees(math.atan2(y, x)) % 360.0
    rho = math.hypot(x, y)
    beta = math.degrees(math.atan2(z, rho))
    return lam, beta

# -------------------- HORIZONS wrappers --------------------
def _vectors_majorbody(naif_id, location, refplane):
    obj = Horizons(id=str(naif_id), id_type='majorbody', location=location,
                   epochs={'start': START, 'stop': STOP, 'step': STEP})
    return obj.vectors(refplane=refplane) if refplane else obj.vectors()

def _vectors_try(key, location, refplane):
    obj = Horizons(id=key, location=location,
                   epochs={'start': START, 'stop': STOP, 'step': STEP})
    return obj.vectors(refplane=refplane) if refplane else obj.vectors()

def _vectors_try_designation(key, location, refplane):
    obj = Horizons(id=key, id_type='designation', location=location,
                   epochs={'start': START, 'stop': STOP, 'step': STEP})
    return obj.vectors(refplane=refplane) if refplane else obj.vectors()

def _vectors_checked(key, location, refplane):
    """Query and verify this is 3I/ATLAS and perihelion distance is plausible."""
    for attempt in (_vectors_try, _vectors_try_designation):
        try:
            tbl = attempt(key, location, refplane)
            tname = str(tbl.meta.get('targetname', ''))
            if not (('3I' in tname) and ('ATLAS' in tname.upper())):
                raise ValueError(f"Unexpected target: {tname}")
            t_utc = horizons_times_to_utc(tbl)
            r = np.sqrt(np.array(tbl['x'], float)**2 +
                        np.array(tbl['y'], float)**2 +
                        np.array(tbl['z'], float)**2)
            peri = pd.Timestamp(PERI, tz='UTC')
            mask = np.abs((pd.to_datetime(t_utc, utc=True) - peri)
                          .astype('timedelta64[s]').astype(float)) < 5*86400
            if mask.any() and np.nanmin(r[mask]) > 3.0:
                raise ValueError("Inconsistent perihelion distance (>3 AU).")
            return tbl
        except Exception:
            continue
    raise RuntimeError(f"HORIZONS failed or wrong body for key '{key}'")

def _probe_key(key):
    try:
        _vectors_checked(key, '@sun', 'ecliptic'); return True
    except Exception:
        return False

def resolve_horizons_key():
    for key in CANDIDATE_KEYS:
        if _probe_key(key):
            print(f"[resolver] Using candidate key: {key}")
            return key, 'candidate'
    try:
        rec = SBDB.query("3I/ATLAS")
        des  = rec.get('object', {}).get('des')
        full = rec.get('object', {}).get('full_name')
        if des and _probe_key(f"DES={des}"):
            key = f"DES={des}"; print(f"[resolver] Using SBDB DES key: {key}")
            return key, 'sbdb-des'
        if full and _probe_key(full):
            print(f"[resolver] Using SBDB full_name key: {full}")
            return full, 'sbdb-name'
    except Exception:
        pass
    raise RuntimeError("No valid HORIZONS key for 3I/ATLAS")

# -------------------- Fetch helpers --------------------
def fetch_heliocentric_xyz(target_key):
    if str(target_key).isdigit():  # planets/major bodies
        tbl = _vectors_majorbody(target_key, '@sun', 'ecliptic')
    else:                          # comet (validated)
        tbl = _vectors_checked(target_key, '@sun', 'ecliptic')
    df = tbl.to_pandas()
    df['t'] = ensure_utc(horizons_times_to_utc(tbl))
    return df[['t','x','y','z']].sort_values('t').reset_index(drop=True)

# -------------------- MPC propagation (comet) --------------------
def rot_z(a): c,s=np.cos(a),np.sin(a); return np.array([[c,-s,0],[s,c,0],[0,0,1]])
def rot_x(a): c,s=np.cos(a),np.sin(a); return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def comet_state_hyperbola(q,e,dt_days,Omega,inc,omega):
    k=0.01720209895; mu=k**2; a=q/(e-1.0)
    def tof(H): return np.sqrt(a**3/mu)*(e*np.sinh(H)-H)
    H = dt_days*np.sqrt(mu/a**3)/max(e-1.0,1e-9)
    for _ in range(80):
        f = tof(H)-dt_days; fp = np.sqrt(a**3/mu)*(e*np.cosh(H)-1.0)
        dH = -f/fp; H += dH
        if abs(dH) < 1e-12: break
    rmag = a*(e*np.cosh(H)-1.0)
    tanv2 = np.sqrt((e+1)/(e-1))*np.tanh(H/2.0)
    nu = 2*np.arctan(tanv2)
    rp = rmag*np.array([np.cos(nu), np.sin(nu), 0.0])
    R = rot_z(OMEG) @ rot_x(INC) @ rot_z(ARGP)
    return R @ rp, rmag

def mpc_comet_track():
    days = pd.date_range(
        start=pd.Timestamp(START).tz_localize("UTC"),
        end=pd.Timestamp(STOP).tz_localize("UTC"),
        freq="1D"
    )
    xs, ys, zs = [], [], []
    for t in days:
        dt = (t.to_pydatetime() - T_PERI).total_seconds() / 86400.0
        r, _ = comet_state_hyperbola(Q_AU, E, dt, OMEG, INC, ARGP)
        xs.append(r[0]); ys.append(r[1]); zs.append(r[2])
    print(f"[info] MPC fallback generated {len(days)} points "
          f"from {days[0].date()} to {days[-1].date()}")
    return pd.DataFrame({'t': days, 'x': xs, 'y': ys, 'z': zs}).sort_values('t').reset_index(drop=True)

# -------------------- Alignment helper --------------------
def asof_align(base: pd.DataFrame, df: pd.DataFrame, cols, suffix):
    left = base[['t']].copy()
    right = df[['t'] + list(cols)].copy().sort_values('t')
    out = pd.merge_asof(left.sort_values('t'), right, on='t',
                        direction='nearest', tolerance=pd.Timedelta('12H'))
    out = out.add_suffix(suffix)
    out.rename(columns={f't{suffix}':'t'}, inplace=True)
    return out

# -------------------- Comet icon --------------------
def comet_marker_path():
    verts = np.array([
        [ 0.60,  0.00],  # nose
        [ 0.10,  0.18],
        [-0.55,  0.00],  # tail
        [ 0.10, -0.18],
        [ 0.60,  0.00],
    ])
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3]
    return Path(verts, codes)

# -------------------- Dataset builder --------------------
def build_dataset():
    base = pd.DataFrame({'t': pd.date_range(pd.Timestamp(START).tz_localize('UTC'),
                                            pd.Timestamp(STOP).tz_localize('UTC'),
                                            freq='1D')})
    earth_hc = fetch_heliocentric_xyz('399')
    mars_hc  = fetch_heliocentric_xyz('499')
    print(f"[earth] rows={len(earth_hc)} span={earth_hc['t'].min()} → {earth_hc['t'].max()}")
    print(f"[mars ] rows={len(mars_hc)}  span={mars_hc['t'].min()} → {mars_hc['t'].max()}")

    if not FORCE_MPC:
        try:
            key, how = resolve_horizons_key()
            print(f"[ok] HORIZONS key via {how}: {key}")
            comet_hc = fetch_heliocentric_xyz(key)
            mode = f"HORIZONS ({how})"
        except Exception as err:
            print(f"[warn] HORIZONS key failed ({err}). Using MPC fallback.")
            comet_hc = mpc_comet_track()
            mode = "MPC fallback (auto); Earth/Mars from HORIZONS"
    else:
        print("[info] FORCING MPC fallback for comet geometry.")
        comet_hc = mpc_comet_track()
        mode = "MPC fallback (forced); Earth/Mars from HORIZONS"

    if comet_hc.empty:
        raise RuntimeError("Comet track returned 0 rows.")
    print(f"[comet] rows={len(comet_hc)} span={comet_hc['t'].min()} → {comet_hc['t'].max()}")

    base = base.sort_values('t').reset_index(drop=True)
    aC = asof_align(base, comet_hc, cols=['x','y','z'], suffix='_C')
    aE = asof_align(base, earth_hc, cols=['x','y','z'], suffix='_E')
    aM = asof_align(base, mars_hc,  cols=['x','y','z'], suffix='_M')
    df = aC.merge(aE, on='t').merge(aM, on='t').dropna().reset_index(drop=True)

    df['sun_range_AU']   = np.sqrt(df['x_C']**2 + df['y_C']**2 + df['z_C']**2)
    df['earth_range_AU'] = np.sqrt((df['x_C']-df['x_E'])**2 + (df['y_C']-df['y_E'])**2 + (df['z_C']-df['z_E'])**2)
    df['mars_range_AU']  = np.sqrt((df['x_C']-df['x_M'])**2 + (df['y_C']-df['y_M'])**2 + (df['z_C']-df['z_M'])**2)

    # Lon/lat + phase & elongation
    lam_list, bet_list, phase_list, elong_list = [], [], [], []
    for xC, yC, zC, xE, yE, zE in zip(df['x_C'], df['y_C'], df['z_C'],
                                       df['x_E'], df['y_E'], df['z_E']):
        lam, bet = lonlat_from_xyz(xC, yC, zC)
        rC = np.array([xC, yC, zC])    # Sun->Comet
        rE = np.array([xE, yE, zE])    # Sun->Earth
        CE = rE - rC                   # Comet->Earth
        CS = -rC                       # Comet->Sun
        ES = -rE                       # Earth->Sun
        EC = rC - rE                   # Earth->Comet
        phase = angle_deg(CS, CE)      # Sun–Comet–Earth
        elong = angle_deg(ES, EC)      # Sun–Earth–Comet
        lam_list.append(lam); bet_list.append(bet)
        phase_list.append(phase); elong_list.append(elong)

    df['lambda_deg'] = lam_list
    df['beta_deg']   = bet_list
    df['phase_deg']  = phase_list
    df['elong_deg']  = elong_list

    print(f"[data] aligned rows={len(df)} span={df['t'].min()} → {df['t'].max()}")
    return df, mode

# -------------------- Plot + slider --------------------
def main():
    df, mode = build_dataset()
    df['t'] = ensure_utc(df['t'])
    mark_dt = pd.Timestamp(MARK, tz='UTC')
    peri_dt = pd.Timestamp(PERI, tz='UTC')

    diffs = (df['t'] - mark_dt).dt.total_seconds().abs()
    i0 = int(diffs.idxmin()) if diffs.notna().any() else len(df)//2

    t0, t1 = df['t'].min(), df['t'].max()
    desired = pd.to_datetime(EXPLICIT_MONTH_LABELS, utc=True)
    month_marks = [d for d in desired if (d >= t0 and d <= t1)]
    weekly_marks = pd.date_range(max(t0, pd.Timestamp("2025-09-01", tz='UTC')), t1, freq='W') if SHOW_WEEKLY_DOTS else []

    # --- Figure 1 ---
    theta = np.linspace(0, 2*np.pi, 1000)
    earth_orbit = np.vstack((np.cos(theta), np.sin(theta)))
    mars_orbit  = np.vstack((1.524*np.cos(theta), 1.524*np.sin(theta)))

    fig1, ax1 = plt.subplots(figsize=(9.5, 9.3))
    # Extra bottom space for slider + labels
    plt.subplots_adjust(left=0.06, right=0.985, top=0.94, bottom=0.20)

    ax1.plot(earth_orbit[0], earth_orbit[1], '--', color='0.75', label='Earth orbit (1 AU)')
    ax1.plot(mars_orbit[0],  mars_orbit[1],  '--', color='0.85', label='Mars orbit (1.524 AU)')
    ax1.scatter([0],[0], s=90, color=SUN_COLOR, edgecolor='k', label='Sun')
    ax1.plot(df['x_C'], df['y_C'], color=COMET_COLOR, linewidth=2.2, label=f'3I/ATLAS trajectory — {mode}')

    comet_mk = comet_marker_path()
    E = ax1.scatter([df['x_E'].iloc[i0]],[df['y_E'].iloc[i0]], s=60, color=EARTH_COLOR, edgecolor='k', label='Earth')
    M = ax1.scatter([df['x_M'].iloc[i0]],[df['y_M'].iloc[i0]], s=60, color=MARS_COLOR,  edgecolor='k', label='Mars')
    C = ax1.scatter([df['x_C'].iloc[i0]],[df['y_C'].iloc[i0]], s=150, color=COMET_COLOR, marker=comet_mk,
                    edgecolor='k', label='3I/ATLAS')

    lineE, = ax1.plot([df['x_E'].iloc[i0], df['x_C'].iloc[i0]],
                      [df['y_E'].iloc[i0], df['y_C'].iloc[i0]],
                      color=EARTH_COLOR, alpha=0.6)
    lineM,  = ax1.plot([df['x_M'].iloc[i0], df['x_C'].iloc[i0]],
                      [df['y_M'].iloc[i0], df['y_C'].iloc[i0]],
                      color=MARS_COLOR, alpha=0.6)

    # Month labels
    for j, d in enumerate(month_marks):
        i = int(np.argmin(np.abs((df['t'] - d).dt.total_seconds())))
        ax1.scatter(df['x_C'].iloc[i], df['y_C'].iloc[i], s=55, color=MONTH_DOT, edgecolor='k', zorder=6)
        yoff = 0.05 if (j % 2 == 0) else -0.07
        ax1.text(df['x_C'].iloc[i], df['y_C'].iloc[i] + yoff,
                 d.strftime('%b\n%Y'), fontsize=8, ha='center', va='center',
                 color='black', zorder=7,
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85, lw=0))

    # Weekly dots
    for d in weekly_marks:
        i = int(np.argmin(np.abs((df['t'] - d).dt.total_seconds())))
        ax1.scatter(df['x_C'].iloc[i], df['y_C'].iloc[i], s=10, color=WEEK_DOT, zorder=4)

    ax1.set_aspect('equal'); ax1.set_xlim(-3, 3); ax1.set_ylim(-3, 3)
    ax1.set_xlabel('x (AU)'); ax1.set_ylabel('y (AU)')
    ax1.set_title('3I/ATLAS — Top-down (ecliptic) with Month Labels')
    ax1.grid(True, ls=':', alpha=0.5)
    ax1.legend(loc='lower right', frameon=True, facecolor='white',
               fontsize=8, borderpad=0.3, labelspacing=0.3, handlelength=1.2)

    # Readout (top-right)
    def make_readout_text(i):
        d  = df['t'].iloc[i].strftime(DATE_FMT)
        sc = df['sun_range_AU'].iloc[i]
        ec = df['earth_range_AU'].iloc[i]
        mc = df['mars_range_AU'].iloc[i]
        lam = df['lambda_deg'].iloc[i]
        bet = df['beta_deg'].iloc[i]
        pha = df['phase_deg'].iloc[i]
        elo = df['elong_deg'].iloc[i]
        return (
            f"Date: {d}\n"
            f"Sun–comet:   {sc:5.2f} AU\n"
            f"Earth–comet: {ec:5.2f} AU\n"
            f"Mars–comet:  {mc:5.2f} AU\n"
            f"λ, β:        {lam:6.2f}°, {bet:6.2f}°\n"
            f"Phase (S–C–E): {pha:5.1f}°   Elong (S–E–C): {elo:5.1f}°"
        )
    tbox1 = ax1.text(0.98, 0.98, make_readout_text(i0),
                     transform=ax1.transAxes, va='top', ha='right',
                     fontsize=9, family='monospace',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, lw=0.5))

    # --- Slider area (clean: no label/right value) + dates + centered counter ---
    ax_slider = plt.axes([0.12, 0.08, 0.76, 0.03])  # lower so x-label is visible
    slider = Slider(ax=ax_slider, label='', valmin=0, valmax=len(df)-1, valinit=i0, valfmt='%0.0f')
    slider.valtext.set_visible(False)  # hide numeric readout on the right

    # Dates pinned to slider ends
    ax_slider.text(0.0, -0.60, pd.Timestamp(START).strftime(DATE_FMT),
                   transform=ax_slider.transAxes, ha='left', va='top', fontsize=9)
    ax_slider.text(1.0, -0.60, pd.Timestamp(STOP).strftime(DATE_FMT),
                   transform=ax_slider.transAxes, ha='right', va='top', fontsize=9)

    # Centered day counter underneath
    ax_counter = plt.axes([0.12, 0.045, 0.76, 0.025]); ax_counter.axis('off')
    ref_dt = (pd.Timestamp(DISCOVERY, tz='UTC') if DISCOVERY
              else pd.Timestamp(START).tz_localize('UTC'))
    def day_count_text(idx: int) -> str:
        days = int((df['t'].iloc[idx] - ref_dt).total_seconds() // 86400)
        label = 'Days since discovery' if DISCOVERY else f"Days since {pd.Timestamp(START).strftime(DATE_FMT)}"
        return f"{label}: {days:+d}"
    counter_txt = ax_counter.text(0.5, 0.5, day_count_text(i0),
                                  ha='center', va='center', fontsize=9)

    # --- Figure 2 ---
    fig2, ax2 = plt.subplots(figsize=(11.3, 6.2))
    plt.subplots_adjust(left=0.06, right=0.985, top=0.94, bottom=0.11)

    days = (df['t'] - peri_dt).dt.total_seconds()/86400.0
    ax2.plot(days, df['sun_range_AU'],   label='Sun–comet (AU)',   color=COMET_COLOR, lw=2.0)
    ax2.plot(days, df['earth_range_AU'], label='Earth–comet (AU)', color=EARTH_COLOR, lw=1.8)
    ax2.plot(days, df['mars_range_AU'],  label='Mars–comet (AU)',  color=MARS_COLOR,  lw=1.8)

    for d in weekly_marks:
        x = (d - peri_dt).total_seconds()/86400.0
        ax2.axvline(x, color=GUIDE_COLOR, lw=0.6, alpha=0.35, zorder=0)
    for d in month_marks:
        x = (d - peri_dt).total_seconds()/86400.0
        ax2.axvline(x, color='0.6', lw=1.2, alpha=0.6, zorder=0)
        ax2.text(x, ax2.get_ylim()[0], d.strftime('%b'), rotation=90,
                 va='bottom', ha='center', fontsize=8, color='0.3', alpha=0.9)

    vline = ax2.axvline((df['t'].iloc[i0] - peri_dt).total_seconds()/86400.0, ls=':', color='k')
    ax2.set_xlabel('Days from perihelion (2025-10-29)')
    ax2.set_ylabel('Distance (AU)')
    ax2.set_title('3I/ATLAS — Distances vs Time (with monthly/weekly guides)')
    ax2.grid(True, ls=':', alpha=0.5)
    ax2.legend(loc='lower right', frameon=True, facecolor='white',
               fontsize=8, borderpad=0.3, labelspacing=0.3, handlelength=1.2)

    # Mirror readout in Fig 2
    tbox2 = ax2.text(0.02, 0.98, make_readout_text(i0),
                     transform=ax2.transAxes, va='top', ha='left',
                     fontsize=9, family='monospace',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, lw=0.5))

    # --- Update + Reset ---
    ax_reset = plt.axes([0.90, 0.08, 0.07, 0.03]); btn = Button(ax_reset, 'Reset')

    def update(_):
        i = int(slider.val)
        # Move bodies
        E.set_offsets([[df['x_E'].iloc[i], df['y_E'].iloc[i]]])
        M.set_offsets([[df['x_M'].iloc[i], df['y_M'].iloc[i]]])
        C.set_offsets([[df['x_C'].iloc[i], df['y_C'].iloc[i]]])
        # Update connectors
        lineE.set_data([df['x_E'].iloc[i], df['x_C'].iloc[i]],
                       [df['y_E'].iloc[i], df['y_C'].iloc[i]])
        lineM.set_data([df['x_M'].iloc[i], df['x_C'].iloc[i]],
                       [df['y_M'].iloc[i], df['y_C'].iloc[i]])
        # Vertical line on Fig 2
        xcur = (df['t'].iloc[i] - peri_dt).total_seconds()/86400.0
        vline.set_xdata([xcur, xcur])
        # Readouts + counter
        txt = make_readout_text(i)
        tbox1.set_text(txt); tbox2.set_text(txt)
        counter_txt.set_text(day_count_text(i))
        fig1.canvas.draw_idle(); fig2.canvas.draw_idle()

    slider.on_changed(update)
    btn.on_clicked(lambda _: slider.reset())
    plt.show()

# -------------------- Entrypoint --------------------
if __name__ == '__main__':
    main()
