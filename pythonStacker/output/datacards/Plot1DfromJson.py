import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
import matplotlib.lines as mlines
from matplotlib.legend import Legend
import json
import argparse
from scipy.signal import savgol_filter
import os

# --- Argument parsing ---
parser = argparse.ArgumentParser(description="Plot 1D WC limits from JSON files")
parser.add_argument("--exp", required=True, help="Expected limits JSON file")
parser.add_argument("--obs", required=True, help="Observed limits JSON file")
parser.add_argument("--WC", required=True, help="Name of the Wilson Coefficient to plot")
parser.add_argument("--xrange", nargs=2, type=float, default=[-2, 2],
                    help="X-axis range: min max")
parser.add_argument("--yrange", nargs=2, type=float, default=[0, 9.2],
                    help="Y-axis range: min max")
parser.add_argument("--smoothing", action='store_true', default=False, help="smooth the shape?")
args = parser.parse_args()

# --- Mapping WC to LaTeX labels ---
wc_label_map = {
    "ctt":   r"$c_{tt}$",
    "cQQ1":   r"$c_{QQ}^{(1)}$",
    "cQt1":   r"$c_{Qt}^{(1)}$",
    "cQt8":   r"$c_{Qt}^{(8)}$",
    "ctHRe": r"$c_{tH}^{Re}$",
    "ctHIm": r"$c_{tH}^{Im}$",
    # add more mappings here
}
x_label = wc_label_map.get(args.WC, args.WC)

# --- Load limits from JSON ---
def load_limits(jsonfile, wc_name):
    with open(jsonfile) as f:
        data = json.load(f)
    for block, limits in data.items():
        if limits.get("WC") == wc_name:
            fixed_arr = np.array(limits["fixed"])
            prof_arr = np.array(limits["profiled"])
            return (fixed_arr[:,0], fixed_arr[:,1],
                    prof_arr[:,0], prof_arr[:,1])
    raise KeyError(f"WC '{wc_name}' not found in {jsonfile}")

x_exp_fixed, y_exp_fixed, x_exp_prof, y_exp_prof = load_limits(args.exp, args.WC)
x_obs_fixed, y_obs_fixed, x_obs_prof, y_obs_prof = load_limits(args.obs, args.WC)

# — apply Savitzky–Golay smoothing to knock down tiny jagged bumps —
if args.smoothing :
    y_exp_fixed  = savgol_filter(y_exp_fixed, 5, 3)
    y_exp_prof   = savgol_filter(y_exp_prof,   5, 3)
    y_obs_fixed  = savgol_filter(y_obs_fixed,  5, 3)
    y_obs_prof   = savgol_filter(y_obs_prof,   5, 3)

# --- Threshold and CL levels ---
threshold_z = 7.5
levels = np.array((1.0, 4.0))

# --- Mask & interpolate values above threshold ---
def mask_and_transform(x, z):
    z = z - np.min(z)
    x_clipped, z_clipped = [], []
    for i in range(len(x) - 1):
        x0, x1 = x[i], x[i+1]
        z0, z1 = z[i], z[i+1]
        if z0 <= threshold_z and z1 <= threshold_z:
            x_clipped.append(x0); z_clipped.append(z0)
        elif z0 <= threshold_z:
            x_clipped.append(x0); z_clipped.append(z0)
            x_i = x0 + (threshold_z - z0) / (z1 - z0) * (x1 - x0)
            x_clipped.append(x_i); z_clipped.append(threshold_z)
        elif z1 <= threshold_z:
            x_i = x0 + (threshold_z - z0) / (z1 - z0) * (x1 - x0)
            x_clipped.append(x_i); z_clipped.append(threshold_z)
        else:
            continue
#            x_clipped.append(np.nan); z_clipped.append(np.nan)
    return np.array(x_clipped), np.array(z_clipped)

def maybe_extrapolate(x, y, x_min, x_max, n=5):
    # if data already reaches or exceeds both ends, do nothing
    if x.min() < x_min and x.max() > x_max:
        return x, y

    #print("we're extrapolating something :)")
    #print(x.min(), "<", x_min)
    #print(x.max(), ">", x_max)
    # otherwise perform linear extrapolation on each end as needed
    # sort by x
    #order = np.argsort(x)
    #x, y = x[order], y[order]

    # fit & extrapolate lower end
    if x.min() > x_min:
        lo_idx = np.arange(min(n, len(x)))
        coef_lo = np.polyfit(x[lo_idx], y[lo_idx], 1)
        y_lo = np.poly1d(coef_lo)(x_min)
        x = np.insert(x, 0, x_min)
        y = np.insert(y, 0, y_lo)
        #print(x,y)

    # fit & extrapolate upper end
    if x.max() < x_max:
        hi_idx = np.arange(len(x)-min(n, len(x)), len(x))
        coef_hi = np.polyfit(x[hi_idx], y[hi_idx], 1)
        y_hi = np.poly1d(coef_hi)(x_max)
        x = np.append(x, x_max)
        y = np.append(y, y_hi)
        #print(x,y)

    return x, y

def find_intersections(x, y, level):
    xs = []
    # 1) any exact matches
    for xi, yi in zip(x, y):
        if np.isclose(yi, level):
            xs.append(xi)

    # 2) zero-crossings between points
    signs = y - level
    idx = np.where(signs[:-1] * signs[1:] < 0)[0]
    for i in idx:
        x0, x1 = x[i], x[i+1]
        y0, y1 = y[i], y[i+1]
        # linear interpolation
        xi = x0 + (level - y0)*(x1 - x0)/(y1 - y0)
        xs.append(xi)

    # sort and return all found crossings
    return sorted(xs)

#    # find segments straddling the horizontal line y=level
#
#    signs = y - level
#    idx = np.where(signs[:-1] * signs[1:] < 0)[0]
#    xs = []
#    for i in idx:
#        x0, x1 = x[i],   x[i+1]
#        y0, y1 = y[i],   y[i+1]
#        xi = x0 + (level-y0)*(x1-x0)/(y1-y0)
#        xs.append(xi)
#    return xs

# --- Plotting ---
fig, ax = plt.subplots(figsize=(8, 8))
# CL lines
ax.axhline(levels[0], color='grey', linestyle='dotted')
ax.axhline(levels[1], color='grey', linestyle='dotted')

#fill the x-range
x_obs_p, y_obs_p = maybe_extrapolate(x_obs_prof, y_obs_prof,args.xrange[0], args.xrange[1])
x_obs_p, y_obs_p = mask_and_transform(x_obs_p, y_obs_p)
x_exp_p, y_exp_p = maybe_extrapolate(x_exp_prof, y_exp_prof,args.xrange[0], args.xrange[1])
x_exp_p, y_exp_p = mask_and_transform(x_exp_p, y_exp_p)
x_obs_f, y_obs_f = maybe_extrapolate(x_obs_fixed, y_obs_fixed,args.xrange[0], args.xrange[1])
x_obs_f, y_obs_f = mask_and_transform(x_obs_f, y_obs_f)
x_exp_f, y_exp_f = maybe_extrapolate(x_exp_fixed, y_exp_fixed,args.xrange[0], args.xrange[1])
x_exp_f, y_exp_f = mask_and_transform(x_exp_f, y_exp_f)


# Expected: dashed (fixed black, profiled light green)
ax.plot(x_exp_f, y_exp_f, color='black', linestyle='dashed', label='Expected (frozen)')
ax.plot(x_exp_p,  y_exp_p,  color='limegreen', linestyle='dashed', label='Expected (profiled)')

# Observed: solid (fixed black, profiled light green)
ax.plot(x_obs_f, y_obs_f, color='black', linestyle='solid', label='Observed (frozen)')
ax.plot(x_obs_p,  y_obs_p,  color='limegreen', linestyle='solid', label='Observed (profiled)')

# Axis labels and ranges
ax.set_xlabel(x_label, fontsize=28)
ax.set_xlim(args.xrange)
ax.set_ylabel(r"$-2\,\Delta\,\mathrm{ln}\,\mathit{L}$", fontsize=28, labelpad=12)
ax.set_ylim(args.yrange)
ax.tick_params(axis='both', labelsize=24)

# Legend
legend = ax.legend(fontsize=18, ncol=2, loc='upper center', frameon=True)
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('none')

# CMS label
hep.cms.label(ax=ax, data=True, label='$~$Preliminary',
              rlabel=r'138$\,$fb$^{-1}$ (13$\,$TeV)', fontsize=26)

plt.tight_layout()
plt.savefig(f"/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v36/1D_limits/{args.WC}_1dlimits.png", bbox_inches='tight')
plt.savefig(f"/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v36/1D_limits/{args.WC}_1dlimits.pdf", bbox_inches='tight')


# --- compute intersections and write LaTeX table ---
curves = {
    'Exp fixed':     (x_exp_f, y_exp_f),
    'Exp profiled':  (x_exp_p,    y_exp_p),
    'Obs fixed':     (x_obs_f,   y_obs_f),
    'Obs profiled':  (x_obs_p,    y_obs_p),
}

# optionally extrapolate here if you need exact xrange coverage
for k, (xv, yv) in curves.items():
    curves[k] = maybe_extrapolate(xv, yv, *args.xrange)

lines = []
for name, (xv, yv) in curves.items():
    i1 = find_intersections(xv, yv, levels[0])
    i4 = find_intersections(xv, yv, levels[1])
    # assume two intersections each
    lines.append(f"{name} & {i1[0]:.3f} & {i1[-1]:.3f} & {i4[0]:.3f} & {i4[-1]:.3f}\\\\\n")

tex_path = f"./v36_IM/json_files/{args.WC}_1dlimits_table.tex"
# --- write Markdown table ---
#md_path = f"{args.WC}_1dlimits_table.md"
header = "| Curve         | best-fit ($y=0$) | $x$ @ 1       | $x$ @ 4       |\n"
sep    = "|---------------|------------------|---------------|---------------|\n"
with open(tex_path, "w") as mf:
    mf.write(header)
    mf.write(sep)
    for name, (xv, yv) in curves.items():
        # intersections
        i1 = find_intersections(xv, yv, levels[0])
        i4 = find_intersections(xv, yv, levels[1])
        # join *all* intersections as comma-separated lists
        i1_str = ", ".join(f"{xi:.3f}" for xi in i1)
        i4_str = ", ".join(f"{xi:.3f}" for xi in i4)
        xb = xv[np.argmin(yv)]
        mf.write(f"| {name:<13} | {xb:>16.3f} | {i1_str:<13} | {i4_str:<13} |\n")
#        x1, x2 = i1[0], i1[-1]
#        x3, x4 = i4[0], i4[-1]
#        # best-fit = x at minimum y
#        mf.write(f"| {name:<13} | {xb:>16.3f} | {x1:>9.3f} | {x2:>9.3f} | {x3:>9.3f} | {x4:>9.3f} |\n")

print(f"Wrote table to {tex_path}")
