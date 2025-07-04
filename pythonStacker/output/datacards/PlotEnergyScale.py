#!/usr/bin/env python3
import ROOT
import math
import math
import numpy as np
import matplotlib.pyplot as plt


# ——— 1) YOUR INPUTS AS A DICTIONARY ——————————————————————————————
# Keys are the Wilson‐coefficient names (these become your x‑axis labels);
# values are the 95% CL bounds on |c_i|/Λ in TeV⁻¹.
wc_bounds = {
    r"$C_{tt}$": 1.30,
    r"$C_{QQ}^{(1)}$": 1.26,
    r"$C_{Qt}^{(1)}$":  1.83,
    r"$C_{Qt}^{(8)}$":  4.15,
    r"$C_{tH}^{Re}$":  8.06,
    r"$C_{tH}^{Im}$":  19.09,
}

# ┌────────────────────────────────────────────────────────────────────────────┐
# │ 2) BENCHMARK c_i VALUES                                                    │
# │    • Now at 1.0, 10.0, and (4π)^2                                            │
# └────────────────────────────────────────────────────────────────────────────┘
ci_map = {
    r"$1.0$":       1.0,
    r"$10.0$":     10.0,
    r"$(4\pi)^2$": (4*math.pi)**2,
}

# ┌────────────────────────────────────────────────────────────────────────────┐
# │ 2a) CUSTOM COLOR SELECTION                                                  │
# │    • Three matte green shades: dark, medium, light                          │
# └────────────────────────────────────────────────────────────────────────────┘
colors = {
    r"$1.0$":       "#1b5e20",  # matte dark green
    r"$10.0$":     "#388e3c",  # matte medium green
    r"$(4\pi)^2$": "#81c784",  # matte light green
}

# ┌────────────────────────────────────────────────────────────────────────────┐
# │ 3) COMPUTE Λ = |c_i| / bound FOR EACH CI                                   │
# └────────────────────────────────────────────────────────────────────────────┘
wc_names    = list(wc_bounds.keys())
bounds_arr  = np.array(list(wc_bounds.values()))
lambda_vals = {
    ci_label: np.array([ci_val / b for b in bounds_arr])
    for ci_label, ci_val in ci_map.items()
}

# ┌────────────────────────────────────────────────────────────────────────────┐
# │ 4) PRINT A TABLE OF RESULTS                                                │
# └────────────────────────────────────────────────────────────────────────────┘
ci_labels = list(ci_map.keys())
# Header
print(f"{'WC':<25}", end="")
for ci in ci_labels:
    header = ci.strip('$').replace("\\pi","π").replace("^2","²")
    print(f"{header:>12}", end="")
print("\n" + "-"*(25 + 12*len(ci_labels)))
# Rows
for i, wc in enumerate(wc_names):
    label_clean = wc.strip('$')
    print(f"{label_clean:<25}", end="")
    for ci in ci_labels:
        print(f"{lambda_vals[ci][i]:12.2f}", end="")
    print()

# ┌────────────────────────────────────────────────────────────────────────────┐
# │ 5) MAKE THE STACKED BAR PLOT                                                │
# └────────────────────────────────────────────────────────────────────────────┘
plt.figure(figsize=(9,6))
x = np.arange(len(wc_names))
bottom = np.zeros_like(bounds_arr)

for ci in ci_labels:
    heights = lambda_vals[ci]
    plt.bar(x, heights, bottom=bottom,
            color=colors[ci], label=f"c = {ci}")
    bottom += heights

# Axes and ticks
plt.yscale('log')
plt.xticks(x, wc_names, rotation=0, ha='center', fontsize=15)
plt.gca().tick_params(axis='x', which='both', length=0)
plt.ylabel(r'$\Lambda\, [\mathrm{TeV}]$', fontsize=15)
plt.xlabel('Wilson coefficient', fontsize=15, loc='right')
plt.ylim(1e-2, bottom.max()*1.2)

# ┌────────────────────────────────────────────────────────────────────────────┐
# │ 6) CMS STYLE ANNOTATIONS                                                    │
# └────────────────────────────────────────────────────────────────────────────┘
ax = plt.gca()
# point y‑ticks inwards and draw them on both left+right
ax.tick_params(axis='y', which='both',
               direction='in',
               left=True, right=True)
ax.text(0.0, 1.02, 'CMS',             transform=ax.transAxes,
        fontweight='bold', fontsize=20)
ax.text(0.12, 1.02, 'Preliminary',      transform=ax.transAxes,
        style='italic', fontsize=14)
#ax.text(0.50, 1.02, r'$95\%\ \mathrm{CL}\ \Lambda$ limits',transform=ax.transAxes, ha='center', fontsize=14)
ax.text(1.0, 1.02, r'138 fb$^{-1}$ (13 TeV)',
        transform=ax.transAxes, ha='right', fontsize=14)

plt.legend(frameon=False, loc="upper right", fontsize=15)
plt.tight_layout()
plt.savefig("/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v36/Energy_limits.pdf", dpi=300)
plt.show()


#c.SaveAs("/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v24/Energy_limits.pdf")

