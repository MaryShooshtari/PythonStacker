import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D
import os,shutil

import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator


import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.path import Path
from matplotlib.patches import PathPatch


def load_data(path):
    """
    Load JSON file containing 2D limit contours and best fit.
    Expected structure:
    {
      "1sigma": [{"x": [...], "y": [...]}],
      "2sigma": [{"x": [...], "y": [...]}],
      "bestFit": {"x": ..., "y": ...}
    }
    """
    with open(path, 'r') as f:
        return json.load(f)

from matplotlib.patches import FancyArrowPatch

def _connect_arc(ax, p1, p2, *, color='black', lw=1, ls='solid', rad=0.2):
    """
    Draw a smooth circular arc from p1 → p2.
    'rad' is the curvature: positive bows one way, negative the other.
    """
    arc = FancyArrowPatch(
        p1, p2,
        arrowstyle='-',
        connectionstyle=f'arc3,rad={rad}',
        color=color, linewidth=lw, linestyle=ls
    )
    ax.add_patch(arc)


def _connect_bezier(ax, p1, p2, *, color='black', lw=1, ls='solid', frac=0.1):
    """
    Draw a smooth quadratic Bezier from p1→p2 by placing the control
    point at the midpoint plus a small perpendicular offset.
    """
    x1, y1 = p1
    x2, y2 = p2
    mx, my = 0.5*(x1+x2), 0.5*(y1+y2)
    dx, dy = x2-x1, y2-y1
    nx, ny = -dy, dx
    norm = np.hypot(nx, ny) or 1.0
    nx, ny = nx/norm, ny/norm
    offset = frac * np.hypot(dx, dy)
    cx, cy = mx + nx*offset, my + ny*offset

    verts = [(x1, y1), (cx, cy), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    path = Path(verts, codes)
    patch = PathPatch(path, facecolor='none', edgecolor=color,
                      lw=lw, linestyle=ls)
    ax.add_patch(patch)

from scipy.optimize import linear_sum_assignment  # SciPy

def global_close(ax, ends, *, color, lw, ls, rad):
    """
    ends: list of (x,y) endpoints to be paired off 1↔2, 3↔4, …
    Draws a circular‐arc between each matched pair.
    """
    pts = np.array(ends)                         # shape (N,2)
    N = len(ends)
    # build cost matrix of squared distances
    D2 = np.sum((pts[:,None,:] - pts[None,:,:])**2, axis=2)
    np.fill_diagonal(D2, np.inf)                 # forbid self‐match

    # solve assignment: one row→one col, total sum minimized
    row_idx, col_idx = linear_sum_assignment(D2)

    # only take each pair once (i<j)
    for i, j in zip(row_idx, col_idx):
        if i < j:
            x1, y1 = pts[i]
            x2, y2 = pts[j]
            # draw an arc with FancyArrowPatch
            arc = FancyArrowPatch(
                (x1,y1), (x2,y2),
                arrowstyle='-',
                connectionstyle=f'arc3,rad={rad}',
                color=color, linewidth=lw, linestyle=ls
            )
            ax.add_patch(arc)



def draw_and_close_contours(ax, segments, color, style, connect_kwargs=None):
    """
    Plot each segment and close all open ends by connecting nearest neighbors.
    """
    # draw raw segments
    for seg in segments:
        ax.plot(seg['x'], seg['y'], color=color,
                lw=style['lw'], ls=style['ls'])

    # collect all endpoints
    ends = []
    for seg in segments:
        if not seg['x']: continue
        ends += [(seg['x'][0], seg['y'][0]), (seg['x'][-1], seg['y'][-1])]

#    #    globally match & draw
#        global_close(ax, ends,
#             color=color, 
#             lw=style['lw'], ls=style['ls'],
#             rad=connect_kwargs.get('rad', 0.1))

    # pair and connect closest ends
    get_frac = lambda: (connect_kwargs or {}).get('frac', 0.1)
    while ends:
        x1, y1 = ends.pop(0)
        # find nearest remaining end
        idx, (x2, y2) = min(
            enumerate(ends),
            key=lambda iv: (iv[1][0]-x1)**2 + (iv[1][1]-y1)**2
        )
        ends.pop(idx)
#        _connect_bezier(
#            ax, (x1, y1), (x2, y2),
#            color=color,
#            lw=style['lw'], ls=style['ls'],
#            frac=get_frac()
#        )
        _connect_arc(
            ax, (x1, y1), (x2, y2),
            color=color,
            lw=style['lw'], ls=style['ls'],
            rad=connect_kwargs.get('rad', -0.08)
        )


def mathify(lbl):
    lbl = lbl.strip()
    return lbl if lbl.startswith('$') and lbl.endswith('$') else f'${lbl}$'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot 2D limits from JSON with closed contours and CMS styling.')
    parser.add_argument('--obs', required=True, help='Observed JSON file')
    parser.add_argument('--exp', required=True, help='Expected JSON file')
    parser.add_argument('--xlabel', default=r'$\kappa_t$', help='X-axis label')
    parser.add_argument('--ylabel', default=r'$\tilde{\kappa}_t$', help='Y-axis label')
    parser.add_argument('--xrange', nargs=2, type=float, default=[-2.4, 2.4],
                        help='X-axis limits')
    parser.add_argument('--yrange', nargs=2, type=float, default=[-2.4, 2.4],
                        help='Y-axis limits')
    parser.add_argument('--output', default='yukawa_limits.png',
                        help='Output filename')
    args = parser.parse_args()

    args.xlabel = mathify(args.xlabel)
    args.ylabel = mathify(args.ylabel)

    obs_data = load_data(args.obs)
    exp_data = load_data(args.exp)

    plt.style.use(hep.style.CMS)
    fig, ax = plt.subplots(figsize=(8,8))

    # observed in limegreen
    draw_and_close_contours(ax, obs_data['1sigma'], 'limegreen',
                            style={'lw':2,'ls':'solid'},
                            connect_kwargs={'frac':0.05})
    draw_and_close_contours(ax, obs_data['2sigma'], 'limegreen',
                            style={'lw':2,'ls':'dashed'},
                            connect_kwargs={'frac':0.05})

    # expected in black
    draw_and_close_contours(ax, exp_data['1sigma'], 'black',
                            style={'lw':2,'ls':'solid'},
                            connect_kwargs={'frac':0.05})
    draw_and_close_contours(ax, exp_data['2sigma'], 'black',
                            style={'lw':2,'ls':'dashed'},
                            connect_kwargs={'frac':0.05})

    # best-fit markers
    bf_obs = obs_data.get('bestFit', {})
    ax.plot(bf_obs.get('x'), bf_obs.get('y'), marker='x', color='limegreen',
            markersize=10, linestyle='None')
    bf_exp = exp_data.get('bestFit', {})
    ax.plot(bf_exp.get('x'), bf_exp.get('y'), marker='+', color='black',
            markersize=12, markeredgewidth=2, linestyle='None')

    # axes
    ax.set_xlim(args.xrange)
    ax.set_ylim(args.yrange)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
    ax.set_xlabel(args.xlabel, fontsize=24)
    ax.set_ylabel(args.ylabel, fontsize=24)
    ax.tick_params(direction='in', which='both', top=False, right=False,
                   labelsize=14)
    ax.tick_params(labeltop=False, labelright=False)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)



    # 2) first legend (Expected → q levels + SM), top-left
    expected_handles = [
        Line2D([0],[0], color='black', lw=2, ls='solid', label=r'$q<2.3$'),
        Line2D([0],[0], color='black', lw=2, ls='dashed', label=r'$q<6.2$'),
        Line2D([0],[0], marker='+', color='black', ls='None', markersize=12, label='SM'),
        Line2D([0], [0], marker='x', color='limegreen', lw=0, markersize=10, label='bestfit'),
    ]
    leg1 = ax.legend(
        handles=expected_handles,
        #title='Profiled',
        #title_fontsize=18,
        loc='upper center',
        bbox_to_anchor=(0.5, 0.965),  # tweak as needed
        columnspacing=1.0,   # ← default is 2.0, make it smaller
        handletextpad=0.4,   # ← space between the line marker and its text
        frameon=False,
        fontsize=18,
        ncol=4,
    )
    ax.add_artist(leg1)  # <— lock this legend in place
    
    # 3) second legend (hypotheses → colored solid lines), bottom-center
    hypo_handles = [
        Line2D([0], [0], color='limegreen', lw=2, linestyle='solid', 
               label='observed'),
        Line2D([0], [0], color='black', lw=2, linestyle='solid', 
               label='expected'),
    ]
    ax.legend(
        handles=hypo_handles,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.035),
        frameon=False,
        fontsize=18,
        ncol=2,
    )
    


    # CMS Preliminary label
    hep.cms.label(ax=ax, data=True,
                  label='Preliminary',
                  rlabel='138 fb$^{-1}$ (13 TeV)')


    # Improve layout and save
    # Save figure
    plot_directory_ = "/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v36/2D_limits/"
    if not os.path.exists(plot_directory_):
       try:
           os.makedirs(plot_directory_)
       except OSError: # Resolve rare race condition
           pass
    
    shutil.copyfile('/user/mshoosht/public_html/index.php', os.path.join( plot_directory_, 'index.php' ) )

    plt.tight_layout()
    fig.savefig(plot_directory_+args.output+".pdf", bbox_inches='tight', dpi=300)
    fig.savefig(plot_directory_+args.output+".png", bbox_inches='tight', dpi=300)
    print("saved plot to : ", plot_directory_+args.output+".pdf")

