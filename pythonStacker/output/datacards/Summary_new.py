#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Improved script to generate the constraint plot from ROOT files (Python 2 version).
This script processes both "profiled" and "frozen" fits from ROOT trees,
extracts crossing points (for 1σ and 2σ intervals), and plots the results
with annotations. Adjust parameters as needed for flexibility.
"""

from __future__ import print_function  # optional: enables print() function in Python 2
import numpy as np
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend for headless environments
import matplotlib.pyplot as plt
import ROOT as root
#from root_numpy import array2hist, hist2array  # Keeping import in case you need them
from scipy.interpolate import CubicSpline
import CombineHarvester.CombineTools.plotting as plot  # Imported if needed in future
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle
import ctypes
from ctypes import byref


def find_crossings(x, y, value):
    """
    Find approximate locations where the curve crosses a given 'value'.
    
    Uses a simple midpoint estimate between points that straddle the 'value'.
    """
    crossings = []
    for i in range(len(x) - 1):
        if (y[i] < value and y[i+1] > value) or (y[i] > value and y[i+1] < value):
            crossings.append(0.5 * (x[i] + x[i+1]))

    # remove near-duplicates (within 0.1) but always keep the smallest and largest
    if crossings:
        crossings = sorted(crossings)
        filtered = [crossings[0]]
        for c in crossings[1:-1]:
            if abs(c - filtered[-1]) >= 0.1:
                filtered.append(c)
        filtered.append(crossings[-1])
        return np.array(filtered)

    return np.array(crossings)


import numpy as np

import numpy as np

crossings = {
    'cQQ1': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-0.563,  0.700]),
                'two': np.array([-0.860,  1.000]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-0.611,  0.709]),
                'two': np.array([-0.917,  1.023]),
            },
        },
        'data': {
            'frozen': {
                'min': 0.711,
                'one': np.array([-0.785, -0.102,  0.166,  1.006]),
                'two': np.array([-1.082,  1.260]),
            },
            'profiled': {
                'min': 0.650,
                'one': np.array([-0.789,  0.977]),
                'two': np.array([-1.131,  1.261]),
            },
        },
    },

    'cQt1': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-1.134,  0.976]),
                'two': np.array([-1.641,  1.477]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-1.159,  1.075]),
                'two': np.array([-1.684,  1.592]),
            },
        },
        'data': {
            'frozen': {
                'min': 0.850,
                'one': np.array([-1.539,  1.389]),
                'two': np.array([-1.991,  1.834]),
            },
            'profiled': {
                'min': 0.150,
                'one': np.array([-1.197,  1.231]),
                'two': np.array([-1.912,  1.835]),
            },
        },
    },

    'cQt8': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-1.949,  2.411]),
                'two': np.array([-2.981,  3.448]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-2.122,  2.417]),
                'two': np.array([-3.204,  3.466]),
            },
        },
        'data': {
            'frozen': {
                'min': -1.679,
                'one': np.array([-2.790,  3.219]),
                'two': np.array([-3.706,  4.156]),
            },
            'profiled': {
                'min': -0.300,
                'one': np.array([-2.461,  2.503]),
                'two': np.array([-3.770,  3.962]),
            },
        },
    },

    'ctHIm': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-11.847, 12.214]),
                'two': np.array([-17.196, 17.543]),
            },
            'profiled': {
                'min': -12.500,
                'one': np.array([-20.929, 14.547]),
                'two': np.array([-25.952, 21.649]),
            },
        },
        'data': {
            'frozen': {
                'min': -6.800,
                'one': np.array([-14.091, 13.646]),
                'two': np.array([-19.128, 19.098]),
            },
            'profiled': {
                'min': -12.500,
                'one': np.array([-22.000, 15.500]),
                'two': np.array([-26.000, 21.500]),
            },
        },
    },

    'ctHRe': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-2.607,  3.229]),
                'two': np.array([-5.037,  7.867]),
            },
            'profiled': {
                'min': 2.500,
                'one': np.array([-2.134, 11.129]),
                'two': np.array([-5.171, 30.277]),
            },
        },
        'data': {
            'frozen': {
                'min': -0.188,
                'one': np.array([-2.962,  3.136]),
                'two': np.array([-5.590,  8.060]),
            },
            'profiled': {
                'min': 2.500,
                'one': np.array([-3.000, 14.500]),
                'two': np.array([-6.000, 31.000]),
            },
        },
    },

    'ctt': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-0.624,  0.739]),
                'two': np.array([-0.948,  1.067]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-0.679,  0.753]),
                'two': np.array([-1.012,  1.090]),
            },
        },
        'data': {
            'frozen': {
                'min': 0.660,
                'one': np.array([-0.872,  1.010]),
                'two': np.array([-1.172,  1.301]),
            },
            'profiled': {
                'min': 0.250,
                'one': np.array([-0.753,  0.888]),
                'two': np.array([-1.188,  1.262]),
            },
        },
    },
}



def lookup_crossings(crossings_dict, coeff, dataset, mode):
    """
    Retrieve crossings for a given coefficient (e.g. 'ctt'),
    dataset ('data' or 'Asimov'), mode ('profiled' or 'frozen'),
    and level ('one', 'two', 'min'). Returns a numpy array.
    """
    try:
        vals = crossings_dict[coeff][dataset][mode]
    except KeyError:
        return {}
    return vals


def tree_to_asymm_graph(tree, x_branch):
    """
    Convert a ROOT tree to a TGraphAsymmErrors object using the provided x_branch.
    
    Only points with 'quantileExpected' > -1.5 are processed. The y values are doubled.
    """
    graph = root.TGraphAsymmErrors()
    point_index = 0
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        quantile = getattr(tree, "quantileExpected")
        if quantile > -1.5:
            x_value = getattr(tree, x_branch)
            y_value = getattr(tree, "deltaNLL")
            graph.SetPoint(point_index, x_value, 2 * y_value)
            point_index += 1
    graph.Sort()
    return graph


def get_graph_data(filename, x_branch):
    """
    Load a ROOT file, convert its 'limit' tree to a graph, and extract data.
    
    Also computes crossing points for y = 1 and y = 4 as well as the minimum.
    
    Returns:
        tuple: A dictionary containing crossing points, and numpy arrays of x and y values.
               The dictionary keys are:
                  "one"  : crossings at y = 1.0,
                  "two"  : crossings at y = 4.0,
                  "min"  : x value corresponding to the minimum y.
    """
    # open file and make graph
    f    = root.TFile.Open(filename)
    tree = f.Get("limit")
    graph = tree_to_asymm_graph(tree, x_branch)

    # number of points
    n = graph.GetN()

    # Access the internal X/Y arrays directly
    x_arr = graph.GetX()   # pointer-like sequence of length n
    y_arr = graph.GetY()

    # Build numpy arrays in one line each
    x_values = np.array([x_arr[i] for i in range(n)])
    y_values = np.array([y_arr[i] for i in range(n)])

    # find crossings and minimum
    crossings = {
        "one": find_crossings(x_values, y_values, 1.0),
        "two": find_crossings(x_values, y_values, 4.0),
        "min": x_values[np.argmin(y_values)],
    }

    return crossings, x_values, y_values


def plot_crossings(ax, offset, index, num_coeffs,
                   crossings_profiled, crossings_frozen,
                   label_prefix, color_profiled, color_frozen):
    """
    Plot crossing lines (for 1σ and 2σ) for both profiled and frozen fits.
    """
    # Define labels for profiled/frozen lines keyed by "one" vs "two"
    profiled_labels = {
        "one": "q < 1 observed(profiled)",
        "two": "q < 4 observed(profiled)"
    }
    frozen_labels = {
#        "one": "q < 1 observed(frozen)",
#        "two": "q < 4 observed(frozen)"
        "one": "q < 1 observed",
        "two": "q < 4 observed"
    }

    # Define settings for 1σ ("one") and 2σ ("two")
    for key, linestyle, linewidth in [("one", "solid", 5), ("two", "dashed", 3)]:
        # Plot profiled crossings
        n_cross = len(crossings_profiled[key]) // 2
        for j in range(n_cross):
            label = "__no_legend__"#profiled_labels[key] if (index == 0 and j == 0) else None 
            ax.hlines(
                offset * (num_coeffs - index - 1),
                crossings_profiled[key][2 * j],
                crossings_profiled[key][2 * j + 1],
                color=color_profiled,
                linestyle=linestyle,
                linewidth=linewidth,
                label=label
            )
        # Plot frozen crossings
        n_cross_frozen = len(crossings_frozen[key]) // 2
        for j in range(n_cross_frozen):
            label = frozen_labels[key] if (index == 0 and j == 0) else None 
            ax.hlines(
                offset * (num_coeffs - index - 1) - 0.5,
                crossings_frozen[key][2 * j],
                crossings_frozen[key][2 * j + 1],
                color=color_frozen,
                linestyle=linestyle,
                linewidth=linewidth,
                label=label
            )

def plot_constraints(base_path, wc_list, wc_tex, factors, output_file,
                     asimov_path=None, asimov_filenames=None, asimov_frozen_filenames=None):
    """
    Generates and saves the constraints plot.
    
    For each Wilson coefficient provided, the function retrieves the corresponding
    profiled and frozen fit data from ROOT files, scales the crossing points by the
    provided factors, and creates a plot with annotations.
    """
    # Increase vertical space
    fig, ax = plt.subplots(figsize=(25, 16))
    plt.subplots_adjust(left=0.17, right=0.98, top=0.95, bottom=0.1)
    
    # Plot vertical reference lines
    ax.vlines(0, -2, 14, color="grey", linestyle="dashed", linewidth=2)
    ax.vlines(2, -2, 14, color="grey", linewidth=2)
    
    offset = 2.5
    num_coeffs = len(wc_list)
    
    # Colors for Asimov boxes (profiled / frozen, 2σ / 1σ)
    color_profiled_2sigma = "palegreen"
    color_profiled_1sigma = "limegreen"#"mediumspringgreen"
    color_frozen_2sigma   = "lightgray"
    color_frozen_1sigma   = "darkgray"
    
    # Half-heights for the vertical boxes, for alignment
    box_half_height_profiled = 0.25
    box_half_height_frozen   = 0.25


    for i, wc in enumerate(wc_list):
        # Build file paths for profiled and frozen fits using .format
        profiled_path = "./v36_IM/Allops/higgsCombine_{0}_combinedFit_profiled_data.root".format(wc)#(wc.split("_")[1])
        frozen_path   = "./v36_IM/Allops/higgsCombine_{0}_combinedFit_data.MultiDimFit.mH120.root".format(wc)#(wc.split("_")[1])
        
        # Retrieve graph crossing data from ROOT files
        crossings_profiled, _, _ = get_graph_data(profiled_path, wc)
        crossings_frozen, _, _   = get_graph_data(frozen_path, wc)
        
        # Load Asimov expected crossings if given
        crossings_asimov = None
        crossings_asimov_frozen = None
        if asimov_path and asimov_filenames and asimov_frozen_filenames:
            asimov_file       = "{0}/{1}".format(asimov_path, asimov_filenames[i])
            asimov_frozen_file = "{0}/{1}".format(asimov_path, asimov_frozen_filenames[i])
            crossings_asimov, _, _        = get_graph_data(asimov_file, wc)
            crossings_asimov_frozen, _, _ = get_graph_data(asimov_frozen_file, wc)
        
        print("Coefficient {}: Asimov profiled crossings: {}".format(wc, crossings_asimov))
        print("Coefficient {}: Asimov frozen crossings: {}".format(wc, crossings_asimov_frozen))

        # now replace them 
        crossings_profiled = lookup_crossings(crossings, wc, 'data', 'profiled') 
        crossings_frozen = lookup_crossings(crossings, wc, 'data', 'frozen') 
        crossings_asimov = lookup_crossings(crossings, wc, 'Asimov', 'profiled') 
        crossings_asimov_frozen = lookup_crossings(crossings, wc, 'Asimov', 'frozen') 
        
        # Debug output of crossing values
        print("Coefficient {}: data profiled crossings: {}".format(wc, crossings_profiled))
        print("Coefficient {}: data frozen crossings: {}".format(wc, crossings_frozen))
        #print("Coefficient {}: Asimov profiled crossings: {}".format(wc, crossings_asimov))
        #print("Coefficient {}: Asimov frozen crossings: {}".format(wc, crossings_asimov_frozen))

        # Vertical position for the current coefficient
        y_position = offset * (num_coeffs - i - 1)
        
        # Annotate the floating (1σ) crossing information
        if crossings_profiled["one"].size >= 2:
            profiled_text = "[{}, {}]".format(
                round(crossings_profiled["one"][0], 2),
                round(crossings_profiled["one"][1], 2)
            )
        else:
            profiled_text = "N/A"
        ax.text(2.4, y_position, profiled_text, fontsize=28, color='black', ha='center')
        
        # Annotate the 2σ crossing information
        if crossings_profiled["two"].size >= 4:
            text_two = "[{}, {}]U[{}, {}]".format(
                round(crossings_profiled["two"][0], 2),
                round(crossings_profiled["two"][1], 2),
                round(crossings_profiled["two"][2], 2),
                round(crossings_profiled["two"][3], 2)
            )
            ax.text(3.4, y_position, text_two, fontsize=20, color='black', ha='center')
        elif crossings_profiled["two"].size >= 2:
            text_two = "[{}, {}]".format(
                round(crossings_profiled["two"][0], 2),
                round(crossings_profiled["two"][1], 2)
            )
            ax.text(3.4, y_position, text_two, fontsize=28, color='black', ha='center')
        else:
            text_two = "N/A"
            ax.text(3.4, y_position, text_two, fontsize=28, color='black', ha='center')
        
        # Static labels for 1σ and 2σ at the top of the plot
        ax.text(2.4, 14, r'1 $\sigma$ profiled', fontsize=28, color='black', ha='center')
        ax.text(3.4, 14, r'2 $\sigma$ profiled', fontsize=28, color='black', ha='center')
        
        # Scale the crossing values by the provided factor for the current coefficient
        for key in crossings_profiled.keys():
            crossings_profiled[key] = factors[i] * crossings_profiled[key]
            crossings_frozen[key] = factors[i] * crossings_frozen[key]
        for key in crossings_asimov:
            crossings_asimov[key] = factors[i] * crossings_asimov[key]
            crossings_asimov_frozen[key] = factors[i] * crossings_asimov_frozen[key]
        
        # Draw Asimov expected limits as background boxes using Rectangle patches
        if crossings_asimov:
            for key, color, alpha in [("two", color_profiled_2sigma, 0.4), ("one", color_profiled_1sigma, 0.5)]:
                for j in range(len(crossings_asimov[key]) // 2):
                    x_min = crossings_asimov[key][2*j]
                    x_max = crossings_asimov[key][2*j+1]
                    rect = Rectangle(
                        (x_min, y_position - box_half_height_profiled),  # lower left corner
                        x_max - x_min,                                    # width
                        2 * box_half_height_profiled,                     # height
                        facecolor=color,
                        alpha=alpha,
                        zorder=0
                    )
                    
                    ax.add_patch(rect)

        if crossings_asimov_frozen:
            for key, color, alpha in [("two", color_frozen_2sigma, 0.4), ("one", color_frozen_1sigma, 0.5)]:
                for j in range(len(crossings_asimov_frozen[key]) // 2):
                    x_min = crossings_asimov_frozen[key][2*j]
                    x_max = crossings_asimov_frozen[key][2*j+1]
                    rect = Rectangle(
                        (x_min, (y_position - 0.5) - box_half_height_frozen),  # lower left corner for frozen (offset by 0.5)
                        x_max - x_min,                                         # width
                        2 * box_half_height_frozen,                           # height
                        facecolor=color,
                        alpha=alpha,
                        zorder=0
                    )
                    
                    ax.add_patch(rect)


        # Plot best-fit markers for profiled and frozen fits
        label_profiled = "Best fit (profiled)" if i == 0 else None
        ax.plot([crossings_profiled["min"]], [y_position],
                marker="o", markersize=12, markerfacecolor="forestgreen",
                markeredgecolor="forestgreen", linestyle='', label=label_profiled)
        
        label_frozen = "Best fit (frozen)" if i == 0 else None
        ax.plot([crossings_frozen["min"]], [y_position - 0.5],
                marker="o", markersize=12, markerfacecolor='black',
                markeredgecolor="black", linestyle='', label=label_frozen)
        
        # Plot crossing lines for both profiled and frozen fits
        plot_crossings(ax, offset, i, num_coeffs,
                       crossings_profiled,
                       crossings_frozen,
                       label_prefix="q < 1",
                       color_profiled="forestgreen",
                       color_frozen="black")
        
        # Annotate the Wilson coefficient name and scaling factor
        ax.text(-2.90, y_position - 0.25, wc_tex[i], fontsize=38, color='black')
        ax.text(-2.55, y_position - 0.25, "[x{}]".format(factors[i]), fontsize=30, color='black')
    

    # Get already existing handles and labels (if any)
    handles, labels = ax.get_legend_handles_labels()
    
    # Create custom patches for the Asimov boxes.
    # Note: use the same colors and alpha values as used in your loops.
    custom_handles = [
#        Patch(facecolor=color_profiled_1sigma, alpha=0.5, label="q < 1 expected(profiled)"),
#        Patch(facecolor=color_frozen_1sigma,   alpha=0.5, label="q < 1 expected(frozen)"),
#        Patch(facecolor=color_profiled_2sigma, alpha=0.4, label="q < 4 expected(profiled)"),
#        Patch(facecolor=color_frozen_2sigma,   alpha=0.4, label="q < 4 expected(frozen)")
        Patch(facecolor=color_frozen_1sigma,   alpha=0.5, label="q < 1 expected"),
        Patch(facecolor=color_frozen_2sigma,   alpha=0.4, label="q < 4 expected")
    ]
    
    # Append these custom patches to your existing legend entries.
    handles.extend(custom_handles)
    labels.extend([h.get_label() for h in custom_handles])


    # Create legend
    legend = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.0, 0.98),
                       ncol=3, fontsize=38)
    legend.get_frame().set_facecolor('none')
    legend.get_frame().set_linewidth(0)
    
    # Format axes and ticks
    ax.tick_params(left=False, labelleft=False)
    ax.set_xlim(-2, 4)
    ax.set_xlabel("Wilson coefficient value", fontsize=38)
    ax.xaxis.set_label_coords(0.84, -0.05)  # Adjust the y coordinate as needed
    ax.set_ylim(-2, 20)
    
    # Customize x-axis ticks
    x_ticks = np.arange(-2, 2.2, 0.2)
    x_tick_labels = [""] * len(x_ticks)
    if len(x_ticks) > 0:
        x_tick_labels[0] = str(round(x_ticks[0], 2))
    if len(x_ticks) > 5:
        x_tick_labels[5] = str(round(x_ticks[5], 2))
    if len(x_ticks) > 10:
        x_tick_labels[10] = "0"
    if len(x_ticks) > 15:
        x_tick_labels[15] = str(round(x_ticks[15], 2))
    if len(x_ticks) > 20:
        x_tick_labels[20] = str(round(x_ticks[20], 2))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.tick_params(axis='both', direction='in', which='both', length=10, width=2)
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='x', labelsize=30)
    
    # Add CMS label in the top-right corner with 'Preliminary' beneath it.
    ax.text(-2.1, 19, 'CMS', #transform=ax.transAxes,
            horizontalalignment='right', verticalalignment='top',
            fontsize=90, fontweight='bold', ha='right')
    ax.text(-2.1, 17, 'Preliminary', #transform=ax.transAxes,
            horizontalalignment='right', verticalalignment='top',
            fontsize=35, style='italic', ha='right')
    
    # Save and close figure
    fig.savefig(output_file)
    plt.close(fig)
    print("Plot saved to", output_file)

def main():
    """Main entry point: sets parameters and creates the plot."""
    print("Plotting Constraints")
    
    base_path = "v36_Signal/Allops"
    wc_list = ["ctt", "cQQ1", "cQt1", "cQt8", "ctHRe", "ctHIm"]
    factors = [1, 1, 1, 0.4, 0.05, 0.05]
    wc_tex = [r'$c_{tt}$', r'$c_{QQ}^{(1)}$', r'$c_{Qt}^{(1)}$',
              r'$c_{Qt}^{(8)}$', r'$c_{tH}^{Re}$', r'$c_{tH}^{Im}$']
    output_file = "/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v36/SummaryScan_new.pdf"
    
    asimov_path = "v36_IM/Allops"
    asimov_filenames = [
        "higgsCombine_ctt_combinedFit_profiled.root",     
        "higgsCombine_cQQ1_combinedFit_profiled.root",     
        "higgsCombine_cQt1_combinedFit_profiled.root",     
        "higgsCombine_cQt8_combinedFit_profiled.root",     
        "higgsCombine_ctHRe_combinedFit_profiled.root",     
        "higgsCombine_ctHIm_combinedFit_profiled.root"
    ]
    
    asimov_frozen_filenames = [
        "higgsCombine_ctt_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_cQQ1_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_cQt1_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_cQt8_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_ctHRe_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_ctHIm_combinedFit.MultiDimFit.mH120.root"
    ]
    
    plot_constraints(base_path, wc_list, wc_tex, factors, output_file,
                     asimov_path=asimov_path,
                     asimov_filenames=asimov_filenames,
                     asimov_frozen_filenames=asimov_frozen_filenames)

if __name__ == "__main__":
    main()

