#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Improved script to generate the constraint plot from ROOT files.
This script processes both "profiled" and "frozen" fits from ROOT trees,
extracts crossing points (for 1σ and 2σ intervals), and plots the results
with annotations. Adjust parameters as needed for flexibility.
"""

import numpy as np
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend for headless environments
import matplotlib.pyplot as plt
import ROOT as root
#from root_numpy import array2hist, hist2array  # Keeping import in case you need them
from scipy.interpolate import CubicSpline
import CombineHarvester.CombineTools.plotting as plot  # Imported if needed in future

def find_crossings(x, y, value):
    """
    Find approximate locations where the curve crosses a given 'value'.
    
    Uses a simple midpoint estimate between points that straddle the 'value'.
    
    Parameters:
        x (np.ndarray): Array of x values.
        y (np.ndarray): Array of y values.
        value (float): The threshold value to find crossings.
        
    Returns:
        np.ndarray: Array of x locations where the crossing occurs.
    """
    crossings = []
    for i in range(len(x) - 1):
        if (y[i] < value and y[i+1] > value) or (y[i] > value and y[i+1] < value):
            crossings.append(0.5 * (x[i] + x[i+1]))
    return np.array(crossings)

def tree_to_asymm_graph(tree, x_branch):
    """
    Convert a ROOT tree to a TGraphAsymmErrors object using the provided x_branch.
    
    Only points with 'quantileExpected' > -1.5 are processed. The y values are doubled.
    
    Parameters:
        tree (ROOT.TTree): The ROOT tree object.
        x_branch (str): Name of the branch used for x values.
        
    Returns:
        ROOT.TGraphAsymmErrors: Graph containing the selected points.
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
    
    Parameters:
        filename (str): Path to the ROOT file.
        x_branch (str): Name of the branch to use for x values.
        
    Returns:
        tuple: A dictionary containing crossing points, and numpy arrays of x and y values.
            The dictionary keys are:
                "one"  : crossings at y = 1.0,
                "two"  : crossings at y = 4.0,
                "min"  : x value corresponding to the minimum y.
    """
    file = root.TFile(filename)
    tree = file.Get("limit")
    graph = tree_to_asymm_graph(tree, x_branch)

    n_points = graph.GetN()
    x_values = []
    y_values = []
    for i in range(n_points):
        x_val = root.Double(0)
        y_val = root.Double(0)
        graph.GetPoint(i, x_val, y_val)
        x_values.append(x_val)
        y_values.append(y_val)
        
    x_values = np.array(x_values)
    y_values = np.array(y_values)
    
    crossings = {
        "one": find_crossings(x_values, y_values, 1.0),
        "two": find_crossings(x_values, y_values, 4.0),
        "min": x_values[np.argmin(y_values)]
    }
    
    return crossings, x_values, y_values

def plot_crossings(ax, offset, index, num_coeffs,
                   crossings_profiled, crossings_frozen,
                   label_prefix, color_profiled, color_frozen):
    """
    Plot crossing lines (for 1σ and 2σ) for both profiled and frozen fits.
    
    Parameters:
        ax (matplotlib.axes.Axes): Axis on which to plot.
        offset (float): Vertical offset multiplier.
        index (int): Index of the current coefficient.
        num_coeffs (int): Total number of coefficients.
        crossings_profiled (dict): Crossing points dictionary for the profiled fit.
        crossings_frozen (dict): Crossing points dictionary for the frozen fit.
        label_prefix (str): Text prefix for the legend.
        color_profiled (str): Color for the profiled crossing lines.
        color_frozen (str): Color for the frozen crossing lines.
    """
    # Define settings for 1σ ("one") and 2σ ("two")
    for key, linestyle, linewidth in [("one", "solid", 5), ("two", "dashed", 3)]:
        # Plot profiled crossings
        n_cross = len(crossings_profiled[key]) // 2
        for j in range(n_cross):
            label = "{} (profiled)".format(label_prefix) if (index == 0 and j == 0) else None
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
            label = "{} (frozen)".format(label_prefix) if (index == 0 and j == 0) else None
            ax.hlines(
                offset * (num_coeffs - index - 1) - 0.5,
                crossings_frozen[key][2 * j],
                crossings_frozen[key][2 * j + 1],
                color=color_frozen,
                linestyle=linestyle,
                linewidth=linewidth,
                label=label
            )

def plot_constraints(base_path, wc_list, wc_tex, factors, output_file,asimov_path=None, asimov_filenames=None, asimov_frozen_filenames=None):
    """
    Generates and saves the constraints plot.
    
    For each Wilson coefficient provided, the function retrieves the corresponding
    profiled and frozen fit data from ROOT files, scales the crossing points by the
    provided factors, and creates a plot with annotations.
    
    Parameters:
        base_path (str): Base directory for the frozen fit ROOT files.
        wc_list (list of str): List of Wilson coefficient branch names.
        wc_tex (list of str): List of LaTeX-formatted names for the coefficients.
        factors (list of float): Scaling factors for each coefficient.
        output_file (str): Full path of the output PDF file.
    """
    fig, ax = plt.subplots(figsize=(25, 15))
    plt.subplots_adjust(left=0.17, right=0.98, top=0.95, bottom=0.1)
    
    # Plot vertical reference lines
    ax.vlines(0, -2, 14, color="grey", linestyle="dashed", linewidth=2)
    ax.vlines(2, -2, 14, color="grey", linewidth=2)
    
    offset = 2.5
    num_coeffs = len(wc_list)
    
    for i, wc in enumerate(wc_list):
        # Build file paths for the profiled and frozen fits
        profiled_path = "v24_Signal/debug/higgsCombine_{}_combinedFit_profiled.root".format(wc.split("_")[1])
        frozen_path = "{}/higgsCombine_{}_combinedFit.MultiDimFit.mH120.root".format(base_path, wc.split("_")[1])
        
        # Retrieve graph crossing data from ROOT files
        crossings_profiled, _, _ = get_graph_data(profiled_path, wc)
        crossings_frozen, _, _ = get_graph_data(frozen_path, wc)
        
        # Load Asimov profiled & frozen crossings if given
        crossings_asimov = None
        crossings_asimov_frozen = None
        if asimov_path and asimov_filenames and asimov_frozen_filenames:
            asimov_file = "{}/{}".format(asimov_path, asimov_filenames[i])
            asimov_frozen_file = "{}/{}".format(asimov_path, asimov_frozen_filenames[i])
            crossings_asimov, _, _ = get_graph_data(asimov_file, wc)
            crossings_asimov_frozen, _, _ = get_graph_data(asimov_frozen_file, wc)
            for key in crossings_asimov:
                crossings_asimov[key] = factors[i] * crossings_asimov[key]
                crossings_asimov_frozen[key] = factors[i] * crossings_asimov_frozen[key]


        # Debug output of crossing values
        print("Coefficient {}: profiled crossings: {}".format(wc, crossings_profiled))
        print("Coefficient {}: frozen crossings: {}".format(wc, crossings_frozen))
        
        # Vertical position for the current coefficient
        y_position = offset * (num_coeffs - i - 1)
        
        # Annotate the floating (1σ) crossing information
        if crossings_profiled["one"].size >= 2:
            profiled_text = "[{}, {}]".format(
                round(crossings_profiled["one"][0], 3),
                round(crossings_profiled["one"][1], 3)
            )
        else:
            profiled_text = "N/A"
        ax.text(2.2, y_position, profiled_text, fontsize=22, color='black')
        
        # Annotate the 2σ crossing information
        if crossings_profiled["two"].size >= 4:
            text_two = "[{}, {}]U[{}, {}]".format(
                round(crossings_profiled["two"][0], 3),
                round(crossings_profiled["two"][1], 3),
                round(crossings_profiled["two"][2], 3),
                round(crossings_profiled["two"][3], 3)
            )
        elif crossings_profiled["two"].size >= 2:
            text_two = "[{}, {}]".format(
                round(crossings_profiled["two"][0], 3),
                round(crossings_profiled["two"][1], 3)
            )
        else:
            text_two = "N/A"
        ax.text(2.8, y_position, text_two, fontsize=22, color='black')
        
        # Static labels for 1σ and 2σ at the top of the plot
        ax.text(2.2, 14, r'1 $\sigma$ profiled', fontsize=22, color='black')
        ax.text(3.2, 14, r'2 $\sigma$ profiled', fontsize=22, color='black')
        
        # Scale the crossing values by the provided factor for the current coefficient
        for key in crossings_profiled.keys():
            crossings_profiled[key] = factors[i] * crossings_profiled[key]
            crossings_frozen[key] = factors[i] * crossings_frozen[key]

        # Draw Asimov expected limits as background bands (profiled: blue, frozen: gray)
        if crossings_asimov:
            for key, color, alpha in [("two", "lightskyblue", 0.25), ("one", "dodgerblue", 0.35)]:
                for j in range(len(crossings_asimov[key]) // 2):
                    ax.axvspan(
                        crossings_asimov[key][2*j],
                        crossings_asimov[key][2*j+1],
                        ymin=(y_position - 0.6)/20.0,
                        ymax=(y_position + 0.6)/20.0,
                        facecolor=color,
                        alpha=alpha,
                        zorder=0
                    )
        if crossings_asimov_frozen:
            for key, color, alpha in [("two", "lightgray", 0.25), ("one", "gray", 0.35)]:
                for j in range(len(crossings_asimov_frozen[key]) // 2):
                    ax.axvspan(
                        crossings_asimov_frozen[key][2*j],
                        crossings_asimov_frozen[key][2*j+1],
                        ymin=(y_position - 1.1)/20.0,
                        ymax=(y_position - 0.1)/20.0,
                        facecolor=color,
                        alpha=alpha,
                        zorder=0
                    )

        
        # Plot best-fit marker for profiled and frozen fits
        label_profiled = "Best profiled fit" if i == 0 else None
        ax.plot([crossings_profiled["min"]], [y_position],
                marker="o", markersize=10, markerfacecolor='red',
                markeredgecolor="black", linestyle='', label=label_profiled)
        
        label_frozen = "Best frozen fit" if i == 0 else None
        ax.plot([crossings_frozen["min"]], [y_position - 0.5],
                marker="o", markersize=10, markerfacecolor='maroon',
                markeredgecolor="black", linestyle='', label=label_frozen)
        
        # Plot crossing lines for both profiled and frozen fits
        plot_crossings(ax, offset, i, num_coeffs,
                       crossings_profiled,
                       crossings_frozen,
                       label_prefix="q < 1",
                       color_profiled="forestgreen",
                       color_frozen="black")
        
        # Annotate the Wilson coefficient name and scaling factor
        ax.text(-2.90, y_position - 0.25, wc_tex[i], fontsize=35, color='black')
        ax.text(-2.55, y_position - 0.25, "[x{}]".format(factors[i]), fontsize=30, color='black')
    
    # Create legend (only non-None labels are included)
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(-0.0, 0.9),
                       ncol=3, fontsize=30)
    legend.get_frame().set_facecolor('none')
    legend.get_frame().set_linewidth(0)
    
    # Format axes and ticks
    ax.tick_params(left=False, labelleft=False)
    ax.set_xlim(-2, 4)
    ax.set_xlabel("Wilson coefficient value", fontsize=35)
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
    ax.text(0.95, 0.95, 'CMS', transform=ax.transAxes,
            horizontalalignment='right', verticalalignment='top',
            fontsize=90, fontweight='bold')
    ax.text(0.95, 0.83, 'Preliminary', transform=ax.transAxes,
            horizontalalignment='right', verticalalignment='top',
            fontsize=40,style='italic')
    

    # Save the figure and clean up
    fig.savefig(output_file)
    plt.close(fig)
    print("Plot saved to", output_file)

def main():
    """Main entry point: sets parameters and creates the plot."""
    print("Plotting Constraints")
    
    # Define base directory, Wilson coefficients, scaling factors, and labels
    base_path = "v24_Signal/Allops_data"
    wc_list = ["k_ctt", "k_cQQ1", "k_cQt1", "k_cQt8", "k_ctHRe", "k_ctHIm"]
    factors = [1, 1, 1, 0.4, 0.05, 0.05]
    wc_tex = [r'$c_{tt}$', r'$c_{QQ}^{1}$', r'$c_{Qt}^{1}$',
              r'$c_{Qt}^{8}$', r'$c_{tH}^{\Re}$', r'$c_{tH}^{\Im}$']
    output_file = "/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v24/SummaryScan_new.pdf"
    
    asimov_path = "v24_Signal/Allops"

    asimov_filenames = [
        "higgsCombine_ctt_combinedFit_profiled.MultiDimFit.mH120.root",     
        "higgsCombine_cQQ1_combinedFit_profiled.MultiDimFit.mH120.root",     
        "higgsCombine_cQt1_combinedFit_profiled.MultiDimFit.mH120.root",     
        "higgsCombine_cQt8_combinedFit_profiled.MultiDimFit.mH120.root",     
        "higgsCombine_ctHRe_combinedFit_profiled.MultiDimFit.mH120.root",     
        "higgsCombine_ctHIm_combinedFit_profiled.MultiDimFit.mH120.root",     
    ]

    asimov_frozen_filenames = [
        "higgsCombine_ctt_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_cQQ1_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_cQt1_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_cQt8_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_ctHRe_combinedFit.MultiDimFit.mH120.root",     
        "higgsCombine_ctHIm_combinedFit.MultiDimFit.mH120.root",     
            ]

    plot_constraints(base_path, wc_list, wc_tex, factors, output_file, asimov_path=asimov_path, asimov_filenames=asimov_filenames,asimov_frozen_filenames=asimov_frozen_filenames)

if __name__ == "__main__":
    main()

