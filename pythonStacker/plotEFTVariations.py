import numpy as np
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))
import awkward as ak
import json
import argparse
import matplotlib.pyplot as plt
import os

from src.variables.variableReader import VariableReader, Variable
from src.configuration import load_channels
from src.histogramTools import HistogramManager

import src.plotTools.figureCreator as fg
from src import generate_binning
from plotHistogramsRedo import modify_yrange_updown, generate_outputfolder, copy_index_html

import plugins.eft as eft

import src.arguments as arguments


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process command line arguments.')

    arguments.add_settingfiles(parser)
    arguments.select_specifics(parser)
    arguments.add_tmp_storage(parser)
    arguments.add_plot_output(parser)

    args = parser.parse_args()
    return args


def lin_quad_plot_EFT(variable: Variable, plotdir: str, histograms, process_info: dict, plotlabel: str, processname: str):
    fig, (ax_main, ax_ratio_one, ax_ratio_two) = fg.create_multi_ratioplot(n_subplots=2)
    binning = generate_binning(variable.range, variable.nbins)
    # first plot nominal, then start adding variations
    print(variable.name)
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    #stat_unc_var = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name]["stat_unc"])) / nominal_content, nan=0.)

    # pretty_name = generate_process_name("SM", info)
    nominal_weights = np.ones(len(nominal_content))
    ax_main.hist(binning[:-1], binning, weights=nominal_weights, histtype="step", color="k", label="SM")

    eft_variations = eft.getEFTVariationsLinear()
    minim = 1.
    maxim = 1.
    for eft_var in eft_variations:
        wc_factor = 1
        if "ctHRe" in eft_var:
            wc_factor = 20
        if "ctHIm" in eft_var:
            wc_factor = 20
        lin_name = "EFT_" + eft_var
        quad_name = lin_name + "_" + eft_var

        current_variation = nominal_content
        current_variation = current_variation + wc_factor * np.array(ak.to_numpy(histograms[variable.name][lin_name]["Up"]))
        current_variation = current_variation + wc_factor * wc_factor * np.array(ak.to_numpy(histograms[variable.name][quad_name]["Up"]))

        current_variation = np.nan_to_num(current_variation / nominal_content, nan=1.)

        minim = min(minim, np.min(current_variation))
        maxim = max(maxim, np.max(current_variation))
        pretty_eft_name = eft_var + f" = {wc_factor}"
        ax_main.hist(binning[:-1], binning, weights=current_variation, histtype="step", label=pretty_eft_name)

        lin_ratio = np.nan_to_num(wc_factor * np.array(ak.to_numpy(histograms[variable.name][lin_name]["Up"])) / nominal_content, nan=0.)
        quad_ratio = np.nan_to_num(wc_factor * wc_factor * np.array(ak.to_numpy(histograms[variable.name][quad_name]["Up"])) / nominal_content, nan=0.)
        ax_ratio_one.hist(binning[:-1], binning, weights=lin_ratio, histtype="step")
        ax_ratio_two.hist(binning[:-1], binning, weights=quad_ratio, histtype="step")
    
    mix_list = ["cQQ1_cQt1","cQQ1_cQt8","cQQ1_ctHIm","cQQ1_ctHRe","cQQ1_ctt",
		    "cQQ8_cQQ1","cQQ8_cQt1",
		    #"cQQ8_cQt8","cQQ8_ctHIm","cQQ8_ctHRe","cQQ8_ctt","cQt1_cQt8","cQt1_ctHIm","cQt1_ctHRe","cQt1_ctt","cQt8_ctHIm","cQt8_ctHRe","ctHRe_ctHIm","ctt_cQt8","ctt_ctHIm","ctt_ctHRe"
		    ]
    for eft_var in mix_list:
        mix_name = "EFT_"+eft_var
        mix_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][mix_name]["Up"])) / nominal_content, nan=0.)
        #ax_ratio_two.hist(binning[:-1], binning, weights=mix_ratio, histtype="step")
 

    ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=np.ones(len(nominal_content)), ecolor='k', label="stat unc.")

    ax_main.set_xlim(variable.range)
    ax_main.set_ylabel("SM + EFT / SM")
    modify_yrange_updown(ax_main, (minim, maxim), up_scale=1.2)
    ax_main.legend(ncol=2)
    ax_main.text(0.059, 0.74, plotlabel, transform=ax_main.transAxes)

    ax_ratio_one.set_xlim(variable.range)
    ax_ratio_two.set_xlim(variable.range)

    ax_ratio_one.set_ylabel("Lin / SM")
    ax_ratio_two.set_ylabel("quad / SM")
    ax_ratio_two.set_xlabel(variable.axis_label)

    # fix output name
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_ratios.png"))
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_ratios.pdf"))
    plt.close(fig)


def main_plot_EFT(variable: Variable, plotdir: str, histograms, process_info: dict, plotlabel: str, processname: str):
    fig, ax_main = fg.create_singleplot()

    binning = generate_binning(variable.range, variable.nbins)
    # first plot nominal, then start adding variations
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    stat_unc_var = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name]["stat_unc"])) / nominal_content, nan=0.)
    # pretty_name = generate_process_name("SM", info)
    nominal_weights = np.ones(len(nominal_content))
    ax_main.hist(binning[:-1], binning, weights=nominal_weights, histtype="step", color="k", label="SM")
    print(nominal_content)

    eft_variations = eft.getEFTVariationsLinear()
    minim = 1.
    maxim = 1.
    for eft_var in eft_variations:
        wc_factor = 1
        if "ctHRe" in eft_var:
            wc_factor = 20
        if "ctHIm" in eft_var:
            wc_factor = 20
        lin_name = "EFT_" + eft_var
        quad_name = lin_name + "_" + eft_var

        current_variation = nominal_content
        current_variation = current_variation + wc_factor * np.array(ak.to_numpy(histograms[variable.name][lin_name]["Up"]))
        current_variation = current_variation + wc_factor * wc_factor * np.array(ak.to_numpy(histograms[variable.name][quad_name]["Up"]))

        current_variation = np.nan_to_num(current_variation / nominal_content, nan=1.)

        minim = min(minim, np.min(current_variation))
        maxim = max(maxim, np.max(current_variation))
        pretty_eft_name = eft_var + f" = {wc_factor}"
        ax_main.hist(binning[:-1], binning, weights=current_variation, histtype="step", label=pretty_eft_name)


    ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=np.ones(len(nominal_content)), yerr=stat_unc_var, ecolor='k', label="stat unc.")
    ax_main.set_xlim(variable.range)
    ax_main.set_ylabel("SM + EFT / SM")
    modify_yrange_updown(ax_main, (minim, maxim), up_scale=1.4)
    ax_main.legend(ncol=2)
    ax_main.set_xlabel(variable.axis_label)
    ax_main.text(0.049, 0.77, plotlabel, transform=ax_main.transAxes)

    # fix output name
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}.png"))
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}.pdf"))
    plt.close(fig)

def mix_plot_EFT(variable: Variable, plotdir: str, histograms, process_info: dict, plotlabel: str, processname: str):
    fig, ax_main = fg.create_singleplot()

    binning = generate_binning(variable.range, variable.nbins)
    # first plot nominal, then start adding variations
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    #stat_unc_var = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name]["stat_unc"])) / nominal_content, nan=0.)
    # pretty_name = generate_process_name("SM", info)
    nominal_weights = np.ones(len(nominal_content))
    #ax_main.hist(binning[:-1], binning, weights=nominal_weights, histtype="step", color="k", label="SM")

    #eft_variations = eft.getEFTVariationsLinear()
    minim = 1.
    maxim = 1.
        
    mix_list = [#"cQQ1_cQt1","cQQ1_cQt8",
		    "cQQ1_ctHIm","cQQ1_ctHRe",
		    #"cQQ1_ctt",
		    #"cQQ8_cQQ1","cQQ8_cQt1","cQQ8_cQt8","cQQ8_ctHIm","cQQ8_ctHRe",
		    #"cQQ8_ctt",
		    #"cQt1_cQt8",
		    "cQt1_ctHIm","cQt1_ctHRe",
		    #"cQt1_ctt",
		    "cQt8_ctHIm","cQt8_ctHRe",
		    "ctHRe_ctHIm",
		    #"ctt_cQt8",
		    "ctt_ctHIm","ctt_ctHRe"
		    ]
    for eft_var in mix_list:
        mix_name = "EFT_"+eft_var
        mix_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][mix_name]["Up"])) , nan=0.)
        ax_main.hist(binning[:-1], binning, weights=mix_ratio, histtype="step" , label=eft_var)
        minim = min(minim, np.min(mix_ratio))
        maxim = max(maxim, np.max(mix_ratio))

    #ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=np.ones(len(nominal_content)), yerr=stat_unc_var, ecolor='k', label="stat unc.")
    ax_main.set_xlim(variable.range)
    ax_main.set_ylabel("mixed term")
    modify_yrange_updown(ax_main, (minim, maxim), up_scale=-0.95)
    ax_main.legend(ncol=2)
    ax_main.set_xlabel(variable.axis_label)
    ax_main.text(0.049, 0.77, plotlabel, transform=ax_main.transAxes)

    # fix output name
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_mix.png"))
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_mix.pdf"))
#    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_mix.root"))
    plt.close(fig)

def mix_SM_plot_EFT(variable: Variable, plotdir: str, histograms, process_info: dict, plotlabel: str, processname: str):
    fig, ax_main = fg.create_singleplot()

    binning = generate_binning(variable.range, variable.nbins)
    # first plot nominal, then start adding variations
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    #stat_unc_var = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name]["stat_unc"])) / nominal_content, nan=0.)
    # pretty_name = generate_process_name("SM", info)
    nominal_weights = np.ones(len(nominal_content))
    #ax_main.hist(binning[:-1], binning, weights=nominal_weights, histtype="step", color="k", label="SM")

    #eft_variations = eft.getEFTVariationsLinear()
    minim = 1.
    maxim = 1.
        
    mix_list = [#"cQQ1_cQt1","cQQ1_cQt8",
		    "cQQ1_ctHIm","cQQ1_ctHRe",
		    #"cQQ1_ctt","cQQ8_cQQ1","cQQ8_cQt1","cQQ8_cQt8",
		    #"cQQ8_ctHIm","cQQ8_ctHRe",
		    #"cQQ8_ctt",
		    #"cQt1_cQt8",
		    "cQt1_ctHIm","cQt1_ctHRe",
		    #"cQt1_ctt",
		    "cQt8_ctHIm","cQt8_ctHRe",
		    "ctHRe_ctHIm",
		    #"ctt_cQt8",
		    "ctt_ctHIm","ctt_ctHRe"
		    ]
    for eft_var in mix_list:
        mix_name = "EFT_"+eft_var
        wc1 = eft_var.split("_")[0]
        wc2 = eft_var.split("_")[1]
        lin_name1 = "EFT_" + wc1
        lin_name2 = "EFT_" + wc2
        quad_name1 = "EFT_" + wc1 +"_"+ wc1
        quad_name2 = "EFT_" + wc2 +"_"+ wc2
        mix_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][mix_name]["Up"])) , nan=0.)
        lin1_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][lin_name1]["Up"])) / nominal_content,nan=0.)
        lin2_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][lin_name2]["Up"])) / nominal_content,nan=0.)
        quad1_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][quad_name1]["Up"])) / nominal_content,nan=0.)
        quad2_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][quad_name2]["Up"])) / nominal_content,nan=0.)
        mix_sm = nominal_content + 2*mix_ratio +  lin1_ratio + lin1_ratio + quad1_ratio + quad2_ratio
        ax_main.hist(binning[:-1], binning, weights=mix_sm, histtype="step" , label=eft_var)
        minim = min(minim, np.min(mix_sm))
        maxim = max(maxim, np.max(mix_sm))

    #ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=np.ones(len(nominal_content)), yerr=stat_unc_var, ecolor='k', label="stat unc.")
    ax_main.set_xlim(variable.range)
    ax_main.set_ylabel("SM_lin_quad_mixed term")
    modify_yrange_updown(ax_main, (minim, maxim), up_scale=2)
    ax_main.legend(ncol=2)
    ax_main.set_xlabel(variable.axis_label)
    ax_main.text(0.049, 0.77, plotlabel, transform=ax_main.transAxes)

    # fix output name
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_mix_sm.png"))
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_mix_sm.pdf"))
#    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_mix.root"))
    plt.close(fig)

def ttH_EFT(variable: Variable, plotdir: str, histograms, process_info: dict, plotlabel: str, processname: str):
    fig, ax_main = fg.create_singleplot()

    binning = generate_binning(variable.range, variable.nbins)
    # first plot nominal, then start adding variations
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    stat_unc_var = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name]["stat_unc"])) / nominal_content, nan=0.)
    #pretty_name = generate_process_name("SM", info)
    nominal_weights = np.ones(len(nominal_content))
    #ax_main.hist(binning[:-1], binning, weights=nominal_weights, histtype="step", color="k", label="SM")
    ax_main.hist(binning[:-1], binning, weights=nominal_content , histtype="step", label="sm")
    ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=nominal_content,fmt='none', yerr=stat_unc_var, label="stat unc.")

    
    eft_variations = eft.getEFTVariationsLinear()
    minim = 1.
    maxim = 1.
    for eft_var in eft_variations:
        wc_factor = 1
        if "ctHRe" in eft_var:
            wc_factor = 5
        elif "ctHIm" in eft_var:
            wc_factor = 15
        else :
            continue
        lin_name = "EFT_" + eft_var
        quad_name = lin_name + "_" + eft_var

        current_variation = nominal_content + wc_factor  * np.array(ak.to_numpy(histograms[variable.name][lin_name]["Up"])) + wc_factor * wc_factor * np.array(ak.to_numpy(histograms[variable.name][quad_name]["Up"]))
        quad =  wc_factor * wc_factor *np.array(ak.to_numpy(histograms[variable.name][quad_name]["Up"]))
        pretty_eft_name = eft_var + f" = {wc_factor}"
        ax_main.hist(binning[:-1], binning, weights=current_variation, histtype="step", label=pretty_eft_name)
        ax_main.hist(binning[:-1], binning, weights=quad, histtype="step", label=eft_var + f"_quad = {wc_factor}")
        for thing1 in current_variation : 
           if thing1<=0 :print("sm_lin_quad : ",current_variation)
        for thing2 in quad:
           if thing2<=0 :print("quad : ",quad)

    #eft_variations = eft.getEFTVariationsLinear()
    minim = 1.
    maxim = 1.
        
    mix_list = [#"cQQ1_cQt1","cQQ1_cQt8",
    #"cQQ1_ctHIm","cQQ1_ctHRe",
    #"cQQ1_ctt","cQQ8_cQQ1","cQQ8_cQt1","cQQ8_cQt8",
    #"cQQ8_ctHIm","cQQ8_ctHRe",
    #"cQQ8_ctt",
    #"cQt1_cQt8",
    #"cQt1_ctHIm","cQt1_ctHRe",
    #"cQt1_ctt",
    #"cQt8_ctHIm","cQt8_ctHRe",
    "ctHRe_ctHIm",
    #"ctt_cQt8",
    #"ctt_ctHIm","ctt_ctHRe"
    ]
    for eft_var in mix_list:
        mix_name = "EFT_"+eft_var
        wc1 = eft_var.split("_")[0]
        wc2 = eft_var.split("_")[1]
        lin_name1 = "EFT_" + wc1
        lin_name2 = "EFT_" + wc2
        quad_name1 = "EFT_" + wc1 +"_"+ wc1
        quad_name2 = "EFT_" + wc2 +"_"+ wc2
        mix_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][mix_name]["Up"])) , nan=0.)
        lin1_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][lin_name1]["Up"])) ,nan=0.)
        lin2_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][lin_name2]["Up"])) ,nan=0.)
        quad1_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][quad_name1]["Up"])) ,nan=0.)
        quad2_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][quad_name2]["Up"])),nan=0.)
        mix_sm = nominal_content + 2*mix_ratio +  lin1_ratio + lin1_ratio + quad1_ratio + quad2_ratio
        ax_main.hist(binning[:-1], binning, weights=mix_sm, histtype="step" , label=eft_var)
        minim = min(minim, np.min(mix_sm))
        maxim = max(maxim, np.max(mix_sm))
        for thing3 in mix_sm:
           if thing3<=0 : print("mix : " , mix_sm)

    ax_main.set_xlim(variable.range)
    ax_main.set_ylabel("rate")
    modify_yrange_updown(ax_main, (minim, maxim), up_scale=2)
    ax_main.legend(ncol=1)
    ax_main.set_xlabel(variable.axis_label)
    ax_main.text(0.049, 0.77, plotlabel, transform=ax_main.transAxes)

    # fix output name
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_all.png"))
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_all.pdf"))
    plt.close(fig)

import os
import numpy as np
import awkward as ak

# Make sure these functions/modules are available in your context:
# - generate_binning(range, nbins)
# - eft.getEFTVariationsLinear()

import os
import numpy as np
import awkward as ak

# Ensure these functions or modules are available in your environment:
# - generate_binning(range, nbins)
# - eft.getEFTVariationsLinear()

def main_save_EFT(variable, plotdir, histograms, process_info, plotlabel, processname):
    """
    Computes EFT weight variations relative to the nominal histogram content and saves 
    all relevant arrays to a single compressed .npz file.

    Parameters
    ----------
    variable : object
        Object with attributes:
          - range: tuple (min, max) defining the histogram range.
          - nbins: number of bins.
          - name : key name to access the histogram content in `histograms`.
    plotdir : str
        Directory where the output file will be saved.
    histograms : dict
        Dictionary containing the nominal histogram and EFT variations stored as Awkward arrays.
        Expected keys:
          - histograms[variable.name]["nominal"]
          - histograms[variable.name][<key>]["Up"]
    process_info : dict
        Dictionary with process-related information (not used directly in the computation).
    plotlabel : str
        A label string (unused in saving but kept for interface consistency).
    processname : str
        Process name used in the naming of the output file.

    Returns
    -------
    None
        The computed arrays are saved in a .npz file.
    """
    # Ensure the output directory exists.
    os.makedirs(plotdir, exist_ok=True)
    
    # Create the histogram binning based on the variable's range and number of bins.
    binning = generate_binning(variable.range, variable.nbins)
    
    # Extract nominal histogram content and define nominal weights (array of ones).
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    nominal_weights = np.ones_like(nominal_content)
    
    # Prepare a dictionary to hold all computed arrays.
    results = {
        "binning": binning,
        "nominal": nominal_content,
        "nominal_weights": nominal_weights,
        "plotlabel": plotlabel,  # Stored for reference if needed.
    }
    
    # Retrieve the list of EFT variations.
    eft_variations = eft.getEFTVariationsLinear()
    eft_weights = {}
    eft_value = {}
    
    # Loop through each EFT variation to compute the corresponding weight array.
    for eft_var in eft_variations:
        # Set Wilson coefficient factor for selected EFT variations.
        wc_factor = 10 if ("ctHRe" in eft_var or "ctHIm" in eft_var) else 1
        
        # Construct the keys for the linear and quadratic contributions.
        lin_key = f"EFT_{eft_var}"
        quad_key = f"{lin_key}_{eft_var}"
        
        # Calculate the variation: add the linear and quadratic contributions to the nominal content.
        variation = (
            nominal_content +
            wc_factor * np.array(ak.to_numpy(histograms[variable.name][lin_key]["Up"])) +
            (wc_factor ** 2) * np.array(ak.to_numpy(histograms[variable.name][quad_key]["Up"]))
        )
        # Normalize relative to the nominal content, handling any division by zero issues.
        ratio_variation = np.nan_to_num(variation / nominal_content, nan=1.)
        
        # Save this EFT variation under its name.
        eft_weights[eft_var] = ratio_variation
        eft_value[eft_var] = variation
        sm_value = nominal_content
    
    # Add the EFT variations to the results dictionary.
    results["EFT_variations"] = eft_weights
    results["EFT"] = eft_value
    results["sm"] = sm_value
    
    # Build a file name that embeds the process name and EFT variation names.
    eft_names_str = "_".join(eft_variations)
    filename = f"{processname}_EFT_{eft_names_str}.npz"
    output_filepath = os.path.join(plotdir, filename)
    
    # Save the results to a compressed .npz file.
    np.savez(output_filepath, **results)
    print(f"EFT arrays successfully saved to: {output_filepath}")

if __name__ == "__main__":
    args = parse_arguments()
    np.seterr(divide='ignore', invalid='ignore')

    # load process specifics
    # need a set of processes
    with open(args.processfile, 'r') as f:
        processfile = json.load(f)
        processinfo = processfile["Processes"][args.process]
        subbasedir = processfile["Basedir"].split("/")[-1]

    variables = VariableReader(args.variablefile, args.variable)
    channels = load_channels(args.channelfile)
    storagepath = os.path.join(args.storage, subbasedir)

    outputfolder_base = generate_outputfolder(args.years, args.outputfolder, subbasedir, suffix="_EFT_Variations_v2")

    # first plot nominal, then start adding variations
    # load variables, want to do this for all processes

    # also load channels

    # contrary to plotHistograms: load uncertainties if no args.UseEFT
    # do need a selection somewhere defined for the uncertainties needed/desired
    for channel in channels:
        print(channel)
        if args.channel is not None and channel != args.channel:
            continue

        storagepath_tmp = os.path.join(storagepath, channel)
        systematics = ["nominal", "stat_unc"]

        systematics.extend(eft.getEFTVariationsGroomed())
        outputfolder = os.path.join(outputfolder_base, channel)
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        copy_index_html(outputfolder)

        histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics, args.years[0])
        histograms.load_histograms()

        for _, variable in variables.get_variable_objects().items():
            if not variable.is_channel_relevant(channel):
                continue
#            lin_quad_plot_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)
#            main_plot_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)
#            mix_plot_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)
#            mix_SM_plot_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)
#            ttH_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)
            main_save_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)

        for subchannel in channels[channel].subchannels.keys():
            print(subchannel)
            if not variable.is_channel_relevant(channel + subchannel):
                continue
            storagepath_tmp = os.path.join(storagepath, channel + subchannel)
            
            outputfolder = os.path.join(outputfolder_base, channel, subchannel)
            if not os.path.exists(outputfolder):
                os.makedirs(outputfolder)
            if not os.path.exists(outputfolder):
                 os.makedirs(outputfolder)
            copy_index_html(outputfolder)
            
            # for process, info in processinfo.items():
            histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics, args.years[0])
            histograms.load_histograms()
            for _, variable in variables.get_variable_objects().items():
#                lin_quad_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
#                main_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
#                mix_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
#                mix_SM_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
#                ttH_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
                main_save_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)

