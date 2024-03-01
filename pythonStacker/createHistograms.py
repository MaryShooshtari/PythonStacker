# import numpy as np
# import awkward as ak
import argparse
import sys
import os
import json

import uproot
import awkward as ak
import src
from src.variables.variableReader import VariableReader, Variable
from src.histogramTools import HistogramManager
from submitHistogramCreation import args_add_settingfiles, args_select_specifics
from src.configuration import load_uncertainties, Channel

"""
Script that takes as input a file or set of files, applies cross sections and necessary normalizations (if still needed), and then creates a histogram.
The histogram is taken from settingfiles/Histogramming, so maybe do multiple histograms for the same inputfile? But take range, nbins etc from there as a start.

This is likely code that will be submitted.
Alternatively: create for a single histogram all inputs, so read all files etc.

Must be generic enough: ie. weight variations are easy to support.
More systematic variations is a different question but we can work something out.
"""


def parse_arguments() -> argparse.Namespace:
    """
    For now directly takne from create_histograms in interpretations.
    """
    parser = argparse.ArgumentParser(description='Process command line arguments.')

    # add file arguments
    args_add_settingfiles(parser)
    args_select_specifics(parser)

    parser.add_argument("--storage", dest="storage", type=str,
                        default="Intermediate", help="Path at which the \
                        histograms are stored")

    # Parse arguments
    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)
        exit(1)

    if args.process is None:
        print("Need to specify a process to produce. Exiting...")
        exit(1)

    if args.channel is None:
        print("Need to specify a channel to produce. Exiting...")
        exit(1)

    return args


def prepare_histogram(data, wgts, variable: Variable):
    hist_content, binning, hist_unc = src.histogram_w_unc_flow(data, range=variable.range, wgts=wgts, nbins=variable.nbins)
    return hist_content, binning, hist_unc


# def analysis_histograms(processinfo: dict, processname: str, variables: VariableReader, storagepath: str, selection: str):
#     basedir = processinfo["Basedir"]
#     filelist = processinfo["Processes"][processname]["files"]
#
#     output_histograms = HistogramManager(storagepath, processname, variables)
#     output_content = dict()  # variable: content
#     output_unc = dict()
#     for file in filelist:
#         current_tree = uproot.open(os.path.join(basedir, file))["test"]
#         weights = ak.to_numpy(current_tree.arrays("weights", cut=selection, aliases={"weights": "nominalWeight"}).weights)
#         input()
#         # just first do central, worry about systematics later
#         print(f"number of variables: {variables.number_of_variables()}")
#
#         for variable in variables.get_variables():
#             hist_content, hist_unc = create_histogram(variables.get_properties(variable), current_tree, weights, selection)
#
#             if variable in output_content:
#                 output_content[variable]["nominal"] += hist_content
#                 output_unc[variable]["nominal"] += hist_unc
#             else:
#                 output_content[variable] = dict()
#                 output_unc[variable] = dict()
#                 output_content[variable]["nominal"] = hist_content
#                 output_unc[variable]["nominal"] = hist_unc
#
#         # load weights we need based on either systematics or diffeerent. Write function to get right weights
#         # load central
#         # Data stays the same for each variable in weight variations.
#         # So load data once with correct cut, also need to specify the cut somewhere...
#         # But then load: the data and the central weights, then for each systematic load the weights, but againthey should be stationary.
#         # maybe store the systematic variation weights also somewhere?
#
#     # save output
#     output_histograms.save_all_histograms(output_content)


def get_histogram_data(variable: Variable, tree, channel: Channel):
    method = variable.get_method()
    data = method(tree, variable.branch_name, channel.selection)
    return data


def create_histogram(variable: Variable, data, weights):
    '''
    TODO
    '''
    # get opts somewhere?
    # method = variable.get_method()
    # data = method(tree, variable.branch_name, selection)

    # to histogram:
    hist_content, _, hist_unc = prepare_histogram(data, weights, variable)
    return hist_content, hist_unc


def weight_variations(variable: Variable, tree, channel: Channel, weights, systematics: dict):
    # load data:
    data = get_histogram_data(variable, tree, channel)

    for systematic in systematics:
        hist_content, _, hist_unc = prepare_histogram(data, weights, variable)

        output_histograms[variable][systematic] += hist_content
        if systematic == "nominal":
            output_histograms[variable]["stat_unc"] += hist_unc


if __name__ == "__main__":
    # parse arguments
    args = parse_arguments()

    # make sure we have correct storagepath
    storagepath = args.storage

    # initialize variable class:
    variables = VariableReader(args.variablefile, args.variable)

    # prob do same with systematics
    if args.systematic == "weight" or args.systematic == "shape":
        systematics: dict = load_uncertainties(args.systematicsfile, typefilter=args.systematic, allowflat=False)
    elif args.systematic is not None:
        systematics: dict = load_uncertainties(args.systematicsfile, namefilter=args.systematic, allowflat=False)
    else:
        systematics: dict = dict()

    if args.systematic != "shape":
        systematics["nominal"] = 0
        systematics["stat_unc"] = 0

    # load process list:
    with open(args.processfile, 'r') as f:
        processfile = json.load(f)
        processinfo = processfile["Processes"][args.process]
        basedir = processfile["Basedir"]

    # prepare channels:
    with open(args.channelfile, 'r') as f:
        # A single job per channel, so we only need to load the current one.
        # This ensures memory usage is refreshed at the end.
        channelinfo = json.load(f)
        channel = Channel(channelinfo[args.channel], channelinfo)
        if channel.isSubchannel:
            print("Current channel is a subchannel. Nothing to be done. Exiting...")
            exit(0)

    # Now also do this for systematic variations? How though?
    # Maybe add an argument for the systematic variations to produce, can define these somewhere.
    # Should be weightvariations here? Maybe also systematic variations or something? idk let's see first for weight variations

    # update storagepath to include channel
    storagepath = os.path.join(storagepath, args.channel)
    output_histograms = HistogramManager(storagepath, args.process, variables, list(systematics.keys()))

    # TODO: get files based on process names -> processmanager can return this, depending on the sys unc?
    for filebase in processinfo["files"]:
        # filebase should not include a suffix
        # generate basepath with correct folder, folder has timestamp now
        # then:
        # if args.systematic == "shape":
        # first loop systematics, prob just one at a time -> or create filename with "systematic" tag
        # then loop variables
        # if args.systematic == "weight":
        # loop variables
        filename = os.path.join(basedir, filebase)
        if args.systematic == "weight" or args.systematic is None:
            pass
            # filebase += "_base.root"

        current_tree = uproot.open(filename)["test"]

        # TODO: use masks somewhere
        # generates masks for subchannels
        subchannelmasks, subchannelnames = channel.produce_masks(current_tree)

        # no structure yet to load systematics weights. Weightmanager? Or move to systematics loop?
        # Does imply more overhead in reading, can try here, see what memory effect it has, otherwise move to systematics.
        weights = ak.to_numpy(current_tree.arrays(["weights"], cut=channel.selection, aliases={"weights": "nominalWeight"}).weights)

        for variable in variables:
            # load data:
            data = get_histogram_data(variable, current_tree, channel)

            for systematic in systematics:
                hist_content, _, hist_unc = prepare_histogram(data, weights, variable)

                output_histograms[variable.name][systematic] += hist_content
                if systematic == "nominal":
                    output_histograms[variable.name]["stat_unc"] += hist_unc

    output_histograms.save_histograms()

    exit(0)