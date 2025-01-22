"""
For files with eft weights, prepare and store the weights per event class, since it takes a lot of computing time.
"""

import plugins.eft as eft

import awkward as ak
import argparse
import os

import src


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-f", "--file", dest="inputfile", action="store", required=True,
                        help="Inputfile for which variations must be stored")
    parser.add_argument("-e", "--eventClass", dest="eventclass", action="store", default=11,
                        help="which event class to check")
    parser.add_argument("-p", "--process", dest="process", action="store", default="TTTT_EFT")
    parser.add_argument("--storage", dest="storage", type=str,
                        default="Intermediate", help="Path at which the histograms are stored")

    opts, opts_unknown = parser.parse_known_args()
    return opts


def get_eftvariations_filename(storage: str, filename: str, eventclass) -> str:
    basefilename = filename.split("/")[-1].split(".")[0].split("_localSub")[0]
    if not os.path.exists(os.path.join(storage, "eftWeights")):
        os.makedirs(os.path.join(storage, "eftWeights"))
    return os.path.join(storage, "eftWeights", f"{basefilename}_evClass_{eventclass}_wgts.parquet")


def reweight_and_write(reweighter, eventclass, tree, storage, filename, process):
    arrs = tree.arrays(["eftVariationsNorm"], "eventClass==" + str(eventclass), aliases={"eftVariationsNorm": "eftVariations/eftVariations[:,0]"})
    # should become a record rather than a single array, using the naming scheme I devised in the init

    import numpy as np
    # var = np.array([reweighter.transform_weights(entry[1:]) for entry in arrs.eftVariationsNorm])
    if process == "TTTT_EFT" or process == "TTH_EFT":
        var = np.transpose(np.array(reweighter.transform_weights_parallel(ak.to_numpy(arrs.eftVariationsNorm[:, 1:]))))
        #print("eftVariationsNorm[:, 1:] ", ak.to_numpy(arrs.eftVariationsNorm[:, 1:]))
        print("var not T",np.array(reweighter.transform_weights_parallel(ak.to_numpy(arrs.eftVariationsNorm[:, 1:]))))
        #print("var", var)
    else:
        var = np.transpose(np.array(reweighter.transform_weights_parallel(ak.to_numpy(arrs.eftVariationsNorm))))

    if len(var) == 0:
        return

    # postprocess: get names, then make dict, then make record
    final_prerecord = dict()
    eft_names = ["Central"] + eft.getEFTVariationsGroomed()
    for i, name in enumerate(eft_names):
        final_prerecord[name] = var[:, i]
    outputfile = get_eftvariations_filename(storage, filename, eventclass)
    ak.to_parquet(ak.Record(final_prerecord), outputfile)
    return


if __name__ == "__main__":
    args = parse_arguments()

    reweighter = eft.eft_reweighter()
    # tree = uproot.open(args.inputfile)

    tree = src.get_tree_from_file(args.inputfile, args.process)
    if (int(args.eventclass) == -1):
        for i in range(15):
            reweight_and_write(reweighter, i, tree, args.storage, args.inputfile, args.process)
    else:
        print(f"running weights for class {args.eventclass}")
        reweight_and_write(reweighter, int(args.eventclass), tree, args.storage, args.inputfile, args.process)

    print("Success!")
