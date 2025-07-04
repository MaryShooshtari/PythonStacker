#!/usr/bin/env python3
import ROOT
import json
import argparse
import glob
import uproot

parser = argparse.ArgumentParser(description="Extract limit curves from ROOT files to JSON.")
parser.add_argument("--fixed", required=True, help="ROOT file with fixed scan")
parser.add_argument("--profiled", required=True, help="ROOT file with profiled scan")
parser.add_argument("--WC", required=True, help="Name of the Wilson Coefficient (WC) to extract")
parser.add_argument("--output", default="limits.json", help="Output JSON file")
args = parser.parse_args()

# Function to extract (WC, deltaNLL) pairs from a ROOT file
#def extract_points(filename, wc_name):
#    file = ROOT.TFile.Open(filename)
#    tree = file.Get("limit")
#    points = []
#    for entry in tree:
#        wc_value = getattr(entry, wc_name, None)
#        delta_nll = getattr(entry, "deltaNLL", None)
#        if wc_value is not None and delta_nll is not None:
#            points.append([float(wc_value), float(2 * delta_nll)])
#    return points

def extract_points(filename, wc_name):
    file = ROOT.TFile.Open(filename)
    tree = file.Get("limit")
    wc_to_min_nll = {}
    for entry in tree:
        wc_value = getattr(entry, wc_name, None)
        delta_nll = getattr(entry, "deltaNLL", None)
        if wc_value is not None and delta_nll is not None:
            wc_value = float(wc_value)
            delta_nll = float(2 * delta_nll)
            if wc_value not in wc_to_min_nll or delta_nll < wc_to_min_nll[wc_value]:
                wc_to_min_nll[wc_value] = delta_nll
    # Convert the dictionary to a list of [wc_value, deltaNLL] pairs
    points = [[wc, nll] for wc, nll in sorted(wc_to_min_nll.items())]
    return points

def extract_points_many_files(file_pattern: str, wc_name: str) -> list[list[float]]:
    """
    Scan all ROOT files matching `file_pattern`, read the TTree "limit",
    and return a sorted list of [wc_value, deltaNLL] pairs.

    Parameters
    ----------
    file_pattern : str
        Glob pattern for the ROOT files, e.g. "higgsCombine_*_profiled.*.root"
    wc_name : str
        Name of the branch in the TTree corresponding to the Wilson coefficient.

    Returns
    -------
    points : list of [float, float]
        Sorted by the WC value (ascending).
    """
    # find & sort files
    files = sorted(glob.glob(file_pattern))
    points = []
    for fname in files:
        with uproot.open(fname) as f:
            tree = f["limit"]
            # read arrays
            wc_vals  = tree[wc_name].array(library="np")
            nll_vals = tree["deltaNLL"].array(library="np")
            # collect pairs
            for x, y in zip(wc_vals, nll_vals):
                points.append([float(x), float(y)])
    # determine minimum deltaNLL
    min_nll = min(y for _, y in points)
    # normalize and sort
    normalized = [[x, 2*(y - min_nll)] for x, y in points]
    normalized.sort(key=lambda p: p[0])
    return normalized
#    # sort by the WC value
#    points.sort(key=lambda p: p[0])
#    return points

# Extract points from both files
fixed_points = extract_points(args.fixed, args.WC)
profiled_points = extract_points(args.profiled, args.WC)
#profiled_points = extract_points_many_files(args.profiled, args.WC)


# Sort by WC value
fixed_points.sort(key=lambda x: x[0])
profiled_points.sort(key=lambda x: x[0])

# Create final JSON structure
out_data = {
    "tttt_ttt_ttH_exp": {
        "WC": args.WC,
        "fixed": fixed_points,
        "profiled": profiled_points
    }
}

# Write to file
with open(args.output, "w") as f:
    json.dump(out_data, f, indent=2)

print(f"Saved JSON to {args.output}")

