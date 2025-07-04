#!/usr/bin/env python3
import os
import argparse
import numpy as np
import ROOT

# Set ROOT to batch mode to avoid opening interactive windows
ROOT.gROOT.SetBatch(True)

# --- Internal library: fixed colors and pretty labels for EFT coefficients --- #
DEFAULT_COEFFICIENT_COLORS = {
    "ctt":  ROOT.TColor.GetColor("#5790fc"),  
    "cQQ1": ROOT.TColor.GetColor("#f89c20"),  
    #"cQQ8": ROOT.TColor.GetColor("#e42536"),  
    "cQt1": ROOT.TColor.GetColor("#964a8b"),  
    "cQt8": ROOT.TColor.GetColor("#9c9ca1"),  
    "ctHRe": ROOT.TColor.GetColor("#e42536"), 
    "ctHIm": ROOT.TColor.GetColor("#7a21dd"), 
		}

DEFAULT_COEFFICIENT_PRETTY = {
    "ctt": "#it{c}_{tt} = 1",
    "cQQ1": "#it{c}_{QQ}^{(1)} = 1",
    #"cQQ8": "#it{c}_{QQ}^{(8)} = 1",
    "cQt1": "#it{c}_{Qt}^{(1)} = 1",
    "cQt8": "#it{c}_{Qt}^{(8)} = 1",
    "ctHRe": "#it{c}_{tH}^{Re} = 10",
    "ctHIm": "#it{c}_{tH}^{Im} = 10",
}

FALLBACK_COLORS = [ROOT.kOrange, ROOT.kCyan, ROOT.kViolet, ROOT.kGray, ROOT.kSpring]

# --- Helper function to create a histogram from binning and values --- #
def create_histogram(hist_name, binning, values):
    """
    Creates a TH1D histogram with the provided variable binning and fills in its bin contents.

    Parameters:
      hist_name : Name of the histogram.
      binning   : Iterable of bin edges (length = nbins+1).
      values    : Iterable of bin contents (length = nbins).

    Returns:
      A ROOT.TH1D histogram.
    """
    nbins = len(binning) - 1
    # Convert binning to a numpy array of floats (needed for variable-width histograms)
    arr = np.array(binning, dtype=float)
    h = ROOT.TH1D(hist_name, "", nbins, arr)
    for i in range(1, nbins + 1):
        h.SetBinContent(i, values[i - 1])
    h.SetLineWidth(2)
    return h

# --- Main function --- #
def main():
    parser = argparse.ArgumentParser(
        description="Plot EFT variations (as histograms) from an .npz file using ROOT. "
                    "Output is saved in the same directory as the input file."
    )
    parser.add_argument("--indir", required=True,
                        help="Directory containing input .npz files.")
    parser.add_argument("--files", nargs='+', required=True,
                        help="List of input .npz file names (relative to --indir).")
    parser.add_argument("--output", "-o", type=str, required=True,
                        help="Name of the output plot file (saved in the same directory as the input files).") 
    args = parser.parse_args()

    # Determine the directory of the input file and build the output file path.
    input_dir = os.path.abspath(args.indir)
    file_paths = [os.path.join(input_dir, f+"_EFT_EFT_cQQ8_cQQ1_cQt1_ctt_cQt8_ctHRe_ctHIm.npz") for f in args.files]
    output_filepath = os.path.join(input_dir, args.output)

    # Initialize accumulators for EFT variations.
    sum_eft_variations = {}  # key: operator, value: summed array
    binning = None
    plotlabel = ""
    n_files = 0


    # Initialize accumulators.
    sum_sm = None              # to accumulate SM arrays over files
    sum_eft = {}               # dictionary, key: operator, value: accumulated EFT arrays


    # Loop over input files and accumulate EFT variations.
    for file_path in file_paths:
        try:
            data = np.load(file_path, allow_pickle=True)
        except Exception as e:
            print("Error loading file", file_path, ":", e)
            continue
        
        # For the first file, store binning and plotlabel.
        if binning is None:
            binning = data["binning"]
            # Use plotlabel if available.
            if "plotlabel" in data:
                plotlabel = str(data["plotlabel"])
            else:
                plotlabel = "EFT Variations"
        else:
            # Optionally: check if binning is consistent.
            if not np.array_equal(binning, data["binning"]):
                print("Warning: Binning does not match in", file_path)
        
        n_files += 1

        # Accumulate the SM arrays.
        file_sm = np.array(data["sm"])  # use the key "sm"
        if sum_sm is None:
            sum_sm = file_sm.copy()
        else:
            sum_sm += file_sm
    
        # Retrieve the EFT dictionary from the file.
        # Expected structure: data["EFT"] is a dict mapping operator names to arrays.
        eft_dict = data["EFT"].item()
        data.close()
    
        order = list(DEFAULT_COEFFICIENT_COLORS.keys())
        
        for key in order:
            if key in eft_dict:
                # pop it out and re‑insert → moves key to the end, in the order you loop
                eft_dict[key] = eft_dict.pop(key)
        # For each operator, accumulate its EFT array.
        for op, arr in eft_dict.items():
            arr_np = np.array(arr)  # convert to NumPy array
            if op in sum_eft:
                sum_eft[op] += arr_np
            else:
                sum_eft[op] = arr_np.copy()
    
        if n_files == 0:
            print("No files loaded successfully.")
            return
    
    # After processing all files, compute the overall variation for each operator:
    final_variation = {}
    for op in sum_eft:
        # Compute the ratio: total EFT divided by total SM, elementwise.
        final_variation[op] = sum_eft[op] / sum_sm
        print(op ,final_variation[op])
    
    print("Summed EFT variations for {} files.".format(n_files))


    # --- Begin Plotting with ROOT --- #
    # --- Keep references to objects to prevent garbage collection ---
    keepers = []    # list for histograms, TLatex, legend, etc.

    # --- Set up the ROOT canvas ---
    c = ROOT.TCanvas("c", "EFT Variations", 600, 600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.15)

    # --- Draw the SM baseline histogram ---
    nbins = len(binning) - 1
    sm_values = [1.0] * nbins
    h_SM = create_histogram("h_SM", binning, sm_values)
    h_SM.SetLineColor(ROOT.kBlack)
    h_SM.GetXaxis().SetTitle("H_{T} [GeV]")
    h_SM.GetXaxis().SetTitleSize(0.048)
    h_SM.GetYaxis().SetTitleSize(0.048)
    h_SM.GetXaxis().SetLabelSize(0.045)
    h_SM.GetYaxis().SetLabelSize(0.045)
    h_SM.GetYaxis().SetTitle("SM + EFT / SM")
    h_SM.SetStats(0)  # disable stats box
    h_SM.SetMinimum(0.4)
    h_SM.SetMaximum(3.0)
    h_SM.Draw("hist")
    keepers.append(h_SM)

    # --- Create and configure the legend ---
    legend = ROOT.TLegend(0.45, 0.7, 0.95, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.05)
    try:
        legend.SetNColumns(2)
    except AttributeError:
        pass  # older ROOT versions may not support this
    #legend.AddEntry(h_SM, "SM", "l")
    keepers.append(legend)

    fallback_index = 0
    # Loop over each EFT operator and draw its summed histogram.
    for op, summed_arr in final_variation.items():
        if op=="cQQ8": continue
        h_name = "h_" + op
        h_eft = create_histogram(h_name, binning, summed_arr)
        # Use the color from our palette if available.
        if op in DEFAULT_COEFFICIENT_COLORS:
            h_eft.SetLineColor(DEFAULT_COEFFICIENT_COLORS[op])
        else:
            h_eft.SetLineColor(FALLBACK_COLOR)
        h_eft.Draw("hist same")
        keepers.append(h_eft)
        # Use pretty label if available.
        label = DEFAULT_COEFFICIENT_PRETTY.get(op, op)
        legend.AddEntry(h_eft, label, "l")

    legend.Draw("same")

    # --- Draw CMS and Simulation labels (top-left) ---
    cms_label = ROOT.TLatex()
    cms_label.SetNDC()
    cms_label.SetTextFont(61)  # Bold for "CMS"
    cms_label.SetTextSize(0.09)
    cms_label.DrawLatex(0.20, 0.8, "CMS")
    keepers.append(cms_label)

    cms_sim_label = ROOT.TLatex()
    cms_sim_label.SetNDC()
    cms_sim_label.SetTextFont(52)  # Lighter italic for "Simulation"
    cms_sim_label.SetTextSize(0.045)
    cms_sim_label.DrawLatex(0.20, 0.74, "Simulation")
    keepers.append(cms_sim_label)

    # --- Draw additional label (e.g. plot category) in top-right ---
    extra_label = ROOT.TLatex()
    extra_label.SetNDC()
    extra_label.SetTextFont(42)
    extra_label.SetTextSize(0.04)
    extra_label.DrawLatex(0.7, 0.96, r'138 fb^{-1}(13TeV)')
    keepers.append(extra_label)

    c.Update()
    for ext in [".pdf",".png"]:
        c.Print(output_filepath+ext)
    print("Plot successfully saved to:", output_filepath)

if __name__ == '__main__':
    main()

