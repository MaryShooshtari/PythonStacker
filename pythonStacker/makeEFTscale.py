import os
import json
import pandas as pd
import argparse
from collections import defaultdict
import math

def generate_multiprocess_json(coefficients, input_location, channel_map, year_list, processes, output_path):
    mix_list = [
        "cQQ1_cQt1", "cQQ1_cQt8", "cQQ1_ctHIm", "cQQ1_ctHRe", "cQQ1_ctt", "cQQ8_cQQ1", "cQQ8_cQt1",
        "cQQ8_cQt8", "cQQ8_ctHIm", "cQQ8_ctHRe", "cQQ8_ctt", "cQt1_cQt8", "cQt1_ctHIm", "cQt1_ctHRe",
        "cQt1_ctt", "cQt8_ctHIm", "cQt8_ctHRe", "ctHRe_ctHIm", "ctt_cQt8", "ctt_ctHIm", "ctt_ctHRe"
    ]

    output_list = []
    for year in year_list:
        for channel, (variable ,prettyname) in channel_map.items():
            for process in processes:
                sm_path = os.path.join(input_location, channel, variable, f"{process}_EFT_{year}_EFT_ctt.parquet")
                #print(f"[INFO] Looking for SM file: {sm_path}")
                if not os.path.exists(sm_path):
                    print(f"SM file not found: {sm_path}")
                    continue
    
                sm_df = pd.read_parquet(sm_path)
                nbins = len(sm_df['Up'][0])
                print(f"[INFO] Found SM file for {channel}/{process} with {nbins} bins")
    
                # Build list of dicts per bin mapping monomial term -> value
                process_bins = []
                for i in range(nbins):
                    sm_val = sm_df['Down'][0][i]
                    process_bins.append({"r_SM*r_SM": sm_val})
    
                # Interference terms: r_SM * coef
                for coef in coefficients:
                    int_path = os.path.join(input_location, channel, variable, f"{process}_EFT_{year}_EFT_{coef}.parquet")
                    #print(f"[INFO] Checking interference file: {int_path}")
                    if os.path.exists(int_path):
                        int_df = pd.read_parquet(int_path)
                        for i in range(nbins):
                            val = int_df['Up'][0][i]
                            sm_down_val = sm_df['Down'][0][i]
                            if val != sm_down_val:
                                term = f"r_SM*{coef}"
                                #print(f"[INFO] Bin {i}: adding {term} = {val}")
				# divide by 2
                                process_bins[i][term] = val/2
                    else:
                        print(f"[WARNING] Interference file not found: {int_path}")
    
                # Quadratic and cross terms
                for coef1 in coefficients:
                    # Quadratic
                    quad_path = os.path.join(input_location, channel, variable, f"{process}_EFT_{year}_EFT_{coef1}_{coef1}.parquet")
                    #print(f"[INFO] Checking quadratic file: {quad_path}")
                    if os.path.exists(quad_path):
                        quad_df = pd.read_parquet(quad_path)
                        for i in range(nbins):
                            val = quad_df['Up'][0][i]
                            sm_down_val = sm_df['Down'][0][i]
                            if val != sm_down_val:
                                term = f"{coef1}*{coef1}"
                                #print(f"[INFO] Bin {i}: adding {term} = {val}")
                                process_bins[i][term] = val
                    else:
                        print(f"[WARNING] Quadratic file not found: {quad_path}")
    
                    # 3b) Cross terms: look up exactly in mix_list
                    for coef2 in coefficients:
                        if coef1 >= coef2:
                            continue
                        key1 = f"{coef1}_{coef2}"
                        key2 = f"{coef2}_{coef1}"
                        if key1 in mix_list:
                            mix_key = key1
                        elif key2 in mix_list:
                            mix_key = key2
                        else:
                            # not a valid cross‐term according to mix_list
                            continue
    
                        cross_path = os.path.join(
                            input_location, channel, variable,
                            f"{process}_EFT_{year}_EFT_{mix_key}.parquet"
                        )
                        #print(f"[INFO] Checking cross file: {cross_path}")
                        if os.path.exists(cross_path):
                            cross_df = pd.read_parquet(cross_path)
                            for i in range(nbins):
                                val = cross_df["Up"][0][i]
                                sm_down_val = sm_df["Down"][0][i]
                                if val != sm_down_val:
                                    # ALWAYS name the term exactly as mix_key.replace("_","*")
                                    term = mix_key.replace("_", "*")
                                    #print(f"[INFO] Bin {i}: adding {term} = {val}")
				    #add 2* Mix term (or divide??) https://indico.cern.ch/event/1349900/contributions/5682809/attachments/2762295/4810682/EFTForum-interference.pdf
                                    process_bins[i][term] = val*2/2
                        else:
                            print(f"[WARNING] Cross file not found: {cross_path}")
    
                #
                # 4) Build a fixed monomial_order that matches mix_list ordering:
                #
                monomial_order = ["r_SM*r_SM"]
#                for coef in coefficients:
                for i, ci in enumerate(coefficients):
                    monomial_order.append(f"r_SM*{ci}")
                    print(f"Found key: SM*{ci}")
                    for j, cj in enumerate(coefficients):
                        if j > i:
                            continue  # Skip duplicates like cj_ci when ci_cj was already checked 
                        if ci == cj:
                            #print(f"Identical coefficients: {ci}")
                            key = f"{ci}*{cj}"
                        else:
                            key = next((k for k in (f"{ci}_{cj}", f"{cj}_{ci}") if k in mix_list), None)
                        if key:
                            print(f"Found mix key: {key}")
                            monomial_order.append(key.replace("_", "*"))
                
    
    
                # Parameter names: ['cSM'] + coefficients
                operator_syntax = {
                          "cSM":   "cSM[1]",
                          "ctt":   "ctt[0,-3,3]",
                          "cQQ1":   "cQQ1[0,-3,3]",
                          "cQt1":   "cQt1[0,-5,5]",
                          "cQt8":   "cQt8[0,-6,6]",
                          "ctHRe":   "ctHRe[0,-20,40]",
                          "ctHIm":   "ctHIm[0,-40,40]",
                          }
                bare_ops = ["cSM"] + coefficients
                parameter_names = [operator_syntax[op] for op in bare_ops if op in operator_syntax]
    
                # Build scaling: each bin -> list of monomial values in monomial_order
                scaling = []
                for bin_dict in process_bins:
                    sm_val = bin_dict.get("r_SM*r_SM", 1.0)
                    row = [1.0]
                    for m in monomial_order[1:]:
                        val = bin_dict.get(m, 0.0)
                        row.append(val / sm_val)
#                    row    = [ ((bin_dict.get(m, 0.0)+sm_val )/ sm_val) for m in monomial_order[1:] ]
                    scaling.append(row)
    
                has_bad = any(
                    (val is None) or (val == 0) or (isinstance(val, float) and math.isnan(val))
                    for row in scaling
                    for val in row
                ) 
                if has_bad:
                    print("There is at least one entry that is None, 0, or NaN.")
                    #print(scaling)
                    continue

                pro_names = {
                              "TTTT" : "tttt",
                              "TTTJ" : "tttj",
                              "TTTW" : "tttW",
                              "TTH" : "ttH",
				}
                output_list.append({
                    "channel": "DC_"+year+"_"+prettyname,
                    "process": pro_names[process],
                    "parameters": parameter_names,
                    "scaling": scaling
                })
                
                
    print(f"[INFO] Writing JSON to {output_path}")
    with open(output_path, 'w') as f:
        json.dump(output_list, f, indent=2)
    print(f"[DONE] JSON saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate multiprocess JSON for Combine EFT model using a channel config JSON.")
    parser.add_argument("--config", required=True, help="Path to JSON config containing channelcontent")
    parser.add_argument("--coefficients", nargs="+", required=True, help="List of EFT coefficients (e.g., ctt cQt8)")
    parser.add_argument("--input-location", required=True, help="Base input directory containing the parquet files")
    parser.add_argument("--year", nargs="+", required=True, help="Data-taking year (e.g., 2018)")
    parser.add_argument("--processes", nargs="+", required=True, help="List of physics processes (e.g., TTTT TTZ)")
    parser.add_argument("--output", required=True, help="Path to output JSON file")

    args = parser.parse_args()

    # Read channel config JSON
    with open(args.config, 'r') as cf:
        cfg = json.load(cf)
    channel_map = {ch: (info["variable"], info["prettyname"]) for ch, info in cfg.get('channelcontent', {}).items()}

    generate_multiprocess_json(
        coefficients=args.coefficients,
        input_location=args.input_location,
        channel_map=channel_map,
        year_list=args.year,
        processes=args.processes,
        output_path=args.output
    )

if __name__ == "__main__":
    main()

