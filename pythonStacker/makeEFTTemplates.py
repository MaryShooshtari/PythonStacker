'''
eft_yield_builder.py

Reads SM, interference, quadratic, and cross-term Parquet files per channel/process/bin,
builds the EFT prediction histogram given operator values, optionally plots yields in histogram step style,
and prints the full analytic formula per bin.

Interference and quadratic extraction assume each Up file directly holds the monomial:
  Up_α = I_α,  Up_{αβ} = Q_{αβ}.

**Formula** (per bin i):
    N_i = N_SM,i
          + ∑_α (c_α/Λ²)·I_α,i
          + ∑_{α≤β} (c_α c_β/Λ⁴)·Q_{αβ},i

**Usage**:
    python eft_yield_builder.py \
        --coeff ctt=1 cQQ1=0 ... \
        --scale 1.0 \
        --input-dir /path/to/parquets \
        --config channels.json \
        --year 2018 \
        --processes TTTT TTZ \
        --output yields.json [--plot]
'''
import os,shutil
import json
import argparse
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt


def load_coeff_dict(coeff_list):
    coeffs = {}
    print("[INFO] Parsing coefficients:")
    for kv in coeff_list:
        if '=' not in kv:
            raise ValueError(f"Coefficient spec '{kv}' must be name=value")
        name, val = kv.split('=', 1)
        coeffs[name] = float(val)
        print(f"  - {name} = {coeffs[name]}")
    return coeffs


def compute_bin_yield(sm_val, int_vals, quad_vals, coeffs, scale):
    # Build symbolic terms list for formula printout
    terms = ["SM"]
    contributions = []

    # Linear interference terms
    for a, I in int_vals.items():
        c = coeffs.get(a, 0.0)
        contrib = (c / scale**2) * I
        contributions.append(contrib)
        terms.append(f"{c}*I_{a}")
        #print(f"      [DEBUG] Linear term: c[{a}]={c}, I[{a}]={I} => +{contrib}")

    # Quadratic and cross terms
    for (a, b), Q in quad_vals.items():
        ca = coeffs.get(a, 0.0)
        cb = coeffs.get(b, 0.0)
        contrib = (ca * cb / scale**4) * Q
        contributions.append(contrib)
        terms.append(f"{ca}*{cb}*Q_{a}{b}")
        #print(f"      [DEBUG] Quad term: c[{a}]={ca}, c[{b}]={cb}, Q[{a},{b}]={Q} => +{contrib}")

    # Total yield
    y = sm_val + sum(contributions)

    # Print full symbolic formula
    formula = ' + '.join(terms)
    #print(f"    [FORMULA] N = {formula}")
    #print(f"    [DEBUG] Final N = {y}")
    return y


def main():
    parser = argparse.ArgumentParser(description="Generate yields + optional histogram-step plots with diagnostics")
    parser.add_argument('--coeff', nargs='+', required=True, help='List of c_alpha=value')
    parser.add_argument('--scale', type=float, default=1.0, help='NP scale Λ')
    parser.add_argument('--input-dir', required=True, help='Base directory with parquet files')
    parser.add_argument('--config', required=True, help='JSON with channelcontent')
    parser.add_argument('--year', required=True, nargs='+', help='Year(s)')
    parser.add_argument('--processes', required=True, nargs='+', help='Processes')
    parser.add_argument('--output', required=True, help='Output JSON for yields')
    parser.add_argument('--plot', action='store_true', help='Plot yields as histogram steps + SM')
    args = parser.parse_args()

    #print(f"[INFO] Loading config from {args.config}")
    coeffs = load_coeff_dict(args.coeff)
    # build a suffix string like "ctt1_cQQ11_ctHRe1"
    coeff_str = "_".join(
        f"{name}{int(val) if float(val).is_integer() else val}"
        for name, val in sorted(coeffs.items())
    )
    with open(args.config) as f:
        cfg = json.load(f)
    channel_map = {k: (v['variable'], v['prettyname']) for k, v in cfg.get('channelcontent', {}).items()}

    results = []
    # key = (channel, var, process)
    plot_data = defaultdict(list)  # yields per bin
    plot_sm   = defaultdict(list)  # SM per bin

    for year in args.year:
        for ch_key, (var, pretty) in channel_map.items():
            for proc in args.processes:
                proc_key = proc.lower()
                #print(f"[INFO] Processing Year={year}, Channel={ch_key}, Var={var}, Process={proc}")
                base = os.path.join(args.input_dir, ch_key, var)

                # Load SM base
                sm_path = os.path.join(base, f"{proc}_EFT_{year}_EFT_ctt.parquet")
                if not os.path.exists(sm_path):
                    #print(f"  [WARN] SM file missing: {sm_path}")
                    continue
                sm_df = pd.read_parquet(sm_path)
                nbins = len(sm_df['Down'][0])

                # Initialize plot lists
                key = (year, ch_key, var, proc_key)
                if args.plot:
                    plot_data[key] = []
                    plot_sm[key]   = []

                for i in range(nbins):
                    #print(f"  [INFO] Bin {i}:")
                    sm_val = sm_df['Down'][0][i]
		    #print(f"    [DEBUG] SM value = {sm_val}")
                    if args.plot:
                        plot_sm[key].append(sm_val)

                    # Interference I_alpha
                    int_vals = {}
                    for a in coeffs:
                        ip = os.path.join(base, f"{proc}_EFT_{year}_EFT_{a}.parquet")
                        if os.path.exists(ip):
                            df = pd.read_parquet(ip)
                            I = df['Up'][0][i]
                            int_vals[a] = I
                            #print(f"      [DEBUG] Loaded I[{a}] = {I}")

                    # Quadratic/cross Q_alpha_beta
                    quad_vals = {}
                    names = list(coeffs)
                    for ii, a in enumerate(names):
                        for b in names[ii:]:
                            qp = os.path.join(base, f"{proc}_EFT_{year}_EFT_{a}_{b}.parquet")
                            if os.path.exists(qp):
                                df = pd.read_parquet(qp)
                                Q = df['Up'][0][i]
                                quad_vals[(a, b)] = Q
                                #print(f"      [DEBUG] Loaded Q[{a},{b}] = {Q}")

                    N = compute_bin_yield(sm_val, int_vals, quad_vals, coeffs, args.scale)
                    results.append({
                        'channel': f"DC_{year}_{pretty}",
                        'process': proc_key,
                        'bin': i,
                        'yield': N
                    })
                    if args.plot:
                        plot_data[key].append(N)

    # Write JSON
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"[INFO] Wrote yields for {len(results)} bins to {args.output}")

    # Plot if requested
    if args.plot:
        plot_directory_ = f"/user/mshoosht/public_html/Interpretations/Plots/EFTyields/{coeff_str}" 
        if not os.path.exists(plot_directory_):
           try:
                os.makedirs(plot_directory_)
           except OSError: # Resolve rare race condition
                pass
        
        shutil.copyfile('/user/mshoosht/public_html/index.php', os.path.join( plot_directory_, 'index.php' ) )

        for (year, ch_key, var, proc), yvals in plot_data.items():
            sm_vals = plot_sm[(year,ch_key, var, proc)]
            bins = range(len(sm_vals))
            print(f"[INFO] Plotting: Channel={ch_key}, Var={var}, Process={proc}")
            print(f"  [DEBUG] SM vals ({len(sm_vals)}): {sm_vals}")
            print(f"  [DEBUG] Yields ({len(yvals)}): {yvals}")

            plt.figure(figsize=(8, 5))
            plt.step(bins, sm_vals, where='mid', marker='x', linestyle='--', label='SM')
            plt.step(bins, yvals, where='mid', marker='o', label=proc)
            plt.xlabel('Bin')
            plt.ylabel('Yield')
            plt.title(f"Yields vs. Bin: {ch_key}/{var} ({proc})")
            plt.legend()
            plt.tight_layout()
            fname = plot_directory_ + f"/yields_{ch_key}_{var}_{proc}_{year}.png"
            plt.savefig(fname)
            print(f"[INFO] Saved plot: {fname}")


        # --- combined over all years AND all processes per channel/var ---
        chans = set((ck, v) for (_, ck, v, _) in plot_sm.keys())
        for ch_key, var in chans:
            # pick any representative key to get nbins
            sample = next(k for k in plot_sm if k[1]==ch_key and k[2]==var)
            nbins = len(plot_sm[sample])
            sm_totals  = [0.]*nbins
            eft_totals = [0.]*nbins
            # sum over every (year, proc)
            for (yr, ck, v, proc), sm_vals in plot_sm.items():
                if ck!=ch_key or v!=var: continue
                yvals = plot_data[(yr,ck,v,proc)]
                sm_totals  = [s + x for s, x in zip(sm_totals, sm_vals)]
                eft_totals = [e + y for e, y in zip(eft_totals, yvals)]
            plt.figure(figsize=(8,5))
            plt.step(range(nbins), sm_totals,  where='mid', marker='x', linestyle='--', label='SM total')
            plt.step(range(nbins), eft_totals, where='mid', marker='o', label='EFT total')
            plt.xlabel('Bin'); plt.ylabel('Yield')
            plt.title(f"Total Yields AllYears+AllProcs: {ch_key}/{var}")
            plt.legend(); plt.tight_layout()
            fname = f"{plot_directory_}/allYearsAllProcs_yields_{ch_key}_{var}.png"
            plt.savefig(fname)
            print(f"[INFO] Saved combined all‐years+processes plot: {fname}")


if __name__ == '__main__':
    main()

