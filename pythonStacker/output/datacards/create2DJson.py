#!/usr/bin/env python3
import numpy as np
import uproot
import json
import matplotlib.pyplot as plt
from scipy.stats import chi2

def sigma_to_delta_nll(sigma, dof=2):
    """
    Convert sigma confidence level to delta NLL for given degrees of freedom.
    """
    return chi2.ppf(chi2.cdf(sigma**2, 1), dof)


def extract_contours_from_ttree(root_filename, x_branch='kappa', y_branch='tildekappa', delta_nll_branch='deltaNLL'):
    """
    Extract 1σ and 2σ 2D contours from a ROOT file TTree named 'limit'.

    Returns:
        dict: { '1sigma': [...], '2sigma': [...] }
    """
    # Load branches
    with uproot.open(root_filename) as f:
        tree = f['limit']
        x = tree[x_branch].array(library='np')
        y = tree[y_branch].array(library='np')
        delta_nll = tree[delta_nll_branch].array(library='np')

    # Define contour levels
    level_map = {'1sigma': 2.3, '2sigma': 6.2}

    # Build regular grid
    xi = np.unique(x)
    yi = np.unique(y)
    X, Y = np.meshgrid(xi, yi)
    Z = np.full_like(X, np.nan, dtype=float)

    for xv, yv, zn in zip(x, y, delta_nll):
        ix = np.where(xi == xv)[0][0]
        iy = np.where(yi == yv)[0][0]
        Z[iy, ix] = zn

    # Extract contours
    fig, ax = plt.subplots()
    contours_out = {}
    for label, lvl in level_map.items():
        cs = ax.contour(X, Y, Z, levels=[lvl])
        paths = cs.collections[0].get_paths()
        contours = []
        for path in paths:
            verts = path.vertices
            contours.append({'x': verts[:, 0].tolist(), 'y': verts[:, 1].tolist()})
        contours_out[label] = contours
    plt.close(fig)
    return contours_out


def extract_bestfit_from_ttree(root_filename, x_branch='kappa', y_branch='tildekappa', delta_nll_branch='deltaNLL'):
    """
    Find the best-fit point (minimum deltaNLL) in the TTree.

    Returns:
        dict: {'x': float, 'y': float}
    """
    with uproot.open(root_filename) as f:
        tree = f['limit']
        x = tree[x_branch].array(library='np')
        y = tree[y_branch].array(library='np')
        delta_nll = tree[delta_nll_branch].array(library='np')

    # Identify index of minimum deltaNLL (should be zero for observed)
    idx = np.nanargmin(delta_nll)
    return {'x': float(x[idx]), 'y': float(y[idx])}


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Extract 2D contours and best-fit point from ROOT TTree.')
    parser.add_argument('input_root', help='Input ROOT file containing TTree "limit"')
    parser.add_argument('-x', '--xbranch', default='kappa', help='Name of the x variable branch')
    parser.add_argument('-y', '--ybranch', default='tildekappa', help='Name of the y variable branch')
    parser.add_argument('-d', '--delta', '--delta_nll_branch', dest='delta_branch',
                        default='deltaNLL', help='Name of the delta NLL branch')
    parser.add_argument('-o', '--output_json', required=True, help='Name of the output JSON file')
    args = parser.parse_args()

    # Extract contours and best-fit
    contours = extract_contours_from_ttree(
        args.input_root,
        x_branch=args.xbranch,
        y_branch=args.ybranch,
        delta_nll_branch=args.delta_branch
    )
    best_fit = extract_bestfit_from_ttree(
        args.input_root,
        x_branch=args.xbranch,
        y_branch=args.ybranch,
        delta_nll_branch=args.delta_branch
    )

    # Combine and save
    output = contours.copy()
    output['bestFit'] = best_fit
    with open(args.output_json, 'w') as f:
        json.dump(output, f, indent=2)

