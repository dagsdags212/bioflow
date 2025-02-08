#!/usr/bin/env python

import argparse
import sys
from io import StringIO
from pathlib import Path

import pandas as pd
import seaborn.objects as so
from seaborn import axes_style


def load_coverage_data(inp):
    """Load input from stdin or file as a dataframe."""
    rstr = StringIO(inp)
    names = ["accession", "position", "coverage"]
    return pd.read_csv(rstr, sep="\t", header=None, names=names)

def plot_coverage(df):
    """Generate a coverage plot using seaborn."""
    # Subset dataframe.
    cov_df = df[["position", "coverage"]]

    # Arrays for settings axes limits.
    positions = cov_df["position"].to_numpy()
    coverage = cov_df["coverage"].to_numpy()

    # Reference file accession.
    acc = df["accession"][0]

    # Set theme config.
    theme_dict = {
        **axes_style("white"), 
        "font.size": 12.0,
        "axes.labelpad": 8.0,
        "axes.labelweight": 500,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.spines.left": False,
        "axes.spines.bottom": True,
        "xtick.labelsize": 10,
        "xtick.bottom": True,
        "xtick.major.size": 3.5,
        "ytick.labelsize": 10,
    }

    # Generate plot.
    buffer = 150
    covplot = \
    so.Plot(cov_df, x="position", y="coverage") \
        .add(so.Area()) \
        .limit(
            x=(positions.min() - buffer, positions.max() + buffer),
            y=(coverage.min() - buffer, coverage.max() + buffer * 3)
        ) \
        .label(
            x="Position", y="Coverage",
            title=f"Coverage of reads mapped against {acc}"
        ) \
        .theme(theme_dict)
    
    return covplot

def main():
    parser = argparse.ArgumentParser(
        prog="plot_coverage.py",
        description="Generate a plot for visualizing read coverage",
        usage="%(prog)s <covfile> [options]",
        epilog="Written by Jan Samson"
    )
    
    parser.add_argument(
        "covfile",
        help="data for read coverage at each genome position",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    parser.add_argument(
        "-o", "--output",
        help="filepath for storing the plot",
        default=Path("./covplot.png")
    )
    
    args = parser.parse_args()
    
    if args.covfile.isatty():
        parser.print_help()
        return -1

    df = load_coverage_data(args.covfile.read())
    covplot = plot_coverage(df)
    covplot.save(args.output, dpi=240)

if __name__ == "__main__":
    sys.exit(main())
