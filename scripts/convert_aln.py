#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import AlignIO


# support alignment file formats
SUPPORTED_FORMATS = [
    "clustal",
    "phylip",
    "phylip-relaxed",
    "stockholm",
    "fasta",
]


class ConversionError(Exception):
    pass


def convert_aln(filename: Path, input_format: str, output_format: str) -> None:
    """Convert an alignment file into a specific format"""
    alignment = AlignIO.parse(filename, input_format)
    if output_format == "phylip-relaxed":
        output_filename = f"{filename.stem}.phylip"
    else:
        output_filename = f"{filename.stem}.{output_format}"
    output_path = filename.parent / output_filename
    write_success = AlignIO.write(alignment, output_path, output_format)

    if write_success:
        print(f"Converted {filename} to {output_format.upper()} format")
    else:
        raise ConversionError(f"Failed to convert {filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="convert_aln",
        description="convert between MSA formats",
    )

    parser.add_argument("filename", help="MSA file to convert", type=Path)
    parser.add_argument(
        "infmt", help="input file format", type=str, choices=SUPPORTED_FORMATS
    )
    parser.add_argument(
        "outfmt", help="target file format", type=str, choices=SUPPORTED_FORMATS
    )
    parser.add_argument(
        "-l", "--formats", help="list supported formats", action="store_true"
    )

    args = parser.parse_args()

    convert_aln(args.filename, args.infmt, args.outfmt)
