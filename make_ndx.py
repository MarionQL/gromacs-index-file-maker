#!/usr/bin/env python3
"""
make_ndx.py

Generate GROMACS index (.ndx) files using MDAnalysis selection syntax.

This script is designed to be:
- robust (validates inputs and selections)
- generalizable (works for any system: protein, membrane, nucleic acids)
- CLI-friendly (fits into HPC / SLURM workflows)
- reproducible (explicit indexing control and deterministic output)
"""

import argparse
import sys
import MDAnalysis as mda

def parse_args():
    """
    Define and parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed CLI arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Generate GROMACS index (.ndx) files from MDAnalysis selections.\n"
            "Supports arbitrary atom selections and produces properly formatted output."
        )
    )

    # --- Required inputs ---
    parser.add_argument(
        "-s", "--structure",
        required=True,
        help="Structure file (e.g., .gro, .pdb, .tpr)"
    )

    parser.add_argument(
        "--selection",
        required=True,
        help=(
            "MDAnalysis selection string.\n"
            "Example: 'segid MEMB and not resname CHL1'"
        )
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output .ndx file"
    )

    # --- Optional inputs ---
    parser.add_argument(
        "-f", "--trajectory",
        default=None,
        help="Optional trajectory file (.xtc, .trr). Not required for static selections."
    )

    parser.add_argument(
        "--name",
        default="GROUP",
        help="Name of the index group (default: GROUP)"
    )

    parser.add_argument(
        "--one-based",
        action="store_true",
        help=(
            "Convert indices to 1-based indexing (required for GROMACS tools).\n"
            "Strongly recommended."
        )
    )

    return parser.parse_args()

def write_ndx(filename, groups):
    """
    Write GROMACS index file.

    Parameters
    ----------
    filename : str
        Output file path.
    groups : dict[str, list[int]]
        Dictionary mapping group names → atom indices.

    Notes
    -----
    - GROMACS expects 1-based indexing.
    - Standard formatting is ~15 indices per line.
    - This function does NOT enforce indexing convention;
      caller must ensure correctness.
    """
    try:
        with open(filename, "w") as f:
            for name, indices in groups.items():

                # Write group header
                f.write(f"[ {name} ]\n")

                # Write indices with fixed-width formatting
                for i, idx in enumerate(indices):
                    f.write(f"{idx:8d}")

                    # Insert newline every 15 entries for readability
                    if (i + 1) % 15 == 0:
                        f.write("\n")

                f.write("\n\n")

    except IOError as e:
        raise RuntimeError(f"Failed to write index file: {e}")

def main():
    """
    Main execution pipeline:
    1. Load system (structure + optional trajectory)
    2. Apply selection
    3. Validate selection
    4. Convert indexing if necessary
    5. Write .ndx file
    """

    args = parse_args()

    try:
        if args.trajectory:
            # Full trajectory context (not required for static selections)
            u = mda.Universe(args.structure, args.trajectory)
        else:
            # Structure-only mode (faster, sufficient for index generation)
            u = mda.Universe(args.structure)

    except Exception as e:
        print(f"[ERROR] Failed to load system:\n{e}", file=sys.stderr)
        sys.exit(1)

    try:
        ag = u.select_atoms(args.selection)

    except Exception as e:
        print(f"[ERROR] Invalid selection string:\n{args.selection}", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)

    if ag.n_atoms == 0:
        print("[ERROR] Selection returned 0 atoms.", file=sys.stderr)
        print(f"Selection: {args.selection}", file=sys.stderr)
        sys.exit(1)

    if args.one_based:
        # Convert to GROMACS-compatible indexing
        indices = [atom.index + 1 for atom in ag]
    else:
        # Raw MDAnalysis indexing (0-based)
        indices = [atom.index for atom in ag]

    groups = {args.name: indices}
    try:
        write_ndx(args.output, groups)

    except Exception as e:
        print(f"[ERROR] Failed to write output file:\n{e}", file=sys.stderr)
        sys.exit(1)

    print(f"[OK] Wrote: {args.output}")
    print(f"[INFO] Group name: {args.name}")
    print(f"[INFO] Atom count: {len(indices)}")
    print(f"[INFO] Indexing: {'1-based (GROMACS)' if args.one_based else '0-based'}")

if __name__ == "__main__":
    main()
