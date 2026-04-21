#!/usr/bin/env python3
"""
Inject KOfam score thresholds into HMM profiles as HMMER3 GA (Gathering
Threshold) lines, enabling use of --cut_ga with hmmsearch.

The KOfam ko_list supplies a threshold and score_type (full/domain) for each
KO.  Both GA fields are set to the same threshold value; HMMER --cut_ga
requires both the full-sequence score and the best-domain score to exceed
their respective GA values, so using the KOfam threshold for both is a
conservative approach that works well in practice.

Profiles for KOs with no defined threshold ('-' in ko_list) are passed
through unchanged.
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path


def parse_ko_list(ko_list_path: Path) -> dict[str, float]:
    """Return {knum: threshold} for entries with a numeric threshold."""
    thresholds = {}
    with ko_list_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            knum = row["knum"].strip()
            raw = row["threshold"].strip()
            if raw == "-":
                continue
            try:
                thresholds[knum] = float(raw)
            except ValueError:
                print(f"WARNING: cannot parse threshold '{raw}' for {knum}, skipping", file=sys.stderr)
    return thresholds


def inject_ga(profile_path: Path, threshold: float, output_path: Path) -> None:
    """Write profile to output_path with a GA line inserted after NAME."""
    ga_line = f"GA    {threshold} {threshold};\n"
    inserted = False
    lines = profile_path.read_text().splitlines(keepends=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as out:
        for line in lines:
            out.write(line)
            if not inserted and line.startswith("NAME "):
                out.write(ga_line)
                inserted = True
    if not inserted:
        print(f"WARNING: no NAME line found in {profile_path.name}; GA not injected", file=sys.stderr)


def process(
    ko_list_path: Path,
    profiles_dir: Path,
    output_dir: "Path | None",
    in_place: bool,
    missing_threshold: str,
) -> None:
    thresholds = parse_ko_list(ko_list_path)
    hmm_files = sorted(profiles_dir.glob("*.hmm"))

    if not hmm_files:
        print(f"WARNING: no .hmm files found in {profiles_dir}", file=sys.stderr)
        return

    injected = 0
    passed_through = 0

    for hmm_file in hmm_files:
        knum = hmm_file.stem
        dest = (hmm_file if in_place else output_dir / hmm_file.name)

        if knum in thresholds:
            inject_ga(hmm_file, thresholds[knum], dest)
            injected += 1
        else:
            # No threshold for this KO
            if missing_threshold == "error":
                print(f"ERROR: no threshold for {knum}", file=sys.stderr)
                sys.exit(1)
            elif missing_threshold == "warn":
                print(f"WARNING: no threshold for {knum}, copying unchanged", file=sys.stderr)
            if not in_place:
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(hmm_file, dest)
            passed_through += 1

    # Warn about KOs in ko_list that have no corresponding profile
    profile_knums = {f.stem for f in hmm_files}
    for knum in thresholds:
        if knum not in profile_knums:
            print(f"WARNING: {knum} has a threshold in ko_list but no matching .hmm file", file=sys.stderr)

    print(
        f"Done: {injected} profiles had GA injected, {passed_through} passed through unchanged.",
        file=sys.stderr,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--ko-list", required=True, type=Path, metavar="FILE",
                        help="KOfam ko_list file (TSV with knum/threshold/score_type columns)")
    parser.add_argument("--profiles-dir", required=True, type=Path, metavar="DIR",
                        help="Directory of individual .hmm profile files")

    dest_group = parser.add_mutually_exclusive_group(required=True)
    dest_group.add_argument("--output-dir", type=Path, metavar="DIR",
                            help="Write modified profiles to this directory")
    dest_group.add_argument("--in-place", action="store_true",
                            help="Modify profiles in-place (overwrites source files)")

    parser.add_argument(
        "--missing-threshold",
        choices=["skip", "warn", "error"],
        default="skip",
        help="Behaviour when a profile has no threshold in ko_list "
             "(skip=copy silently, warn=copy with warning, error=abort). Default: skip",
    )

    args = parser.parse_args()

    if not args.ko_list.is_file():
        parser.error(f"ko_list not found: {args.ko_list}")
    if not args.profiles_dir.is_dir():
        parser.error(f"profiles-dir not found: {args.profiles_dir}")
    if args.output_dir is not None and args.output_dir.resolve() == args.profiles_dir.resolve():
        parser.error("--output-dir must differ from --profiles-dir; use --in-place to modify in place")

    process(
        ko_list_path=args.ko_list,
        profiles_dir=args.profiles_dir,
        output_dir=args.output_dir,
        in_place=args.in_place,
        missing_threshold=args.missing_threshold,
    )


if __name__ == "__main__":
    main()
