"""
Tests for inject_kofam_ga_thresholds.py

All fixtures are synthetic minimal HMMER3 headers — no dependency on real
test databases.
"""

import sys
from pathlib import Path

import pytest

# Allow importing the script directly from the parent directory
sys.path.insert(0, str(Path(__file__).parent.parent))
from inject_kofam_ga_thresholds import parse_ko_list, process

# ---------------------------------------------------------------------------
# Shared fixture helpers
# ---------------------------------------------------------------------------

MINIMAL_HMM = """\
HMMER3/f [3.2.1 | June 2018]
NAME  {knum}
LENG  10
ALPH  amino
HMM          A        C
//
"""

KO_LIST_HEADER = "knum\tthreshold\tscore_type\tprofile_type\tdefinition\n"


def make_profiles(directory: Path, knums: list[str]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for knum in knums:
        (directory / f"{knum}.hmm").write_text(MINIMAL_HMM.format(knum=knum))


def make_ko_list(path: Path, rows: list[tuple]) -> None:
    """rows: list of (knum, threshold, score_type)"""
    lines = [KO_LIST_HEADER]
    for knum, threshold, score_type in rows:
        lines.append(f"{knum}\t{threshold}\t{score_type}\tall\t{knum} definition\n")
    path.write_text("".join(lines))


# ---------------------------------------------------------------------------
# parse_ko_list
# ---------------------------------------------------------------------------

def test_parse_ko_list_numeric(tmp_path):
    ko = tmp_path / "ko_list"
    make_ko_list(ko, [("K00001", "345.37", "domain"), ("K00002", "453.33", "full")])
    result = parse_ko_list(ko)
    assert result == {"K00001": 345.37, "K00002": 453.33}


def test_parse_ko_list_skips_dash(tmp_path):
    ko = tmp_path / "ko_list"
    make_ko_list(ko, [("K00001", "345.37", "domain"), ("K00003", "-", "domain")])
    result = parse_ko_list(ko)
    assert "K00003" not in result
    assert "K00001" in result


# ---------------------------------------------------------------------------
# GA injection — output-dir mode
# ---------------------------------------------------------------------------

def test_ga_injected_domain(tmp_path):
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001"])
    make_ko_list(tmp_path / "ko_list", [("K00001", "345.37", "domain")])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="skip")

    lines = (out / "K00001.hmm").read_text().splitlines()
    name_idx = next(i for i, l in enumerate(lines) if l.startswith("NAME "))
    assert lines[name_idx + 1] == "GA    345.37 345.37;"


def test_ga_injected_full(tmp_path):
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00002"])
    make_ko_list(tmp_path / "ko_list", [("K00002", "453.33", "full")])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="skip")

    lines = (out / "K00002.hmm").read_text().splitlines()
    name_idx = next(i for i, l in enumerate(lines) if l.startswith("NAME "))
    assert lines[name_idx + 1] == "GA    453.33 453.33;"


def test_ga_line_immediately_after_name(tmp_path):
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001"])
    make_ko_list(tmp_path / "ko_list", [("K00001", "100.0", "full")])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="skip")

    lines = (out / "K00001.hmm").read_text().splitlines()
    name_idx = next(i for i, l in enumerate(lines) if l.startswith("NAME "))
    assert lines[name_idx + 1].startswith("GA ")


def test_no_ga_for_missing_threshold(tmp_path):
    """Profiles with threshold '-' must not receive a GA line."""
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00003"])
    make_ko_list(tmp_path / "ko_list", [("K00003", "-", "domain")])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="skip")

    content = (out / "K00003.hmm").read_text()
    assert "GA " not in content
    # File should otherwise be identical to source
    assert content == (profiles / "K00003.hmm").read_text()


def test_source_profiles_untouched_in_output_mode(tmp_path):
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001"])
    original = (profiles / "K00001.hmm").read_text()
    make_ko_list(tmp_path / "ko_list", [("K00001", "200.0", "full")])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="skip")

    assert (profiles / "K00001.hmm").read_text() == original


# ---------------------------------------------------------------------------
# In-place mode
# ---------------------------------------------------------------------------

def test_in_place_mode(tmp_path):
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001"])
    make_ko_list(tmp_path / "ko_list", [("K00001", "77.5", "full")])

    process(tmp_path / "ko_list", profiles, output_dir=None, in_place=True, missing_threshold="skip")

    lines = (profiles / "K00001.hmm").read_text().splitlines()
    assert any(l == "GA    77.5 77.5;" for l in lines)


# ---------------------------------------------------------------------------
# Missing profile / warning behaviour
# ---------------------------------------------------------------------------

def test_missing_profile_warns(tmp_path, capsys):
    """KO in ko_list but no .hmm file → warning, no crash."""
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001"])
    make_ko_list(tmp_path / "ko_list", [
        ("K00001", "100.0", "full"),
        ("K99999", "200.0", "full"),  # no matching file
    ])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="skip")

    stderr = capsys.readouterr().err
    assert "K99999" in stderr


def test_missing_threshold_warn_mode(tmp_path, capsys):
    """Profile with no ko_list entry → warning printed when missing_threshold='warn'."""
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001", "K00002"])
    make_ko_list(tmp_path / "ko_list", [("K00001", "100.0", "full")])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="warn")

    stderr = capsys.readouterr().err
    assert "K00002" in stderr


def test_missing_threshold_error_mode(tmp_path):
    """missing_threshold='error' aborts when a profile has no threshold."""
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001", "K00002"])
    make_ko_list(tmp_path / "ko_list", [("K00001", "100.0", "full")])
    out = tmp_path / "out"

    with pytest.raises(SystemExit):
        process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="error")


# ---------------------------------------------------------------------------
# Three-profile integration smoke test
# ---------------------------------------------------------------------------

def test_three_profiles(tmp_path):
    """End-to-end with 3 profiles: domain threshold, full threshold, no threshold."""
    profiles = tmp_path / "profiles"
    make_profiles(profiles, ["K00001", "K00002", "K00003"])
    make_ko_list(tmp_path / "ko_list", [
        ("K00001", "345.37", "domain"),
        ("K00002", "453.33", "full"),
        ("K00003", "-",      "domain"),
    ])
    out = tmp_path / "out"

    process(tmp_path / "ko_list", profiles, out, in_place=False, missing_threshold="skip")

    def ga_line(knum):
        lines = (out / f"{knum}.hmm").read_text().splitlines()
        return next((l for l in lines if l.startswith("GA ")), None)

    assert ga_line("K00001") == "GA    345.37 345.37;"
    assert ga_line("K00002") == "GA    453.33 453.33;"
    assert ga_line("K00003") is None
