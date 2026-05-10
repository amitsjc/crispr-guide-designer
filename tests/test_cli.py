"""Smoke tests for the CLI."""

import subprocess, sys, os

PYTHON = sys.executable
VEGFA  = os.path.join(os.path.dirname(__file__), "..", "data", "vegfa.fasta")
BG     = os.path.join(os.path.dirname(__file__), "..", "data", "genome_sample.fasta")
ROOT   = os.path.join(os.path.dirname(__file__), "..")


def run(*args):
    return subprocess.run(
        [PYTHON, "-m", "crispr_designer", *args],
        capture_output=True, text=True, cwd=ROOT,
    )


def test_design_exits_zero():
    assert run("design", VEGFA).returncode == 0

def test_design_shows_ranked_guides():
    r = run("design", VEGFA)
    assert "Ranked Guide RNAs" in r.stdout

def test_design_with_offtarget():
    r = run("design", VEGFA, "--offtarget", BG)
    assert r.returncode == 0
    assert "Off-Target Analysis" in r.stdout

def test_pamsites_command():
    r = run("pamsites", VEGFA)
    assert r.returncode == 0
    assert "PAM Sites" in r.stdout

def test_strict_flag_accepted():
    r = run("design", VEGFA, "--strict")
    assert r.returncode == 0

def test_missing_file_nonzero():
    r = run("design", "nonexistent.fasta")
    assert r.returncode != 0

def test_help():
    r = run("--help")
    assert r.returncode == 0
    assert "design" in r.stdout
    assert "pamsites" in r.stdout
