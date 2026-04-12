"""Tests for pyclann."""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent.parent  # pyclann/ → clann/
EXAMPLES = REPO_ROOT / "examples"
CLANN_BIN = REPO_ROOT / "clann"

# Skip all tests when the clann binary is absent (e.g. pure Python CI).
pytestmark = pytest.mark.skipif(
    not CLANN_BIN.exists(),
    reason="clann binary not found at repo root",
)


import pyclann  # noqa: E402

pyclann.set_clann_path(str(CLANN_BIN))

SINGLE_PH = str(EXAMPLES / "tutorial_single.ph")
TREES_PH = str(EXAMPLES / "tutorial_trees.ph")
CANDIDATES_PH = str(EXAMPLES / "tutorial_candidates.ph")


# ---------------------------------------------------------------------------
# Parser unit tests (no binary needed)
# ---------------------------------------------------------------------------


class TestParsers:
    def test_parse_scored_tree_file(self):
        from pyclann._parser import parse_scored_tree_file

        text = (
            "((A,B),C);\t[7.297618]\n"
            "((A,C),B);\t[8.123456]\n"
            "\n"
        )
        results = parse_scored_tree_file(text)
        assert len(results) == 2
        assert results[0][0] == "((A,B),C)"
        assert abs(results[0][1] - 7.297618) < 1e-5
        assert results[1][0] == "((A,C),B)"

    def test_parse_plain_tree_file(self):
        from pyclann._parser import parse_plain_tree_file

        text = "((A,B),C);\n\n((X,Y),Z);\n"
        trees = parse_plain_tree_file(text)
        assert len(trees) == 2
        assert trees[0].startswith("(")

    def test_parse_bootstrap_file(self):
        from pyclann._parser import parse_bootstrap_file

        line = "((A,B),C) [1.000000];\t[score = 5.480159]\n"
        results = parse_bootstrap_file(line)
        assert len(results) == 1
        assert abs(results[0][1] - 5.480159) < 1e-5

    def test_extract_stdout_score(self):
        from pyclann._parser import extract_stdout_score

        stdout = "Supertree 1 of 1 score = 9.200000\n"
        score = extract_stdout_score(stdout)
        assert score is not None
        assert abs(score - 9.2) < 1e-5

    def test_extract_stdout_score_no_match(self):
        from pyclann._parser import extract_stdout_score

        assert extract_stdout_score("No score here\n") is None


# ---------------------------------------------------------------------------
# Core unit tests
# ---------------------------------------------------------------------------


class TestCore:
    def test_set_get_clann_path(self):
        old = pyclann.get_clann_path()
        pyclann.set_clann_path("/tmp/fake_clann")
        assert pyclann.get_clann_path() == "/tmp/fake_clann"
        pyclann.set_clann_path(old)

    def test_clann_error_raised_for_missing_binary(self):
        import pyclann._core as core

        saved = core._clann_path
        core._clann_path = "/nonexistent/clann"
        try:
            with pytest.raises(pyclann.ClannError):
                pyclann.hs(SINGLE_PH)
        finally:
            core._clann_path = saved

    def test_clann_error_raised_for_missing_treefile(self):
        with pytest.raises(pyclann.ClannError, match="not found"):
            pyclann.hs("/no/such/file.ph")


# ---------------------------------------------------------------------------
# Integration tests (require clann binary)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not CLANN_BIN.exists(), reason="clann binary not found")
class TestHeuristicSearch:
    def test_hs_returns_result(self):
        result = pyclann.hs(SINGLE_PH, nreps=1, nthreads=1)
        assert result.returncode == 0
        assert result.best_tree is not None
        assert result.best_tree.startswith("(")
        assert result.score is not None
        assert result.score > 0

    def test_hs_tree_list(self):
        result = pyclann.hs(SINGLE_PH, nreps=1, nthreads=1)
        assert len(result.trees) >= 1
        assert len(result.scores) == len(result.trees)

    def test_hs_command_recorded(self):
        result = pyclann.hs(SINGLE_PH, nreps=1, nthreads=1)
        assert "hs" in result.command
        assert "nreps=1" in result.command

    def test_hs_stdout_captured(self):
        result = pyclann.hs(SINGLE_PH, nreps=1, nthreads=1)
        assert len(result.stdout) > 0

    def test_hs_criterion_ml(self):
        result = pyclann.hs(SINGLE_PH, criterion="ml", nreps=1, nthreads=1)
        assert result.returncode == 0
        assert result.best_tree is not None


@pytest.mark.skipif(not CLANN_BIN.exists(), reason="clann binary not found")
class TestNJ:
    def test_nj_returns_result(self):
        result = pyclann.nj(SINGLE_PH)
        assert result.returncode == 0
        assert result.best_tree is not None
        assert result.best_tree.startswith("(")
        assert result.score is None  # NJ has no single score

    def test_nj_tree_has_branch_lengths(self):
        result = pyclann.nj(SINGLE_PH)
        # NJ tree includes branch lengths (:0.xxxxxx)
        assert ":" in result.best_tree


@pytest.mark.skipif(not CLANN_BIN.exists(), reason="clann binary not found")
class TestConsensus:
    def test_consensus_returns_result(self):
        result = pyclann.consensus(SINGLE_PH)
        assert result.returncode == 0
        assert result.best_tree is not None
        assert result.best_tree.startswith("(")

    def test_consensus_majority_rule(self):
        result = pyclann.consensus(SINGLE_PH, percentage=0.5)
        assert result.returncode == 0
        assert result.best_tree is not None


@pytest.mark.skipif(not CLANN_BIN.exists(), reason="clann binary not found")
class TestAllTrees:
    def test_alltrees_small_dataset(self):
        # tutorial_trees.ph has 6 taxa — manageable for alltrees
        result = pyclann.alltrees(TREES_PH)
        assert result.returncode == 0
        assert result.best_tree is not None
        assert result.score is not None


@pytest.mark.skipif(not CLANN_BIN.exists(), reason="clann binary not found")
class TestUserTrees:
    def test_usertrees_returns_scores(self):
        result = pyclann.usertrees(SINGLE_PH, CANDIDATES_PH, criterion="dfit")
        assert result.returncode == 0
        assert result.best_tree is not None
        assert len(result.scores) > 0

    def test_usertrees_missing_candidates(self):
        with pytest.raises(pyclann.ClannError, match="not found"):
            pyclann.usertrees(SINGLE_PH, "/no/such/candidates.ph")


@pytest.mark.skipif(not CLANN_BIN.exists(), reason="clann binary not found")
class TestBootstrap:
    def test_bootstrap_returns_consensus(self):
        result = pyclann.bootstrap(SINGLE_PH, nreps=5, nthreads=1)
        assert result.returncode == 0
        assert result.best_tree is not None
        assert result.best_tree.startswith("(")
        assert len(result.trees) == 5

    def test_bootstrap_scores_parallel_to_trees(self):
        result = pyclann.bootstrap(SINGLE_PH, nreps=3, nthreads=1)
        assert len(result.scores) == len(result.trees)


@pytest.mark.skipif(not CLANN_BIN.exists(), reason="clann binary not found")
class TestRun:
    def test_run_nj(self):
        result = pyclann.run("nj", SINGLE_PH)
        assert result.returncode == 0
        assert len(result.stdout) > 0
        assert result.best_tree is None  # run() doesn't parse output files

    def test_run_consensus(self):
        result = pyclann.run("consensus", SINGLE_PH)
        assert result.returncode == 0
        assert len(result.stdout) > 0
