import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
WATCHER = ROOT / "watch_iqtree_progress_htop.sh"


def generator():
    loader = importlib.machinery.SourceFileLoader("stage12_progress_test", str(ROOT / "stage12_prepare_supermatrix_and_gene_tree_datasets"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@pytest.fixture
def runtime(tmp_path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    iqtree = bindir / "iqtree"
    iqtree.write_text("#!" + sys.executable + "\n" + '''import os, sys
from pathlib import Path
args = sys.argv[1:]
assert Path(args[args.index('-s') + 1]).is_file(), args
prefix = Path(args[args.index('--prefix') + 1])
prefix.with_suffix('.log').write_text('Total wall-clock time used: 1 sec\\n')
if os.environ.get('TEST_IQTREE_FAIL'):
    sys.exit(7)
prefix.with_suffix('.treefile').write_text('(A,B);\\n')
''')
    iqtree.chmod(0o755)
    tmux = bindir / "tmux"
    tmux.write_text("#!" + sys.executable + "\n" + '''import json, os, sys
with open(os.environ['TMUX_CAPTURE'], 'a') as handle:
    handle.write(json.dumps(sys.argv[1:]) + '\\n')
''')
    tmux.chmod(0o755)
    return {**os.environ, "PATH": str(bindir) + os.pathsep + os.environ["PATH"],
            "TMUX": "test-session", "TMUX_CAPTURE": str(tmp_path / "tmux.jsonl"), "TERM": "xterm"}


@pytest.mark.parametrize("failed", [False, True])
def test_supermatrix_launches_one_monitor_and_preserves_exit_status(tmp_path, runtime, failed):
    folder = tmp_path / "analysis with spaces"
    folder.mkdir()
    phy = folder / "matrix.phylip"
    phy.write_text("2 4\nA ACGT\nB AGGT\n")
    charsets = folder / "matrix.nexus"
    charsets.write_text("#nexus\n")
    script = folder / "run.sh"
    generator().write_iqtree_script(script, phy, charsets, "Project name")
    if failed:
        runtime["TEST_IQTREE_FAIL"] = "1"
    result = subprocess.run(["bash", str(script)], cwd=tmp_path, env=runtime, capture_output=True, text=True)
    assert result.returncode == (7 if failed else 0), result.stderr
    state = (folder / "Project name_matrix.progress_status").read_text().strip()
    assert state == ("FAILED (exit 7)" if failed else "FINISHED")
    calls = [json.loads(line) for line in Path(runtime["TMUX_CAPTURE"]).read_text().splitlines()]
    assert len(calls) == 1
    command = shlex.split(calls[0][-1])
    assert Path(command[1]) == folder / WATCHER.name
    assert Path(command[2]) == folder / "Project name_matrix.log"
    assert Path(command[-1]) == folder / "Project name_matrix.progress_status"
    # A failure status takes precedence over a finish-looking log.
    monitor = subprocess.run(command, env=runtime, text=True, capture_output=True, timeout=5)
    assert monitor.returncode == 0
    assert state in monitor.stdout


def test_gene_tree_launches_one_batch_monitor_and_runs_each_locus(tmp_path, runtime):
    folder = tmp_path / "gene trees"
    (folder / "alignments").mkdir(parents=True)
    for locus in ("uce-1", "uce 2"):
        (folder / "alignments" / (locus + ".fasta")).write_text(">A\nACGT\n>B\nAGGT\n")
    generator().write_gene_tree_analysis_scripts(folder, "Project")
    runtime["JOBS"] = "2"
    result = subprocess.run(["bash", str(folder / "01_run_iqtree_gene_trees.sh")], cwd=tmp_path,
                            env=runtime, text=True, capture_output=True, timeout=10)
    assert result.returncode == 0, result.stderr
    assert len(list((folder / "01_iqtree_gene_trees").glob("*.treefile"))) == 2
    calls = Path(runtime["TMUX_CAPTURE"]).read_text().splitlines()
    assert len(calls) == 1
    command = shlex.split(json.loads(calls[0])[-1])
    assert command[4] == "2"
    monitor = subprocess.run(command, env=runtime, text=True, capture_output=True, timeout=5)
    assert "2 / 2 loci finished" in monitor.stdout
    assert "FINISHED" in monitor.stdout


def test_monitor_waits_for_new_log_and_handles_missing_elapsed(tmp_path):
    log = tmp_path / "later.log"
    with subprocess.Popen(["bash", str(WATCHER), str(log), "0.05"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        assert "WAITING" in process.stdout.readline()
        log.write_text("Generating 1000 samples for ultrafast bootstrap\nIteration 2 / LogL: -100\nTotal wall-clock time used: 1 sec\n")
        stdout, stderr = process.communicate(timeout=5)
        assert process.returncode == 0, stderr
        assert "FINISHED" in stdout


def test_disable_automatic_monitor(tmp_path, runtime):
    phy = tmp_path / "matrix.phylip"
    phy.write_text("2 4\nA ACGT\nB AGGT\n")
    script = tmp_path / "run.sh"
    generator().write_iqtree_script(script, phy, tmp_path / "matrix.nexus", "Project")
    runtime["IQTREE_PROGRESS"] = "0"
    subprocess.run(["bash", str(script)], env=runtime, check=True)
    assert not Path(runtime["TMUX_CAPTURE"]).exists()
