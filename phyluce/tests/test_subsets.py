"""Exercise isolated subset creation and real stage 7-12 data processing.

Only MAFFT is replaced by an identity executable on prealigned toy sequences;
these tests do not claim to validate MAFFT or infer biological phylogenies.
"""

import csv
import importlib.machinery
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from phyluce import subsets


ROOT = Path(__file__).resolve().parents[2]


def launcher_module():
    loader = importlib.machinery.SourceFileLoader("subset_launcher_test", str(ROOT / "phyluce_launcher"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@pytest.fixture
def project(tmp_path):
    project = tmp_path / "CompleteProject"
    project.mkdir()
    fake_mafft = tmp_path / "mafft"
    fake_mafft.write_text("#!{}\nimport sys\nfrom pathlib import Path\nprint(Path(sys.argv[-1]).read_text(), end='')\n".format(sys.executable))
    fake_mafft.chmod(0o755)
    (project / "project.info").write_text(
        "project_id\tFull\nthreads\t1\nmafft_threads_per_job\t1\nmax_mafft_jobs\t1\n"
        "mafft_path\t{}\ntrim_mode\tnone\ndataset_modes\tcomplete,incomplete\n"
        "align_datasets\tall\npost_align_thresholds\t1.0,0.75\n".format(fake_mafft)
    )
    fields = ["UCE locus", "Target name", "Target start", "Target end", "Strand", "Probe hit count", "Query names"]
    for sample in "ABCDEF":
        # uce-2 becomes complete when E and F are removed. uce-3 remains
        # incomplete (3/4) in that subset, and is absent in the other one.
        loci = ["uce-1"] + (["uce-2"] if sample in "ABCD" else []) + (["uce-3"] if sample in "ABC" else [])
        hits = project / "05_Filtering_uce_hits_trace" / sample / (sample + ".filtered.loci.tsv")
        hits.parent.mkdir(parents=True)
        with open(hits, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fields)
            for locus in loci:
                writer.writerow([locus, "contig1", 0, 120, "+", 1, locus + "_p1"])
        fasta = project / "06_extract_uce_fasta_from_hits" / sample / (sample + ".filtered.fasta")
        fasta.parent.mkdir(parents=True)
        fasta.write_text("".join(">{}|{}|contig1\n{}\n".format(sample, locus, "ACGT" * 30) for locus in loci))
    for stage in range(1, 13):
        directory = project / launcher_module().STAGE_NAMES[stage]
        directory.mkdir(exist_ok=True)
        (directory / "original_result.txt").write_text("Do not modify stage {}\n".format(stage))
    return project


def cli(project, *options, check=True):
    return subprocess.run(
        [sys.executable, str(ROOT / "phyluce_launcher"), "--project-id", str(project), "--start-stage", "13", *options],
        text=True, capture_output=True, check=check, cwd=project.parent,
    )


def original_hashes(project):
    return {str(path.relative_to(project)): subsets.digest(path)
            for path in project.rglob("*") if path.is_file() and subsets.ROOT_NAME not in path.parts}


def metadata(branch):
    return json.loads((branch / subsets.RECORDS / subsets.INFO).read_text())


def branch_samples(branch, stage, filename):
    with open(branch / launcher_module().STAGE_NAMES[stage] / filename) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_two_independent_subsets_rebuild_loci_and_preserve_original(project):
    before = original_hashes(project)
    cli(project, "--subset-action", "create", "--subset-name", "第一組")
    first = next((project / subsets.ROOT_NAME).glob("第一組_*"))
    assert subsets.read_sample_list(first / subsets.EDITABLE) == list("ABCDEF")
    assert metadata(first)["created_at"].endswith("+08:00")
    subsets.write_samples(first / subsets.EDITABLE, list("ABCD"))
    result = cli(project, "--subset-action", "run", "--subset-dir", str(first))
    assert "Stage 7-12 preparation complete" in result.stdout
    assert metadata(first)["status"] == "prepared"
    assert subsets.read_sample_list(first / subsets.RECORDS / subsets.EXCLUDED) == list("EF")
    rows = branch_samples(first, 7, "Full.complete.manifest.tsv")
    assert {row["uce_locus"] for row in rows} == {"uce-1", "uce-2"}
    filters = branch_samples(first, 11, "all_thresholds_filter_manifest.tsv")
    locus = [row for row in filters if row["dataset_label"] == "incomplete" and row["uce_locus"] == "uce-3"]
    assert {(row["filter_label"], row["status"]) for row in locus} == {("complete_matrix", "dropped"), ("min_taxa_075", "kept")}
    assert all(row["total_taxa"] == "4" for row in filters)

    cli(project, "--subset-action", "create", "--subset-name", "second")
    second = next((project / subsets.ROOT_NAME).glob("second_*"))
    assert subsets.read_sample_list(second / subsets.EDITABLE) == list("ABCDEF")
    subsets.write_samples(second / subsets.EDITABLE, list("CDEF"))
    cli(project, "--subset-action", "run", "--subset-dir", second.name)
    rows = branch_samples(second, 7, "Full.complete.manifest.tsv")
    assert {row["sample_id"] for row in rows} == set("CDEF")
    assert {row["uce_locus"] for row in rows} == {"uce-1"}
    assert original_hashes(project) == before
    for branch, selected in [(first, set("ABCD")), (second, set("CDEF"))]:
        assert not (branch / "06_extract_uce_fasta_from_hits").exists()
        outputs = branch / "12_Analysis_Branches"
        assert (outputs / "supermatrix").is_dir() and (outputs / "gene-tree").is_dir()
        for fasta in outputs.rglob("*.fasta"):
            headers = {line[1:].split("|")[0] for line in fasta.read_text().splitlines() if line.startswith(">")}
            assert headers <= selected
        manifest_text = (outputs / "branch_manifest.tsv").read_text()
        assert str(branch) in manifest_text
        assert str(second if branch == first else first) not in manifest_text
    # Running a prepared branch again must not rewrite its outputs.
    timestamp = (first / "12_Analysis_Branches" / "branch_manifest.tsv").stat().st_mtime_ns
    assert "already prepared" in cli(project, "--subset-action", "run", "--subset-dir", first.name).stdout
    assert (first / "12_Analysis_Branches" / "branch_manifest.tsv").stat().st_mtime_ns == timestamp

    stage11_manifest = first / "11_Filtered_Locus_ALIGNMENTS" / "threshold_manifest.tsv"
    stage11_timestamp = stage11_manifest.stat().st_mtime_ns
    shutil.rmtree(first / "12_Analysis_Branches")
    resumed = cli(project, "--subset-action", "run", "--subset-dir", first.name)
    assert "Stage 12 output is missing or empty" in resumed.stdout
    assert stage11_manifest.stat().st_mtime_ns == stage11_timestamp
    assert (first / "12_Analysis_Branches" / "branch_manifest.tsv").is_file()
    assert metadata(first)["status"] == "prepared"

    subsets.write_samples(first / subsets.EDITABLE, list("ABCDE"))
    error = cli(project, "--subset-action", "run", "--subset-dir", first.name, check=False)
    assert error.returncode != 0 and "Sample list changed" in error.stderr


@pytest.mark.parametrize("content, message", [
    ("A\nB\nA\n", "Duplicate"), ("A\nB\nunknown\n", "Unknown"),
    ("\n", "empty"), ("A\nB\n", "at least 3"), ("A\nB\nC D\n", "Invalid"),
])
def test_invalid_lists_do_not_start_analysis(project, content, message):
    branch = subsets.create_subset(project, "invalid", launcher_module())
    (branch / subsets.EDITABLE).write_text(content)
    result = cli(project, "--subset-action", "run", "--subset-dir", str(branch), check=False)
    assert result.returncode != 0 and message in result.stderr
    assert metadata(branch)["status"] == "draft"
    assert not (branch / "07_Aggregated_loci").exists()
    assert not (branch / subsets.RECORDS / subsets.USED).exists()


def test_copy_name_collision_dry_run_and_source_change(project, monkeypatch):
    launcher = launcher_module()
    instant = subsets.now()
    monkeypatch.setattr(subsets, "now", lambda: instant)
    first = subsets.create_subset(project, "same name", launcher)
    subsets.write_samples(first / subsets.EDITABLE, list("ABCD"))
    second = subsets.create_subset(project, "same name", launcher, copy_from=first)
    assert second.name == first.name + "_02"
    assert subsets.read_sample_list(second / subsets.EDITABLE) == list("ABCD")
    assert set(metadata(second)["sources"]) == set("ABCDEF")
    before = {str(p): subsets.digest(p) for p in project.rglob("*") if p.is_file()}
    cli(project, "--subset-action", "run", "--subset-dir", str(second), "--dry-run")
    after = {str(p): subsets.digest(p) for p in project.rglob("*") if p.is_file()}
    assert before == after
    source = Path(metadata(second)["sources"]["A"]["fasta"])
    source.write_text(source.read_text().replace("ACGT", "TGCA"))
    error = cli(project, "--subset-action", "run", "--subset-dir", str(second), check=False)
    assert error.returncode != 0 and "source changed" in error.stderr
    assert not (second / subsets.RECORDS / subsets.USED).exists()


def test_failed_run_resumes_from_last_completed_stage(project, monkeypatch):
    launcher = launcher_module()
    branch = subsets.create_subset(project, "resume", launcher)
    subsets.write_samples(branch / subsets.EDITABLE, list("ABCD"))
    monkeypatch.setattr(sys, "argv", ["phyluce_launcher"])
    args = launcher.get_args()
    original = subsets.stream_stage
    called = []

    def fail_at_stage9(command, *args, **kwargs):
        called.append(Path(command[1]).name)
        if "stage9_" in str(command[1]):
            raise RuntimeError("simulated interruption")
        return original(command, *args, **kwargs)

    monkeypatch.setattr(subsets, "stream_stage", fail_at_stage9)
    with pytest.raises(RuntimeError, match="simulated interruption"):
        subsets.run_subset(project, branch, args, launcher)
    assert metadata(branch)["last_completed_stage"] == 8
    assert metadata(branch)["status"] == "failed"
    stage7_before = (branch / "07_Aggregated_loci" / "Full.uce.sqlite").stat().st_mtime_ns
    monkeypatch.setattr(subsets, "stream_stage", original)
    subsets.run_subset(project, branch, args, launcher)
    assert metadata(branch)["status"] == "prepared"
    assert (branch / "07_Aggregated_loci" / "Full.uce.sqlite").stat().st_mtime_ns == stage7_before


def test_stage13_is_separate_and_cannot_escape_output_root(project):
    launcher = launcher_module()
    assert launcher.parse_interactive_stage_choice("f") == (1, 12)
    assert launcher.parse_interactive_stage_choice("13") == (13, 13)
    with pytest.raises(ValueError):
        launcher.parse_interactive_stage_choice("7-13")
    with pytest.raises(RuntimeError, match="name"):
        subsets.create_subset(project, "../bad", launcher)
    with pytest.raises(RuntimeError, match="subfolder"):
        subsets.safe_branch(project, project)
    error = cli(project, "--end-stage", "12", check=False)
    assert error.returncode != 0 and "separately" in error.stderr


def test_interrupt_stops_stage_before_releasing_lock(tmp_path, monkeypatch):
    launcher = launcher_module()
    pidfile = tmp_path / "worker.pid"
    code = "import os, time; from pathlib import Path; Path({!r}).write_text(str(os.getpid())); print('ready', flush=True); time.sleep(60)".format(str(pidfile))

    class InterruptedTerminal:
        def write(self, text):
            raise KeyboardInterrupt()

    with open(tmp_path / "run.lock", "a") as lock:
        with monkeypatch.context() as patch:
            patch.setattr(sys, "stderr", InterruptedTerminal())
            with pytest.raises(KeyboardInterrupt):
                subsets.stream_stage([sys.executable, "-c", code], launcher.PipelineLogger(tmp_path), tmp_path / "stage.log", lock.fileno())
    with pytest.raises(ProcessLookupError):
        os.kill(int(pidfile.read_text()), 0)
