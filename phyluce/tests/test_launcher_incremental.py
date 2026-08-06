import importlib.machinery
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]


def load_launcher():
    loader = importlib.machinery.SourceFileLoader(
        "phyluce_launcher_under_test", str(ROOT / "phyluce_launcher")
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def write_project_info(project_dir, project_id="Acropora"):
    project_dir.mkdir(parents=True, exist_ok=True)
    (project_dir / "project.info").write_text(
        "project_id\t{}\n".format(project_id)
    )


def test_stage1_can_reuse_project_with_project_info(tmp_path):
    launcher = load_launcher()
    raw_dir = tmp_path / "RAW"
    raw_dir.mkdir()
    project_dir = tmp_path / "Acropora"
    write_project_info(project_dir)

    project_id, resolved_dir = launcher.resolve_project_for_run(
        tmp_path, "Acropora", 1, raw_dir=raw_dir
    )

    assert project_id == "Acropora"
    assert resolved_dir == project_dir


def test_stage1_can_reuse_legacy_project_with_stage_dirs(tmp_path):
    launcher = load_launcher()
    raw_dir = tmp_path / "RAW"
    raw_dir.mkdir()
    project_dir = tmp_path / "Acropora"
    (project_dir / "01_Trimmed").mkdir(parents=True)

    project_id, resolved_dir = launcher.resolve_project_for_run(
        tmp_path, "Acropora", 1, raw_dir=raw_dir
    )

    assert project_id == "Acropora"
    assert resolved_dir == project_dir


def test_stage1_rejects_unrelated_nonempty_project_dir(tmp_path):
    launcher = load_launcher()
    raw_dir = tmp_path / "RAW"
    raw_dir.mkdir()
    project_dir = tmp_path / "Acropora"
    project_dir.mkdir()
    (project_dir / "notes.txt").write_text("not a launcher project\n")

    with pytest.raises(RuntimeError, match="does not look like a phyluce"):
        launcher.resolve_project_for_run(
            tmp_path, "Acropora", 1, raw_dir=raw_dir
        )


def test_fastp_skips_completed_sample(tmp_path):
    launcher = load_launcher()
    messages = []

    class Logger:
        def info(self, message):
            messages.append(message)

    sample = "sampleA"
    raw_pair = (tmp_path / "sampleA_R1.fastq", tmp_path / "sampleA_R2.fastq")
    for path in raw_pair:
        path.write_text("@read\nACGT\n+\n!!!!\n")
    outputs = launcher.fastp_outputs(tmp_path / "01_Trimmed", sample)
    outputs["sample_dir"].mkdir(parents=True)
    outputs["report_dir"].mkdir()
    outputs["r1"].write_text("@read\nACGT\n+\n!!!!\n")
    outputs["r2"].write_text("@read\nACGT\n+\n!!!!\n")

    ran = launcher.run_fastp(
        sample,
        raw_pair,
        tmp_path / "01_Trimmed",
        "fastp",
        4,
        Logger(),
        tmp_path,
    )

    assert ran is False
    assert messages == [
        "Skipping fastp for sampleA; trimmed FASTQ pair already exists."
    ]
