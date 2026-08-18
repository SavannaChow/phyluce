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


def test_parallel_spades_runs_fstrim_after_concurrent_batch(monkeypatch, tmp_path):
    launcher = load_launcher()
    calls = []

    class Logger:
        def info(self, message):
            calls.append(("log", message))

    def fake_run_spades(sample, *args):
        calls.append(("spades", sample))
        return sample == "new_sample"

    def fake_fstrim(spades_tmp_root, logger, log_root):
        calls.append(("fstrim", str(spades_tmp_root)))
        return True

    monkeypatch.setattr(
        launcher,
        "require_sudo_for_fstrim",
        lambda logger, log_root: calls.append(("sudo", str(log_root))),
    )
    monkeypatch.setattr(
        launcher,
        "start_sudo_keepalive",
        lambda logger, log_root: (
            calls.append(("keepalive_start", str(log_root))) or (object(), object())
        ),
    )
    monkeypatch.setattr(
        launcher,
        "stop_sudo_keepalive",
        lambda stop_event, thread, logger: calls.append(("keepalive_stop", None)),
    )
    monkeypatch.setattr(launcher, "run_spades", fake_run_spades)
    monkeypatch.setattr(launcher, "run_fstrim_for_spades_scratch", fake_fstrim)
    monkeypatch.setattr(
        launcher, "write_spades_summary_table", lambda *args, **kwargs: None
    )

    launcher.run_spades_for_samples(
        {"old_sample": None, "new_sample": None},
        tmp_path / "01_Trimmed",
        tmp_path / "02_Spades_assembly",
        "spades.py",
        72,
        32,
        2,
        110,
        "2",
        tmp_path / "scratch",
        True,
        True,
        Logger(),
        tmp_path / "00_Log",
    )

    fstrim_call = ("fstrim", str(tmp_path / "scratch"))
    assert fstrim_call in calls
    assert calls.index(fstrim_call) > max(
        idx for idx, call in enumerate(calls) if call[0] == "spades"
    )
    assert calls[-1] == ("keepalive_stop", None)
