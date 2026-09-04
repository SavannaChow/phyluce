"""Independent, named stage 7-12 runs over a shared stage 5/6 source.

The launcher dispatches here before opening its normal project logs or writing
project.info. Only this module's 13_Subset_Analyses directory is writable.
"""

import csv
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys


ROOT_NAME = "13_Subset_Analyses"
RECORDS = "analysis_records"
EDITABLE = "original_samples.txt"
USED = "analysis_samples.txt"
EXCLUDED = "excluded_samples_in_this_run.txt"
INFO = "analysis_info.json"
SETTINGS_KEYS = (
    "threads", "dataset_modes", "align_datasets", "align_ambiguous",
    "align_matrix_mode", "mafft_path", "mafft_threads_per_job", "max_mafft_jobs",
    "trim_window", "trim_proportion", "trim_threshold", "trim_max_divergence",
    "trim_min_length", "trim_mode", "gblocks_path", "gblocks_b1", "gblocks_b2",
    "gblocks_b3", "gblocks_b4", "post_align_thresholds", "analysis_branches",
    "missing_character",
)


def now():
    # Fixed UTC+8 works on compute hosts without a timezone database as well.
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8)))


def read_sample_list(path):
    samples = []
    seen = set()
    with open(path, encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            sample = line.strip()
            if not sample:
                continue
            if any(char.isspace() for char in sample) or "|" in sample:
                raise RuntimeError("Invalid sample ID on line {} in {}: {}".format(number, path, sample))
            if sample in seen:
                raise RuntimeError("Duplicate sample ID in {}: {}".format(path, sample))
            seen.add(sample)
            samples.append(sample)
    if not samples:
        raise RuntimeError("Sample list is empty: {}".format(path))
    return samples


def write_samples(path, samples):
    Path(path).write_text("".join(sample + "\n" for sample in samples), encoding="utf-8")


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def save_info(branch, info):
    info["updated_at"] = now().isoformat()
    path = branch / RECORDS / INFO
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_info(branch):
    return json.loads((branch / RECORDS / INFO).read_text(encoding="utf-8"))


def subset_root(project):
    root = project / ROOT_NAME
    if root.is_symlink():
        raise RuntimeError("Subset output directory cannot be a symlink: {}".format(root))
    return root


def safe_branch(project, value):
    root = subset_root(project).resolve()
    value = Path(value).expanduser()
    branch = value if value.is_absolute() or value.is_dir() else root / value
    if branch.is_symlink() or branch.resolve().parent != root:
        raise RuntimeError("Select a direct subfolder of {}".format(root))
    branch = branch.resolve()
    # Never let downstream cleanup/copy operations follow output symlinks.
    if any(path.is_symlink() for path in branch.rglob("*")):
        raise RuntimeError("Subset output contains a symlink: {}".format(branch))
    info = load_info(branch)
    if Path(info["source_project"]).resolve() != project:
        raise RuntimeError("Subset belongs to another source project.")
    return branch


def source_catalog(project, launcher):
    stage_dirs = launcher.get_stage_dirs(project)
    source = stage_dirs[launcher.STAGE_NAMES[6]]
    catalog = {}
    for fasta in sorted(source.glob("*/*.filtered.fasta")):
        sample = fasta.parent.name
        if fasta.name != sample + ".filtered.fasta" or not fasta.stat().st_size:
            continue
        loci = stage_dirs[launcher.STAGE_NAMES[5]] / sample / (sample + ".filtered.loci.tsv")
        if not loci.is_file():
            raise RuntimeError("Missing stage 5 loci records for {}: {}".format(sample, loci))
        with open(fasta) as handle:
            count = sum(line.startswith(">") for line in handle)
        if count:
            catalog[sample] = {
                "fasta": str(fasta), "fasta_sha256": digest(fasta),
                "loci": str(loci), "loci_sha256": digest(loci), "loci_count": count,
            }
    if not catalog:
        raise RuntimeError("No usable stage 6 sample FASTAs found in {}".format(source))
    return catalog


def validate_selected(info, samples):
    unknown = set(samples) - set(info["sources"])
    if unknown:
        raise RuntimeError("Unknown sample IDs: {}".format(", ".join(sorted(unknown))))
    if len(samples) < 3:
        raise RuntimeError("Keep at least 3 samples for the locus alignment workflow.")
    for sample in samples:
        source = info["sources"][sample]
        for key in ("fasta", "loci"):
            path = Path(source[key])
            if not path.is_file() or digest(path) != source[key + "_sha256"]:
                raise RuntimeError("Original source changed or missing for {}: {}. Create a new subset.".format(sample, path))
        # Stage 7 uses hits, stage 8 uses extracted sequences. Require matching
        # locus sets so missing sequence records cannot change the denominator.
        loci = set()
        current = None
        has_sequence = False
        with open(source["fasta"]) as handle:
            for line in handle:
                if line.startswith(">"):
                    if current is not None and not has_sequence:
                        raise RuntimeError("Empty sequence for {} in {}".format(current, source["fasta"]))
                    parts = line[1:].split()[0].split("|")
                    if len(parts) < 2 or parts[0] != sample or parts[1] in loci:
                        raise RuntimeError("Invalid or duplicate sample/locus header in {}".format(source["fasta"]))
                    current = parts[1]
                    loci.add(current)
                    has_sequence = False
                elif line.strip():
                    has_sequence = True
        if not loci or not has_sequence:
            raise RuntimeError("Empty extracted sequence for {}".format(sample))
        with open(source["loci"], newline="") as handle:
            expected = {row["UCE locus"] for row in csv.DictReader(handle, delimiter="\t")}
        if loci != expected:
            raise RuntimeError("Stage 5/6 locus mismatch for {}; check the original extraction.".format(sample))


def refresh_index(project):
    root = subset_root(project)
    fields = ["folder", "name", "created_at", "status", "samples", "excluded",
              "last_completed_stage", "post_align_thresholds"]
    # Atomic replace; metadata files remain authoritative for concurrent runs.
    temporary = root / (".analysis_index.{}.tmp".format(os.getpid()))
    with open(temporary, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for path in sorted(root.glob("*/" + RECORDS + "/" + INFO)):
            info = json.loads(path.read_text(encoding="utf-8"))
            writer.writerow({
                "folder": path.parent.parent.name, "name": info["name"],
                "created_at": info["created_at"], "status": info["status"],
                "samples": len(info["analysis_samples"]) if "analysis_samples" in info else "",
                "excluded": len(info["sources"]) - len(info["analysis_samples"]) if "analysis_samples" in info else "",
                "last_completed_stage": info["last_completed_stage"],
                "post_align_thresholds": info.get("settings", {}).get("post_align_thresholds", ""),
            })
    temporary.replace(root / "analysis_index.tsv")


def create_subset(project, name, launcher, copy_from=None, dry_run=False):
    name = (name or "").strip()
    if not name or name in (".", "..") or re.search(r'[/\\\x00-\x1f\x7f]', name):
        raise RuntimeError("Enter a nonempty analysis name without path separators or control characters.")
    catalog = source_catalog(project, launcher)
    samples = sorted(catalog)
    inherited_settings = {}
    if copy_from:
        previous = load_info(copy_from)
        samples = previous.get("analysis_samples") or read_sample_list(copy_from / EDITABLE)
        if set(samples) - set(catalog):
            raise RuntimeError("Copied sample list contains samples absent from the original source.")
        inherited_settings = previous.get("settings", previous.get("default_settings", {}))
    root = subset_root(project)
    stamp = now()
    base = root / (name + "_" + stamp.strftime("%Y%m%d_%H%M%S"))
    if dry_run:
        print("Would create {} with {} samples; no files written.".format(base, len(samples)))
        return base
    root.mkdir(exist_ok=True)
    branch = base
    counter = 1
    while True:
        try:
            branch.mkdir()
            break
        except FileExistsError:
            counter += 1
            branch = root / (base.name + "_{:02d}".format(counter))
    (branch / RECORDS).mkdir()
    source_info = launcher.read_project_info(project) or {}
    info = {
        "name": name, "created_at": stamp.isoformat(), "source_project": str(project),
        "source_project_id": source_info.get("project_id", project.name),
        "status": "draft", "last_completed_stage": 6, "sources": catalog,
        "default_settings": inherited_settings or source_info,
    }
    save_info(branch, info)
    write_samples(branch / EDITABLE, samples)
    with open(branch / RECORDS / "all_samples.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "extracted_loci"])
        writer.writerows((sample, item["loci_count"]) for sample, item in catalog.items())
    refresh_index(project)
    print("Created: {}\nEdit {}: delete the sample lines you want to exclude.\nThen choose Stage 13 -> run and select this folder.".format(branch, branch / EDITABLE))
    return branch


def build_commands(project, branch, info, settings, launcher):
    s = settings
    scripts = Path(launcher.__file__).resolve().parent
    dirs = {number: branch / launcher.STAGE_NAMES[number] for number in range(7, 13)}
    original = launcher.get_stage_dirs(project)
    prefix = info["source_project_id"]
    script = lambda name: scripts / name
    commands = {
        7: launcher.build_aggregate_command(script("stage7_get_complete_incomplete_dataset_definitions"), original[launcher.STAGE_NAMES[5]], dirs[7], prefix)
           + ["--samples-file", str(branch / RECORDS / USED)],
        8: launcher.build_export_locus_fastas_command(script("stage8_get_complete_incomplete_dataset_fastas"), dirs[7], original[launcher.STAGE_NAMES[6]], dirs[8], prefix, s["dataset_modes"])
           + ["--samples-file", str(branch / RECORDS / USED)],
        9: launcher.build_align_locus_fastas_command(script("stage9_align_locus_fastas_with_mafft"), dirs[8], dirs[9], s["align_datasets"], s["max_mafft_jobs"], s["mafft_threads_per_job"], s["align_ambiguous"], s["align_matrix_mode"], s["mafft_path"]),
        10: launcher.build_trim_locus_alignments_command(script("stage10_trim_alignments_with_edge_internal_methods"), dirs[9], dirs[10], s["align_datasets"], s["max_mafft_jobs"], s["trim_window"], s["trim_proportion"], s["trim_threshold"], s["trim_max_divergence"], s["trim_min_length"], s["trim_mode"], s["gblocks_path"], s["gblocks_b1"], s["gblocks_b2"], s["gblocks_b3"], s["gblocks_b4"]),
        11: launcher.build_filter_locus_alignments_command(script("stage11_get_loci_with_minimum_taxa_completeness"), dirs[10], dirs[9], dirs[11], s["align_datasets"], s["post_align_thresholds"]),
        12: launcher.build_add_missing_data_command(script("stage12_prepare_supermatrix_and_gene_tree_datasets"), dirs[11], dirs[7] / (prefix + ".filtered.loci.tsv"), dirs[12], s["align_datasets"], s["post_align_thresholds"], s["analysis_branches"], s["missing_character"]),
    }
    return dirs, commands


def check_executable(value):
    found = shutil.which(value)
    if not found:
        raise RuntimeError("Required executable not found: {}".format(value))
    return str(Path(found).resolve())


def report_final_samples(branch, requested, launcher):
    stage12 = branch / launcher.STAGE_NAMES[12]
    with open(stage12 / "branch_manifest.tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with open(branch / RECORDS / "final_sample_presence.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["analysis_branch", "threshold_group", "filter_label", "sample_id", "present"])
        for row in rows:
            present = set()
            for fasta in (Path(row["output_dir"]) / "alignments").glob("*.fasta"):
                with open(fasta) as sequences:
                    for line in sequences:
                        if line.startswith(">"):
                            sample = line[1:].split()[0].split("|")[0]
                            while sample.startswith("_R_"):
                                sample = sample[3:]
                            present.add(sample)
            missing = set(requested) - present
            if missing:
                print("No retained loci for {} / {}: {}".format(row["threshold_group"], row["filter_label"], ", ".join(sorted(missing))))
            for sample in requested:
                writer.writerow([row["analysis_branch"], row["threshold_group"], row["filter_label"], sample, int(sample in present)])


def stream_stage(command, logger, log_path, lock_fd):
    """Stop the entire alignment worker group before releasing the run lock."""
    command = [str(part) for part in command]
    logger.command(command, log_path=log_path)
    # The stage process also holds the lock. If the launcher is killed outright,
    # a second run cannot delete outputs while that stage is still writing.
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, bufsize=1, start_new_session=True,
                          pass_fds=(lock_fd,)) as process:
        try:
            with open(log_path, "a") as log:
                for line in process.stdout:
                    sys.stderr.write(line)
                    log.write(line)
                    log.flush()
            returncode = process.wait()
            logger.command_status(returncode, log_path=log_path)
            if returncode:
                raise RuntimeError("Stage command failed; see {}".format(log_path))
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise


def run_subset(project, branch, args, launcher):
    info = load_info(branch)
    samples = sorted(read_sample_list(branch / EDITABLE))
    validate_selected(info, samples)
    if "analysis_samples" in info:
        if samples != info["analysis_samples"] or sorted(read_sample_list(branch / RECORDS / USED)) != samples:
            raise RuntimeError("Sample list changed after analysis started. Copy this subset to a new named branch.")
        if info["status"] == "prepared":
            print("This subset is already prepared. Tree-running scripts are in its stage 12 folder; create/copy a branch for a different analysis.")
            return
        for key in SETTINGS_KEYS:
            if launcher.cli_flag_present("--" + key.replace("_", "-")):
                raise RuntimeError("Resume uses saved settings. Create/copy a branch to change parameters.")
        settings = info["settings"]
    else:
        defaults = info.get("default_settings", {})
        resolved = launcher.resolve_pipeline_settings(args, defaults, set(range(7, 13)), project.parent, branch)
        settings = {key: resolved[key] for key in SETTINGS_KEYS}
        settings["mafft_path"] = check_executable(settings["mafft_path"])
        if settings["trim_mode"] in ("internal", "both"):
            settings["gblocks_path"] = check_executable(settings["gblocks_path"])
    dirs, commands = build_commands(project, branch, info, settings, launcher)
    print("Subset: {}\nOriginal: {} samples; analyze: {}; exclude: {}".format(branch.name, len(info["sources"]), len(samples), len(info["sources"]) - len(samples)))
    if args.dry_run:
        for stage, command in commands.items():
            if stage > info["last_completed_stage"]:
                print("Stage {}: {}".format(stage, launcher.format_command(command)))
        return
    # A leftover lock file is harmless; the OS lock lives only while its owning
    # launcher/stage processes are alive.
    with open(branch / RECORDS / "run.lock", "a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("This subset is already running.")
        if load_info(branch) != info:
            raise RuntimeError("Subset state changed while preparing this run; select it again.")
        log_root = branch / "logs"
        log_root.mkdir(exist_ok=True)
        logger = launcher.PipelineLogger(log_root)
        if "analysis_samples" not in info:
            info["analysis_samples"] = samples
            info["settings"] = settings
            write_samples(branch / RECORDS / USED, samples)
            write_samples(branch / RECORDS / EXCLUDED, sorted(set(info["sources"]) - set(samples)))
            # A branch-specific identity keeps exported final tree names distinct.
            tree_prefix = re.sub(r"[^\w.-]+", "_", info["source_project_id"] + "_" + branch.name)
            launcher.write_project_info(branch, {"project_id": tree_prefix})
        info["status"] = "running"
        info.pop("error", None)
        save_info(branch, info)
        refresh_index(project)
        try:
            for stage, command in commands.items():
                if stage <= info["last_completed_stage"]:
                    if not dirs[stage].is_dir():
                        raise RuntimeError("Completed stage output missing: {}".format(dirs[stage]))
                    continue
                # Only unfinished stage outputs are rebuilt. This also prevents
                # stale Stage 12 FASTAs after an interrupted preparation.
                if dirs[stage].exists():
                    shutil.rmtree(dirs[stage])
                dirs[stage].mkdir()
                logger.info("Subset {}: running stage {}".format(branch.name, stage))
                stream_stage(command, logger, log_root / ("stage{:02d}.log".format(stage)), lock.fileno())
                info["last_completed_stage"] = stage
                save_info(branch, info)
                refresh_index(project)
            report_final_samples(branch, samples, launcher)
            info["status"] = "prepared"
            save_info(branch, info)
        except BaseException as exc:
            info["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
            info["error"] = str(exc)
            save_info(branch, info)
            raise
        finally:
            refresh_index(project)
    print("Stage 7-12 preparation complete: {}\nRun the generated IQ-TREE/ASTRAL scripts inside this branch's stage 12 folder to infer trees.".format(branch))


def choose_subset(project):
    paths = sorted(subset_root(project).glob("*/" + RECORDS + "/" + INFO))
    if not paths:
        raise RuntimeError("No subsets exist yet. Choose create first.")
    for number, path in enumerate(paths, 1):
        print("{}. {} [{}]".format(number, path.parent.parent.name, load_info(path.parent.parent)["status"]))
    while True:
        value = input("Select subset number: ").strip()
        if value.isdigit() and 1 <= int(value) <= len(paths):
            return paths[int(value) - 1].parent.parent
        print("Enter one of the listed numbers.")


def manage_subsets(args, launcher):
    _, project = launcher.resolve_project_path(Path.cwd(), args.project_id)
    project = project.resolve()
    if not project.is_dir() or not (project / launcher.STAGE_NAMES[6]).is_dir():
        raise RuntimeError("Select an original project with completed stage 6 outputs.")
    if ROOT_NAME in project.parts:
        raise RuntimeError("Select the original project, not an existing subset.")
    action = args.subset_action
    if action is None:
        if not sys.stdin.isatty():
            raise RuntimeError("Stage 13 requires --subset-action create/run/copy/list.")
        print("Stage 13: create = new full sample list; run = read edited list; copy = new branch from an existing list; list = show branches")
        action = launcher.prompt_choice_with_default("Subset action", "create", ["create", "run", "copy", "list"])
    branch = None
    if action in ("run", "copy"):
        value = args.subset_dir
        if not value:
            if not sys.stdin.isatty():
                raise RuntimeError("Specify --subset-dir for run/copy.")
            value = choose_subset(project)
        branch = safe_branch(project, value)
    if action in ("create", "copy"):
        name = args.subset_name
        if name is None and sys.stdin.isatty():
            name = input("Enter analysis name (required; Taiwan date/time is appended): ")
        create_subset(project, name, launcher, copy_from=branch, dry_run=args.dry_run)
    elif action == "list":
        for path in sorted(subset_root(project).glob("*/" + RECORDS + "/" + INFO)):
            print("{} [{}]".format(path.parent.parent.name, load_info(path.parent.parent)["status"]))
    else:
        run_subset(project, branch, args, launcher)
