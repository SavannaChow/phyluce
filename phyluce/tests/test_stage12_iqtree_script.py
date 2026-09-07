import csv
import os
import stat
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def write_locus(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(">{}\n{}\n".format(name, seq))


def write_project_info(path, project_id):
    path.write_text("project_id\t{}\n".format(project_id))


def write_executable(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_stage12_writes_named_iqtree_script_for_supermatrix(tmp_path):
    write_project_info(tmp_path / "project.info", "Acropora")
    stage11_combo = tmp_path / "stage11" / "edge" / "incomplete" / "min_taxa_050"
    alignments = stage11_combo / "alignments"
    alignments.mkdir(parents=True)
    write_locus(
        alignments / "uce-1.fasta",
        [
            ("sampleA|uce-1|contigA|10-20|strand_+", "ACGT"),
            ("_R_sampleB|uce-1|contigB|30-40|strand_-", "A-GT"),
        ],
    )
    write_locus(
        alignments / "uce-2.fasta",
        [
            ("sampleA|uce-2|contigC|50-60|strand_+", "TTAA"),
            ("sampleB|uce-2|contigD|70-80|strand_-", "TT-A"),
        ],
    )

    stage11_manifest = tmp_path / "stage11" / "threshold_manifest.tsv"
    with open(stage11_manifest, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "trim_mode",
                "dataset_label",
                "threshold_group",
                "completeness_threshold",
                "filter_label",
                "output_dir",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "trim_mode": "edge",
                "dataset_label": "incomplete",
                "threshold_group": "incomplete",
                "completeness_threshold": "0.50",
                "filter_label": "min_taxa_050",
                "output_dir": str(stage11_combo),
            }
        )

    aggregated = tmp_path / "not_the_project_name.filtered.loci.tsv"
    with open(aggregated, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"sample_id": "sampleA"})
        writer.writerow({"sample_id": "sampleB"})

    output_root = tmp_path / "stage12"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "stage12_prepare_supermatrix_and_gene_tree_datasets"),
            "--input-dir",
            str(tmp_path / "stage11"),
            "--aggregated-loci",
            str(aggregated),
            "--output-dir",
            str(output_root),
            "--analysis-branches",
            "supermatrix",
        ],
        check=True,
    )

    script = (
        output_root
        / "supermatrix"
        / "edge"
        / "incomplete"
        / "min_taxa_050"
        / "analysis_ready"
        / "raxml_iqtree"
        / "run_iqtree_edge-incomplete-min_taxa_050.sh"
    )
    assert script.is_file()
    assert script.stat().st_mode & stat.S_IXUSR
    assert script.read_text().splitlines()[-1:] == [
        (
            "iqtree -s ./edge-incomplete-min_taxa_050.phylip "
            "-p ./edge-incomplete-min_taxa_050.charsets.nexus "
            "-m MFP+MERGE -rcluster 10 -T AUTO -B 1000 "
            "--prefix ./Acropora_edge-incomplete-min_taxa_050"
        ),
    ]
    assert 'cd "${SCRIPT_DIR}"' in script.read_text()
    assert "watch_iqtree_progress_htop.sh --launch ./Acropora_edge-incomplete-min_taxa_050.log" in script.read_text()
    monitor = script.parent / "watch_iqtree_progress_htop.sh"
    assert monitor.read_bytes() == (ROOT / monitor.name).read_bytes()
    assert monitor.stat().st_mode & stat.S_IXUSR
    subprocess.run(["bash", "-n", str(script)], check=True)

    manifest = (output_root / "branch_manifest.tsv").read_text()
    assert "iqtree_script" in manifest
    assert str(script) in manifest


def test_stage12_writes_astral_workflow_for_gene_tree_branch(tmp_path):
    write_project_info(tmp_path / "project.info", "Acropora")
    stage11_combo = tmp_path / "stage11" / "internal" / "incomplete" / "min_taxa_075"
    alignments = stage11_combo / "alignments"
    alignments.mkdir(parents=True)
    write_locus(
        alignments / "uce-1.fasta",
        [
            ("sampleA|uce-1|contigA|10-20|strand_+", "ACGT"),
            ("_R_sampleB|uce-1|contigB|30-40|strand_-", "A-GT"),
            ("sampleEmpty|uce-1|contigE|90-100|strand_+", "????"),
        ],
    )
    write_locus(
        alignments / "uce-2.fasta",
        [
            ("sampleA|uce-2|contigC|50-60|strand_+", "TTAA"),
            ("sampleB|uce-2|contigD|70-80|strand_-", "TT-A"),
        ],
    )

    stage11_manifest = tmp_path / "stage11" / "threshold_manifest.tsv"
    with open(stage11_manifest, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "trim_mode",
                "dataset_label",
                "threshold_group",
                "completeness_threshold",
                "filter_label",
                "output_dir",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "trim_mode": "internal",
                "dataset_label": "incomplete",
                "threshold_group": "incomplete",
                "completeness_threshold": "0.75",
                "filter_label": "min_taxa_075",
                "output_dir": str(stage11_combo),
            }
        )

    aggregated = tmp_path / "not_the_project_name.filtered.loci.tsv"
    with open(aggregated, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"sample_id": "sampleA"})
        writer.writerow({"sample_id": "sampleB"})

    output_root = tmp_path / "stage12"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "stage12_prepare_supermatrix_and_gene_tree_datasets"),
            "--input-dir",
            str(tmp_path / "stage11"),
            "--aggregated-loci",
            str(aggregated),
            "--output-dir",
            str(output_root),
            "--analysis-branches",
            "gene-tree",
        ],
        check=True,
    )

    combo = (
        output_root
        / "gene-tree"
        / "internal"
        / "incomplete"
        / "min_taxa_075"
    )
    expected_scripts = [
        "watch_iqtree_progress_htop.sh",
        "validate_gene_trees.py",
        "01_run_iqtree_gene_trees.sh",
        "02_collect_gene_trees.sh",
        "03_collapse_low_support_branches.sh",
        "04_run_treeshrink.sh",
        "05_run_astral.sh",
    ]
    for script_name in expected_scripts:
        script = combo / script_name
        assert script.is_file()
        assert script.stat().st_mode & stat.S_IXUSR
        if script.suffix == ".py":
            subprocess.run([sys.executable, "-m", "py_compile", str(script)], check=True)
        else:
            subprocess.run(["bash", "-n", str(script)], check=True)

    readme = (combo / "README_ASTRAL.md").read_text()
    assert "Recommended workflow / 建議流程" in readme
    assert "TreeShrink 會依據 branch" in readme

    astral_script = (combo / "05_run_astral.sh").read_text()
    assert "astral-mp" in astral_script
    assert '-i "04_treeshrink/treeshrink.tre"' in astral_script
    assert '-o "05_astral/astral_species.tree"' in astral_script
    assert "ASTRAL_JAR" not in astral_script

    assert (combo / "alignments" / "uce-1.fasta").read_text() == (
        ">sampleA\nACGT\n>sampleB\nA-GT\n>sampleEmpty\n????\n"
    )
    assert (combo / "alignments" / "uce-2.fasta").read_text() == (
        ">sampleA\nTTAA\n>sampleB\nTT-A\n"
    )

    tree_dir = combo / "01_iqtree_gene_trees"
    tree_dir.mkdir()
    (tree_dir / "uce-1.treefile").write_text(
        "(sampleA|uce-1|contigA:0.1,sampleB|uce-1|contigB:0.2);\n"
    )
    old_labels = subprocess.run(
        [
            sys.executable,
            str(combo / "validate_gene_trees.py"),
            "--single",
            "--alignment",
            str(combo / "alignments" / "uce-1.fasta"),
            str(tree_dir / "uce-1.treefile"),
        ],
        capture_output=True,
        text=True,
    )
    assert old_labels.returncode != 0
    assert "tree/alignment taxon mismatch" in old_labels.stderr

    (tree_dir / "uce-1.treefile").write_text("(sampleA,sampleB);\n")
    (tree_dir / "uce-2.treefile").write_text("(sampleA:0.1,sampleB:0.2);\n")
    subprocess.run([str(combo / "02_collect_gene_trees.sh")], check=True)

    assert (combo / "02_gene_trees" / "all_gene_trees.tre").read_text() == (
        "(sampleA,sampleB);\n(sampleA:0.1,sampleB:0.2);\n"
    )

    (tree_dir / "uce-1.treefile").write_text("(sampleA:0.1,sampleB:0.0000")

    fake_iqtree = tmp_path / "fake_iqtree"
    write_executable(
        fake_iqtree,
        """#!/usr/bin/env bash
set -euo pipefail
prefix=""
while (( $# > 0 )); do
  case "$1" in
    --prefix) prefix="$2"; shift 2 ;;
    *) shift ;;
  esac
done
printf '(sampleA:0.1,sampleB:0.2);\\n' > "${prefix}.treefile"
""",
    )
    subprocess.run(
        [str(combo / "01_run_iqtree_gene_trees.sh")],
        check=True,
        env={
            **dict(os.environ),
            "IQTREE_BIN": str(fake_iqtree),
            "IQTREE_PROGRESS": "0",
            "JOBS": "2",
            "THREADS_PER_JOB": "1",
        },
    )
    assert (tree_dir / "uce-1.treefile").read_text().endswith(");\n")
    archived = list((tree_dir / "invalid_previous").glob("uce-1_*/*"))
    assert any(path.name == "uce-1.treefile" for path in archived)
    subprocess.run([str(combo / "02_collect_gene_trees.sh")], check=True)

    fake_nw_ed = tmp_path / "nw_ed"
    write_executable(fake_nw_ed, "#!/usr/bin/env bash\ncat \"$1\"\n")
    subprocess.run(
        [str(combo / "03_collapse_low_support_branches.sh")],
        check=True,
        env={
            **dict(os.environ),
            "PATH": "{}:{}".format(tmp_path, os.environ.get("PATH", "")),
        },
    )
    collapsed = combo / "03_collapsed_gene_trees" / "all_gene_trees.BS10collapsed.tre"
    assert collapsed.is_file()

    fake_treeshrink = tmp_path / "run_treeshrink.py"
    write_executable(
        fake_treeshrink,
        """#!/usr/bin/env bash
set -euo pipefail
input=""; outdir=""; prefix=""
while (( $# > 0 )); do
  case "$1" in
    -t) input="$2"; shift 2 ;;
    -o) outdir="$2"; shift 2 ;;
    -O) prefix="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$outdir"
cp "$input" "$outdir/$prefix.tre"
: > "$outdir/$prefix.txt"
: > "$outdir/$prefix.log"
""",
    )
    subprocess.run(
        [str(combo / "04_run_treeshrink.sh")],
        check=True,
        env={
            **dict(os.environ),
            "PATH": "{}:{}".format(tmp_path, os.environ.get("PATH", "")),
        },
    )
    assert (combo / "04_treeshrink" / "treeshrink.tre").read_text() == collapsed.read_text()

    fake_astral = tmp_path / "astral-mp"
    fake_astral_args = tmp_path / "fake_astral_args.txt"
    write_executable(
        fake_astral,
        """#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$@" > "$FAKE_ASTRAL_ARGS"
output=""
while (( $# > 0 )); do
  case "$1" in
    -o) output="$2"; shift 2 ;;
    *) shift ;;
  esac
done
    printf '(sampleA,sampleB);\\n' > "$output"
""",
    )
    subprocess.run(
        [str(combo / "05_run_astral.sh")],
        check=True,
        env={
            **dict(os.environ),
            "PATH": "{}:{}".format(tmp_path, os.environ.get("PATH", "")),
            "FAKE_ASTRAL_ARGS": str(fake_astral_args),
        },
    )
    astral_arguments = fake_astral_args.read_text().splitlines()
    assert astral_arguments == [
        "-i",
        "04_treeshrink/treeshrink.tre",
        "-o",
        "05_astral/astral_species.tree",
    ]
    assert (combo / "05_astral" / "astral_species.tree").is_file()

    manifest = (output_root / "branch_manifest.tsv").read_text()
    assert "analysis_readme" in manifest
    assert "astral_script" in manifest
    assert str(combo / "05_run_astral.sh") in manifest
