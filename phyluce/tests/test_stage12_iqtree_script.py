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
    assert "`raw`、`collapsed` 或 `treeshrink`" in readme

    astral_script = (combo / "05_run_astral.sh").read_text()
    assert "Choose gene tree input for ASTRAL" in astral_script
    assert "1. raw" in astral_script
    assert "2. collapsed" in astral_script
    assert "3. treeshrink" in astral_script
    assert 'OUTPUT="astral/${PROJECT_PREFIX}_all_gene_trees.${LABEL}.astral.tre"' in astral_script
    assert 'PROJECT_PREFIX="Acropora"' in astral_script
    assert "astral/Acropora_all_gene_trees.raw.astral.tre" in readme
    assert "ASTRAL_THREADS" not in astral_script
    assert 'ASTRAL_ANNOTATION="${ASTRAL_ANNOTATION:-}"' in astral_script
    assert 'ASTRAL_COMMAND+=("-t" "${ASTRAL_ANNOTATION}")' in astral_script

    assert (combo / "alignments" / "uce-1.fasta").read_text() == (
        ">sampleA\nACGT\n>sampleB\nA-GT\n"
    )
    assert (combo / "alignments" / "uce-2.fasta").read_text() == (
        ">sampleA\nTTAA\n>sampleB\nTT-A\n"
    )

    tree_dir = combo / "iqtree_gene_trees"
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

    assert (combo / "gene_trees" / "gene_tree_files.txt").read_text().splitlines() == [
        "iqtree_gene_trees/uce-1.treefile",
        "iqtree_gene_trees/uce-2.treefile",
    ]
    assert (combo / "gene_trees" / "all_gene_trees.tre").read_text() == (
        "(sampleA,sampleB);\n(sampleA:0.1,sampleB:0.2);\n"
    )

    (tree_dir / "uce-1.treefile").write_text("(sampleA:0.1,sampleB:0.0000")
    invalid = subprocess.run(
        [str(combo / "02_collect_gene_trees.sh")],
        capture_output=True,
        text=True,
    )
    assert invalid.returncode != 0
    assert "Invalid or incomplete Newick" in invalid.stderr
    assert (combo / "gene_trees" / "invalid_gene_trees.txt").read_text() == "uce-1\n"

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

    fake_nw_ed = tmp_path / "fake_nw_ed"
    write_executable(fake_nw_ed, "#!/usr/bin/env bash\ncat \"$1\"\n")
    subprocess.run(
        [str(combo / "03_collapse_low_support_branches.sh"), "25"],
        check=True,
        env={**dict(os.environ), "NW_ED_BIN": str(fake_nw_ed)},
    )
    collapsed = combo / "gene_trees" / "all_gene_trees.BS25collapsed.tre"
    assert collapsed.is_file()
    assert (combo / "gene_trees" / "collapsed_gene_trees_path.txt").read_text() == (
        "gene_trees/all_gene_trees.BS25collapsed.tre\n"
    )

    fake_treeshrink = tmp_path / "fake_treeshrink"
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
cp "$input" "$outdir/$prefix.trees"
: > "$outdir/$prefix.txt"
: > "$outdir/$prefix.log"
""",
    )
    subprocess.run(
        [str(combo / "04_run_treeshrink.sh"), "collapsed"],
        check=True,
        env={
            **dict(os.environ),
            "TREESHRINK_BIN": str(fake_treeshrink),
        },
    )
    assert (combo / "treeshrink" / "treeshrink.trees").read_text() == collapsed.read_text()

    fake_java = tmp_path / "fake_java"
    fake_java_args = tmp_path / "fake_java_args.txt"
    write_executable(
        fake_java,
        """#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$@" > "$FAKE_JAVA_ARGS"
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
    fake_jar = tmp_path / "astral.5.7.8.jar"
    fake_jar.write_text("fixture")
    subprocess.run(
        [str(combo / "05_run_astral.sh"), "treeshrink"],
        check=True,
        env={
            **dict(os.environ),
            "JAVA_BIN": str(fake_java),
            "ASTRAL_JAR": str(fake_jar),
            "JAVA_MEMORY": "1G",
            "FAKE_JAVA_ARGS": str(fake_java_args),
        },
    )
    java_arguments = fake_java_args.read_text().splitlines()
    assert java_arguments[0] == "-Xmx1G"
    assert "-t" not in java_arguments
    assert (combo / "astral" / "Acropora_all_gene_trees.treeshrink.astral.tre").is_file()

    manifest = (output_root / "branch_manifest.tsv").read_text()
    assert "analysis_readme" in manifest
    assert "astral_script" in manifest
    assert str(combo / "05_run_astral.sh") in manifest
