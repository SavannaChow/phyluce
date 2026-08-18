import csv
import stat
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def write_locus(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(">{}\n{}\n".format(name, seq))


def test_stage12_writes_named_iqtree_script_for_supermatrix(tmp_path):
    stage11_combo = tmp_path / "stage11" / "edge" / "incomplete" / "min_taxa_050"
    alignments = stage11_combo / "alignments"
    alignments.mkdir(parents=True)
    write_locus(
        alignments / "uce-1.fasta",
        [("sampleA", "ACGT"), ("sampleB", "A-GT")],
    )
    write_locus(
        alignments / "uce-2.fasta",
        [("sampleA", "TTAA"), ("sampleB", "TT-A")],
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

    aggregated = tmp_path / "aggregated.tsv"
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
    assert script.read_text().splitlines() == [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        (
            "iqtree -s ./edge-incomplete-min_taxa_050.phylip "
            "-p ./edge-incomplete-min_taxa_050.charsets.nexus "
            "-m MFP+MERGE -rcluster 10 -T AUTO -B 1000"
        ),
    ]

    manifest = (output_root / "branch_manifest.tsv").read_text()
    assert "iqtree_script" in manifest
    assert str(script) in manifest


def test_stage12_writes_astral_workflow_for_gene_tree_branch(tmp_path):
    stage11_combo = tmp_path / "stage11" / "internal" / "incomplete" / "min_taxa_075"
    alignments = stage11_combo / "alignments"
    alignments.mkdir(parents=True)
    write_locus(
        alignments / "uce-1.fasta",
        [("sampleA", "ACGT"), ("sampleB", "A-GT")],
    )
    write_locus(
        alignments / "uce-2.fasta",
        [("sampleA", "TTAA"), ("sampleB", "TT-A")],
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

    aggregated = tmp_path / "aggregated.tsv"
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

    readme = (combo / "README_ASTRAL.md").read_text()
    assert "Recommended workflow / 建議流程" in readme
    assert "TreeShrink 會依據 branch" in readme
    assert "`raw`、`collapsed` 或 `treeshrink`" in readme

    astral_script = (combo / "05_run_astral.sh").read_text()
    assert "Choose gene tree input for ASTRAL" in astral_script
    assert "1. raw" in astral_script
    assert "2. collapsed" in astral_script
    assert "3. treeshrink" in astral_script

    tree_dir = combo / "iqtree_gene_trees"
    tree_dir.mkdir()
    (tree_dir / "uce-1.treefile").write_text("(sampleA,sampleB);\n")
    (tree_dir / "uce-2.treefile").write_text("((sampleA,sampleB));\n")
    subprocess.run([str(combo / "02_collect_gene_trees.sh")], check=True)

    assert (combo / "gene_trees" / "gene_tree_files.txt").read_text().splitlines() == [
        "iqtree_gene_trees/uce-1.treefile",
        "iqtree_gene_trees/uce-2.treefile",
    ]
    assert (combo / "gene_trees" / "all_gene_trees.tre").read_text() == (
        "(sampleA,sampleB);\n\n((sampleA,sampleB));\n\n"
    )

    manifest = (output_root / "branch_manifest.tsv").read_text()
    assert "analysis_readme" in manifest
    assert "astral_script" in manifest
    assert str(combo / "05_run_astral.sh") in manifest
