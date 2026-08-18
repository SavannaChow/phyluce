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
