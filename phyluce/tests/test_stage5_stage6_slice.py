import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def write_hits(path, rows):
    header = (
        "Query name\tQuery length\tQuery start\tQuery end\tRelative strand\t"
        "Target name\tTarget length\tTarget start\tTarget end\tMatches\t"
        "Alignment block\tIdentity\tCoverage\tMapQ\n"
    )
    path.write_text(header + "".join(rows))


def run_stage5(tmp_path, rows, max_locus_span=5000):
    hits = tmp_path / "hits.tsv"
    filtered = tmp_path / "filtered.tsv"
    trace = tmp_path / "trace.tsv"
    loci = tmp_path / "loci.tsv"
    summary = tmp_path / "summary.txt"
    write_hits(hits, rows)
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "stage5_filter_uce_hits_trace"),
            "--input",
            str(hits),
            "--output",
            str(filtered),
            "--trace-output",
            str(trace),
            "--locus-output",
            str(loci),
            "--summary",
            str(summary),
            "--max-locus-span",
            str(max_locus_span),
        ],
        check=True,
    )
    return filtered, trace, loci, summary


def fasta_sequence_length(path):
    sequence = []
    for line in path.read_text().splitlines():
        if not line.startswith(">"):
            sequence.append(line)
    return len("".join(sequence))


def test_stage5_removes_same_contig_locus_with_multiple_hit_clusters(tmp_path):
    rows = [
        "uce-1_p1\t120\t0\t120\t+\tchr1\t20000\t100\t220\t120\t120\t1.0000\t1.0000\t60\n",
        "uce-1_p2\t120\t0\t120\t+\tchr1\t20000\t10000\t10120\t120\t120\t1.0000\t1.0000\t60\n",
    ]

    filtered, trace, loci, summary = run_stage5(tmp_path, rows, max_locus_span=5000)

    assert len(filtered.read_text().splitlines()) == 1
    assert len(loci.read_text().splitlines()) == 1
    assert "uce_matches_multiple_locus_clusters" in trace.read_text()
    assert "chr1:100-220,chr1:10000-10120" in summary.read_text()


def test_stage5_removes_uce_with_hits_on_multiple_contigs(tmp_path):
    rows = [
        "uce-1_p1\t120\t0\t120\t+\tcontig1\t2000\t100\t220\t120\t120\t1.0000\t1.0000\t60\n",
        "uce-1_p2\t120\t0\t120\t+\tcontig2\t2000\t100\t220\t120\t120\t1.0000\t1.0000\t60\n",
    ]

    filtered, trace, loci, summary = run_stage5(tmp_path, rows, max_locus_span=5000)

    assert len(filtered.read_text().splitlines()) == 1
    assert len(loci.read_text().splitlines()) == 1
    assert "uce_matches_multiple_locus_clusters" in trace.read_text()
    assert "contig1:100-220,contig2:100-220" in summary.read_text()


def test_stage6_slice_uses_locus_span_plus_flank(tmp_path):
    rows = [
        "uce-1_p1\t180\t0\t180\t+\tchr1\t5000\t1000\t1180\t180\t180\t1.0000\t1.0000\t60\n",
        "uce-1_p2\t180\t0\t180\t+\tchr1\t5000\t1120\t1300\t180\t180\t1.0000\t1.0000\t60\n",
    ]
    _, _, loci, _ = run_stage5(tmp_path, rows, max_locus_span=5000)
    contigs = tmp_path / "contigs.fasta"
    contigs.write_text(">chr1\n{}\n".format("A" * 5000))

    subprocess.run(
        [
            sys.executable,
            str(ROOT / "stage6_extract_per_sample_uce_fasta_records"),
            "--hits",
            str(loci),
            "--contigs",
            str(contigs),
            "--sample-prefix",
            "sample",
            "--filename-prefix",
            "sample.filtered",
            "--output-dir",
            str(tmp_path),
            "--mode",
            "slice",
            "--flank",
            "500",
        ],
        check=True,
    )

    assert fasta_sequence_length(tmp_path / "sample.filtered.simple.fasta") == 1300
    name_map = (tmp_path / "sample.filtered.name_map.tsv").read_text()
    assert "\t1000\t1300\t500\t1800\tslice\t500\t" in name_map
