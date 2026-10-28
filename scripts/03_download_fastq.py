#!/usr/bin/env python3
"""
Step 3 - Download FASTQ files from SRA.

For every run listed in docs/metadata/samples.tsv:
  1. prefetch       -> download the .sra file
  2. fasterq-dump   -> convert .sra to .fastq
  3. verify         -> number of reads must match the SRA run table (Bases / AvgSpotLen)
  4. pigz           -> compress to .fastq.gz and test the compressed file
  5. clean up       -> delete the .sra file (only if everything above succeeded)

The script is resumable: runs that are already finished are skipped,
so it can be re-run safely after an interruption.

Usage (from anywhere, with the egf-ngs environment active):
    caffeinate -i python3 scripts/03_download_fastq.py --threads 6
"""

import argparse
import csv
import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths: built from the location of this script, NOT from the current folder.
# scripts/03_download_fastq.py -> parent = scripts/ -> parent = project root
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent.parent
SAMPLES_FILE = PROJECT_DIR / "docs" / "metadata" / "samples.tsv"
RUN_TABLE = PROJECT_DIR / "docs" / "metadata" / "GSE69822_SraRunTable.csv"
SRA_DIR = PROJECT_DIR / "data" / "raw" / "sra"
FASTQ_DIR = PROJECT_DIR / "data" / "raw" / "fastq"
TMP_DIR = PROJECT_DIR / "data" / "tmp"
LOG_DIR = PROJECT_DIR / "results" / "logs"

REQUIRED_TOOLS = ["prefetch", "fasterq-dump", "pigz"]


def setup_logging():
    """Write messages both to the screen and to results/logs/03_download_fastq.log."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / "03_download_fastq.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )
    return log_file


def run_command(cmd):
    """Run a shell command; stop with an error if it fails."""
    logging.info("    $ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def read_runs():
    """Return the list of SRR accessions from samples.tsv."""
    with open(SAMPLES_FILE) as f:
        return [row["run"] for row in csv.DictReader(f, delimiter="\t")]


def read_expected_reads():
    """Return {SRR: expected number of reads} from the SRA run table."""
    expected = {}
    with open(RUN_TABLE) as f:
        for row in csv.DictReader(f):
            expected[row["Run"]] = int(row["Bases"]) // int(row["AvgSpotLen"])
    return expected


def count_reads(fastq):
    """Count reads in an uncompressed FASTQ (4 lines per read)."""
    result = subprocess.run(["wc", "-l", str(fastq)], capture_output=True, text=True, check=True)
    n_lines = int(result.stdout.split()[0])
    return n_lines // 4


def process_run(run, expected_reads, threads):
    """Download, convert, verify, compress and clean up one run."""
    fastq = FASTQ_DIR / f"{run}.fastq"
    fastq_gz = FASTQ_DIR / f"{run}.fastq.gz"
    done_flag = FASTQ_DIR / f"{run}.done"   # empty file marking a finished run
    sra_folder = SRA_DIR / run

    # 0. Already finished in a previous execution?
    if done_flag.exists() and fastq_gz.exists():
        logging.info(f"  {run}: already done, skipping")
        return "skipped"

    # 1-2. Reuse an uncompressed FASTQ if it is already complete (e.g. the pilot run)
    if fastq.exists() and count_reads(fastq) == expected_reads:
        logging.info(f"  {run}: complete FASTQ found, skipping download")
    else:
        logging.info(f"  {run}: downloading .sra")
        run_command(["prefetch", run, "--output-directory", str(SRA_DIR)])

        logging.info(f"  {run}: converting to FASTQ")
        run_command([
            "fasterq-dump", str(sra_folder),
            "--outdir", str(FASTQ_DIR),
            "--temp", str(TMP_DIR),
            "--threads", str(threads),
            "--qual-defline", "+",   # shorter '+' line, saves disk space
            "--force",               # overwrite a partial FASTQ from a failed attempt
        ])

        # 3. Verify the number of reads
        n_reads = count_reads(fastq)
        if n_reads != expected_reads:
            raise RuntimeError(f"read count mismatch: got {n_reads}, expected {expected_reads}")
        logging.info(f"  {run}: {n_reads:,} reads (matches run table)")

    # 4. Compress (pigz deletes the .fastq after creating .fastq.gz) and test the archive
    logging.info(f"  {run}: compressing")
    run_command(["pigz", "--force", "-p", str(threads), str(fastq)])
    run_command(["pigz", "--test", str(fastq_gz)])

    # 5. Clean up: only reached if every step above succeeded
    if sra_folder.exists():
        shutil.rmtree(sra_folder)
        logging.info(f"  {run}: .sra deleted")
    done_flag.touch()
    return "ok"


def main():
    parser = argparse.ArgumentParser(description="Download and verify FASTQ files from SRA.")
    parser.add_argument("--threads", type=int, default=6, help="threads for fasterq-dump and pigz (default: 6)")
    args = parser.parse_args()

    log_file = setup_logging()

    # Check that the tools are available (i.e. the egf-ngs environment is active)
    missing = [tool for tool in REQUIRED_TOOLS if shutil.which(tool) is None]
    if missing:
        logging.error(f"Missing tools: {', '.join(missing)}. Did you run 'conda activate egf-ngs'?")
        sys.exit(1)

    for folder in (SRA_DIR, FASTQ_DIR, TMP_DIR):
        folder.mkdir(parents=True, exist_ok=True)

    runs = read_runs()
    expected = read_expected_reads()
    logging.info(f"Starting: {len(runs)} runs, {args.threads} threads. Log: {log_file}")

    results = {}
    for i, run in enumerate(runs, start=1):
        logging.info(f"[{i}/{len(runs)}] {run}")
        start = time.time()
        try:
            status = process_run(run, expected[run], args.threads)
        except Exception as error:   # keep going with the next run if this one fails
            status = "failed"
            logging.error(f"  {run}: FAILED - {error}")
        minutes = (time.time() - start) / 60
        results[run] = status
        logging.info(f"  {run}: {status} ({minutes:.1f} min)")

    # Summary
    ok = [r for r, s in results.items() if s == "ok"]
    skipped = [r for r, s in results.items() if s == "skipped"]
    failed = [r for r, s in results.items() if s == "failed"]
    logging.info("=" * 50)
    logging.info(f"Done: {len(ok)} ok, {len(skipped)} skipped, {len(failed)} failed")
    if failed:
        logging.info("Failed runs (re-run the script to retry): " + ", ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()
