#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import time
from argparse import ArgumentParser
from collections import defaultdict
from shutil import which
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

# ---- Bitwise SAM flag helpers ----
def is_proper_pair(flag: int) -> bool:      return (flag & 0x2)  != 0
def is_first_in_pair(flag: int) -> bool:    return (flag & 0x40) != 0
def is_secondary(flag: int) -> bool:        return (flag & 0x100)!= 0
def is_unpaired_template(flag: int) -> bool:return (flag & 0x1) == 0
def mate_unmapped(flag: int) -> bool:       return (flag & 0x8)  != 0

# ---- CLI ----
def get_parser():
    p = ArgumentParser(
        description=(
            'Process all BAM files in a directory (non-recursive) in parallel using "samtools view -h". '
            'Counts proper pairs once via first-in-pair (FLAG 0x40) on same reference (RNEXT=="="). '
            'Optionally includes singles. Outputs two CSVs: PREFIX.counts.csv and PREFIX.tpm.csv.'
        )
    )
    p.add_argument('--fasta',    metavar='', type=str, required=True, help='FASTA of reference sequences')
    p.add_argument('--bam_dir',  metavar='', type=str, required=True, help='Directory with BAM files (non-recursive)')
    p.add_argument('--output',   metavar='', type=str, required=True, help='Output prefix (.counts.csv and .tpm.csv)')
    p.add_argument('--decimals', metavar='', type=int, default=4,      help='Decimal places for TPM (default 4)')
    p.add_argument('--single',   action='store_true', default=False,   help='Include single reads (default False)')
    p.add_argument('--threads',  metavar='', type=int, default=10,     help='Parallel workers (default 10)')
    p.add_argument('--progress', action='store_true', default=False,   help='Show progress')
    return p

# ---- FASTA lengths (preserve order) ----
def load_lengths(fasta_path, progress=False):
    lengths = {}
    current = None
    with open(fasta_path, 'r') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                current = line[1:].split()[0]
                lengths.setdefault(current, 0)
            else:
                if current is None:
                    print('Error: sequence before any header in FASTA.', file=sys.stderr)
                    sys.exit(1)
                lengths[current] += len(line)
    if progress:
        print(f'Loaded {len(lengths)} reference entries from FASTA.', file=sys.stderr)
    return lengths

# ---- Discover BAMs ----
def find_bams(bam_dir):
    return sorted(
        f for f in (os.path.join(bam_dir, x) for x in os.listdir(bam_dir))
        if os.path.isfile(f) and f.lower().endswith('.bam')
    )

# ---- Worker ----
def process_bam(bam_path, include_single, q, report_every, length_keys):
    """
    Streams alignments via samtools.
    Sends to parent via q:
      ('tick', sample, align_count), ('done', sample), ('error', sample, message)
    Returns (sample_name, counts_dict).
    """
    sample = os.path.splitext(os.path.basename(bam_path))[0]
    counts = defaultdict(int)
    ticker = 0

    cmd = f"samtools view -h {bam_path}"
    try:
        with os.popen(cmd) as pipe:
            for raw in pipe:
                if not raw or raw[0] == '@':
                    continue
                parts = raw.rstrip('\n').split('\t')
                if len(parts) < 7:
                    continue
                try:
                    flag = int(parts[1])
                except ValueError:
                    continue
                rname = parts[2]
                rnext = parts[6]

                if is_secondary(flag):
                    continue

                # Proper pairs counted once via first-in-pair on same reference
                if rname != '*' and is_proper_pair(flag) and is_first_in_pair(flag) and rnext == '=':
                    if rname in length_keys:
                        counts[rname] += 1

                # Optionally include singles
                elif include_single and rname != '*' and (is_unpaired_template(flag) or mate_unmapped(flag)):
                    if rname in length_keys:
                        counts[rname] += 1

                ticker += 1
                if q is not None and report_every and ticker % report_every == 0:
                    try:
                        q.put(('tick', sample, ticker), block=False)
                    except Exception:
                        pass

        if q is not None:
            try:
                q.put(('done', sample), block=False)
            except Exception:
                pass

    except Exception as e:
        if q is not None:
            try:
                q.put(('error', sample, str(e)), block=False)
            except Exception:
                pass
        raise

    return sample, dict(counts)

# ---- Centralized progress: overall bar + live per-file table ----
def run_with_progress(all_bams, args, lengths):
    # Try rich; else fallback
    try:
        from rich.progress import Progress, BarColumn, TextColumn, MofNCompleteColumn, TimeElapsedColumn
        from rich.table import Table
        from rich.live import Live
        from rich.console import Group
        use_rich = True
    except Exception:
        use_rich = False

    manager = mp.Manager()
    q = manager.Queue() if args.progress else None

    results = []
    length_keys = set(lengths.keys())
    ctx = mp.get_context('spawn')

    if use_rich and args.progress:
        overall = Progress(
            TextColumn("[bold]Overall[/bold]"),
            MofNCompleteColumn(),
            BarColumn(),
            TimeElapsedColumn(),
            transient=False,
            redirect_stdout=False,
            redirect_stderr=False,
        )
        task_overall = overall.add_task("overall", total=len(all_bams))

        per_file_counts = {}   # sample -> last alignment count
        active = set()         # currently running samples

        def render_table():
            table = Table(title="Active files", show_lines=False)
            table.add_column("Sample", no_wrap=True)
            table.add_column("Alignments processed", justify="right", no_wrap=True)
            for s in sorted(active):
                n = per_file_counts.get(s, 0)
                table.add_row(s, f"{n//1000}K")
            return table

        with Live(Group(overall, render_table()), refresh_per_second=10,
                  transient=False, redirect_stdout=False, redirect_stderr=False) as live:

            def refresh():
                live.update(Group(overall, render_table()))

            with ProcessPoolExecutor(max_workers=args.threads, mp_context=ctx) as ex:
                futures = {ex.submit(process_bam, bam, args.single, q, 200_000, length_keys): bam for bam in all_bams}
                pending = set(futures.keys())

                while pending:
                    # Drain queue
                    if q is not None:
                        for _ in range(512):
                            try:
                                msg = q.get_nowait()
                            except Exception:
                                break
                            kind = msg[0]
                            if kind == 'tick':
                                _, sample, n = msg
                                active.add(sample)
                                per_file_counts[sample] = n
                            elif kind == 'done':
                                _, sample = msg
                                active.discard(sample)
                            elif kind == 'error':
                                _, sample, emsg = msg
                                active.discard(sample)
                                print(f"[{sample}] error: {emsg}", file=sys.stderr)
                        refresh()

                    # Collect finished futures
                    done_now = [f for f in pending if f.done()]
                    if not done_now:
                        time.sleep(0.05)
                        continue
                    for f in done_now:
                        sample, counts = f.result()
                        results.append((sample, counts))
                        overall.update(task_overall, advance=1)
                        pending.remove(f)
                        active.discard(sample)
                        refresh()
    else:
        # Fallback: overall prints only
        if args.progress:
            print("rich not available; using simple overall progress.", file=sys.stderr)
        with ProcessPoolExecutor(max_workers=args.threads, mp_context=ctx) as ex:
            futures = {ex.submit(process_bam, bam, args.single, None, 0, set(lengths.keys())): bam for bam in all_bams}
            total = len(futures)
            done = 0
            for fut in as_completed(futures):
                sample, counts = fut.result()
                results.append((sample, counts))
                done += 1
                if args.progress:
                    print(f"[Overall] {done}/{total} files complete.", file=sys.stderr)

    return results

# ---- Write CSVs ----
def write_counts_csv(path, gene_ids, results):
    with open(path, 'w') as out:
        out.write('GeneID')
        for sample, _ in results:
            out.write(',' + sample)
        out.write('\n')
        for gid in gene_ids:
            out.write(gid)
            for _, counts in results:
                out.write(',' + str(counts.get(gid, 0)))
            out.write('\n')

def write_tpm_csv(path, gene_ids, results, lengths, decimals):
    L = {k: float(v) for k, v in lengths.items()}
    total_rpk = {}
    rpks_by_sample = {}
    for sample, counts in results:
        trpk = 0.0
        rpks = {}
        for gid in gene_ids:
            c = counts.get(gid, 0)
            l = L.get(gid, 0.0)
            r = (c / l) if l > 0 else 0.0
            rpks[gid] = r
            trpk += r
        rpks_by_sample[sample] = rpks
        total_rpk[sample] = trpk
    scalers = {s: (total_rpk[s] / 1e6) if total_rpk[s] > 0 else 1.0 for s, _ in results}

    with open(path, 'w') as out:
        out.write('GeneID')
        for sample, _ in results:
            out.write(',' + sample)
        out.write('\n')
        for gid in gene_ids:
            out.write(gid)
            for sample, _ in results:
                tpm = (rpks_by_sample[sample][gid] / scalers[sample]) if scalers[sample] > 0 else 0.0
                out.write(',' + str(round(tpm, decimals)))
            out.write('\n')

# ---- Main ----
def main():
    parser = get_parser()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    if which('samtools') is None:
        print('Error: samtools is not in PATH. Exiting.', file=sys.stderr)
        sys.exit(1)

    lengths = load_lengths(args.fasta, progress=args.progress)
    gene_ids_in_order = list(lengths.keys())

    try:
        all_bams = find_bams(args.bam_dir)
    except FileNotFoundError:
        print(f'Error: BAM directory not found: {args.bam_dir}', file=sys.stderr)
        sys.exit(1)
    if not all_bams:
        print('No BAM files found in the specified directory.', file=sys.stderr)
        sys.exit(1)

    if args.progress:
        print(f'Found {len(all_bams)} BAM files. Using {args.threads} workers.', file=sys.stderr)

    results_pairs = run_with_progress(all_bams, args, lengths)
    results_pairs.sort(key=lambda x: x[0])

    counts_path = args.output + '.counts.csv'
    tpm_path    = args.output + '.tpm.csv'

    if args.progress:
        print(f'Writing counts to {counts_path}', file=sys.stderr)
    write_counts_csv(counts_path, gene_ids_in_order, results_pairs)

    if args.progress:
        print(f'Writing TPMs to {tpm_path}', file=sys.stderr)
    write_tpm_csv(tpm_path, gene_ids_in_order, results_pairs, lengths, args.decimals)

    if args.progress:
        print('All done.', file=sys.stderr)

if __name__ == '__main__':
    main()
