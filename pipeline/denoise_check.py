"""Flag samples that DADA2 handled unusually badly (stage 03, after denoising).

Input: stats.tsv exported from denoising_stats.qza (columns input, filtered, denoised,
merged [paired only], non-chimeric). For each step a sample is flagged when it is both
  - more than --mads robust standard deviations (1.4826 x MAD) worse than the median, and
  - more than --min-points percentage points worse than the median,
so steps where nearly every sample is identical (MAD ~ 0) do not flag trivial differences.
  low_filter_pass   % of input reads passing the quality filter
  low_merge         % of denoised reads that merged (paired-end only)
  high_chimera      % of merged (or denoised) reads removed as chimeras
  low_depth         final reads below --min-depth (when given)
Writes a per-sample table to --out and prints a summary. Flags are warnings: nothing is
removed here (MIN_SAMPLE_DEPTH removes shallow samples in stage 06).
"""
import argparse
import csv
import statistics
import sys


def robust(values):
    med = statistics.median(values)
    mad = statistics.median([abs(v - med) for v in values])
    return med, 1.4826 * mad


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--stats', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--mads', type=float, default=4)
    ap.add_argument('--min-points', type=float, default=5)
    ap.add_argument('--min-depth', type=float)
    a = ap.parse_args()

    with open(a.stats, newline='') as fh:
        rows = [r for r in csv.DictReader(fh, delimiter='\t')
                if not r[next(iter(r))].startswith('#q2:')]
    if not rows:
        sys.exit(f'{a.stats}: no samples')
    sid = next(iter(rows[0]))
    paired = 'merged' in rows[0]
    recs = []
    for r in rows:
        inp, filt, den = (float(r[k]) for k in ('input', 'filtered', 'denoised'))
        joined = float(r['merged']) if paired else den
        final = float(r['non-chimeric'])
        recs.append({
            'sample': r[sid], 'input': inp, 'final_reads': final,
            'pass_filter_pct': 100 * filt / inp if inp else 0.0,
            'merged_of_denoised_pct': 100 * joined / den if (paired and den) else None,
            'chimera_pct': 100 * (1 - final / joined) if joined else 0.0})

    checks = [('low_filter_pass', 'pass_filter_pct', -1), ('high_chimera', 'chimera_pct', +1)]
    if paired:
        checks.insert(1, ('low_merge', 'merged_of_denoised_pct', -1))
    limits = {}
    for flag, col, sign in checks:
        med, rsd = robust([r[col] for r in recs])
        dist = max(a.mads * rsd, a.min_points)
        limits[flag] = (col, med, med + sign * dist, sign)
    for r in recs:
        flags = [flag for flag, (col, med, lim, sign) in limits.items()
                 if (r[col] - lim) * sign > 0]
        if a.min_depth is not None and r['final_reads'] < a.min_depth:
            flags.append('low_depth')
        r['flags'] = ','.join(flags)

    cols = ['sample', 'input', 'pass_filter_pct', 'merged_of_denoised_pct', 'chimera_pct',
            'final_reads', 'flags']
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(cols)
        for r in recs:
            w.writerow([r['sample'], f"{r['input']:.0f}", f"{r['pass_filter_pct']:.2f}",
                        '' if r['merged_of_denoised_pct'] is None else f"{r['merged_of_denoised_pct']:.2f}",
                        f"{r['chimera_pct']:.2f}", f"{r['final_reads']:.0f}", r['flags']])

    print(f'{len(recs)} samples; flagging rule: worse than median by > {a.mads:g} robust SD '
          f'and > {a.min_points:g} points')
    for flag, (col, med, lim, sign) in limits.items():
        n = sum(flag in r['flags'].split(',') for r in recs)
        print(f'  {flag:16s} median {med:6.2f}%  limit {lim:6.2f}%  flagged {n}')
    if a.min_depth is not None:
        n = sum('low_depth' in r['flags'].split(',') for r in recs)
        print(f'  {"low_depth":16s} final reads < {a.min_depth:g}  flagged {n}')
    for r in sorted((r for r in recs if r['flags']), key=lambda r: r['pass_filter_pct']):
        print(f"  {r['sample']}: {r['flags']} (passed {r['pass_filter_pct']:.1f}%, "
              f"final {r['final_reads']:.0f} reads)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
