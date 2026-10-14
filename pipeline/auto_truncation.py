"""Choose DADA2 truncation lengths from QIIME 2 per-position quality summaries.

Used by 03_denoise.sh when DADA2_TRUNC_F/R are "auto". Input: the
forward-/reverse-seven-number-summaries.tsv exported from
02_trim/reads_trimmed_summary.qzv (rows: count, 2%, 9%, 25%, 50%, 75%, 91%, 98%;
one column per read position, 1-based).

Rule, per read direction (the nf-core/ampliseq quality rule):
  quality cut = last position before the first one whose median quality is below
                --min-median-q (ampliseq default 25)
  length cut  = last position reached by at least --min-read-fraction of the reads
                (ampliseq default 0.75); DADA2 discards reads shorter than the cut
  cut         = the smaller of the two
Paired reads with --sum-target: if forward + reverse cut is below the target (longest
insert that must merge + 12 nt overlap), lengthen the read whose next position has
the higher median quality, one position at a time, never past its length cut. If the
target cannot be reached, exit with status 2 and say why; nothing is chosen silently.

Prints "FORWARD REVERSE" (or "FORWARD") on stdout and writes every input and decision
to --out as a two-column setting/value table.
"""
import argparse
import csv
import math
import sys


def load(path):
    with open(path, newline='') as fh:
        rows = [r for r in csv.reader(fh, delimiter='\t') if r]
    data = {}
    for row in rows[1:]:
        values = []
        for v in row[1:]:
            try:
                x = float(v)
            except ValueError:
                x = math.nan
            values.append(x)
        data[row[0]] = values
    for key in ('count', '50%'):
        if key not in data:
            sys.exit(f'{path}: no "{key}" row; is this a seven-number summary?')
    return data


def cuts(data, min_q, min_frac):
    count, median = data['count'], data['50%']
    n = len(count)
    total = count[0] if count and not math.isnan(count[0]) else 0
    q_cut = next((i for i, m in enumerate(median) if math.isnan(m) or m < min_q), n)
    len_cut = next((i for i, c in enumerate(count)
                    if math.isnan(c) or c < min_frac * total), n)
    return q_cut, len_cut


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--forward', required=True)
    ap.add_argument('--reverse')
    ap.add_argument('--min-median-q', type=float, default=25)
    ap.add_argument('--min-read-fraction', type=float, default=0.75)
    ap.add_argument('--sum-target', type=int)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    reads = {'forward': load(a.forward)}
    if a.reverse:
        reads['reverse'] = load(a.reverse)
    info, chosen = {}, {}
    for name, data in reads.items():
        q_cut, len_cut = cuts(data, a.min_median_q, a.min_read_fraction)
        info[name] = (q_cut, len_cut)
        chosen[name] = min(q_cut, len_cut)

    note = {name: ('quality cut' if info[name][0] <= info[name][1] else 'length cut')
            for name in reads}
    status, message, state = 0, '', 'ok'
    for name in reads:
        if chosen[name] == 0:
            status, state = 2, 'no qualifying position'
            message = (f'{name} reads: no truncation length qualifies (median quality '
                       f'below {a.min_median_q:g}, or too few reads, from position 1). '
                       f'Inspect the quality plot and set DADA2_TRUNC_F/R yourself.')
            break

    if a.reverse and a.sum_target and not status:
        f, r = reads['forward']['50%'], reads['reverse']['50%']
        while chosen['forward'] + chosen['reverse'] < a.sum_target:
            can_f = chosen['forward'] < info['forward'][1]
            can_r = chosen['reverse'] < info['reverse'][1]
            if not (can_f or can_r):
                break
            # next position (0-based index = current cut) with the better median quality;
            # ties go to the forward read, which is usually the better one
            qf = f[chosen['forward']] if can_f else -math.inf
            qr = r[chosen['reverse']] if can_r else -math.inf
            if qf >= qr:
                chosen['forward'] += 1
                note['forward'] = 'lengthened for overlap'
            else:
                chosen['reverse'] += 1
                note['reverse'] = 'lengthened for overlap'
        total = chosen['forward'] + chosen['reverse']
        if total < a.sum_target:
            status, state = 2, 'overlap target not reachable'
            message = (
                f'Overlap target not reachable: forward + reverse can reach at most '
                f'{total} nt (read-length limits {info["forward"][1]} + '
                f'{info["reverse"][1]}), but DADA2_TRUNC_SUM_TARGET is {a.sum_target}. '
                f'Inserts longer than {total - 12} nt cannot merge. Choose one: lower '
                f'DADA2_TRUNC_SUM_TARGET (accept losing taxa with longer inserts), set '
                f'DADA2_TRUNC_F/R to numbers yourself, or use the forward reads only '
                f'(LAYOUT="single").')

    rows = [('setting', 'value'),
            ('min_median_q', f'{a.min_median_q:g}'),
            ('min_read_fraction', f'{a.min_read_fraction:g}'),
            ('sum_target', a.sum_target if a.sum_target else 'none')]
    for name, data in reads.items():
        c = chosen[name]
        med = data['50%'][c - 1] if c > 0 else math.nan
        frac = data['count'][c - 1] / data['count'][0] if c > 0 and data['count'][0] else 0
        rows += [(f'{name}_quality_cut', info[name][0]),
                 (f'{name}_length_cut', info[name][1]),
                 (f'{name}_chosen', c),
                 (f'{name}_median_q_at_chosen', f'{med:g}'),
                 (f'{name}_reads_reaching_chosen', f'{frac:.4f}'),
                 (f'{name}_reason', note[name])]
    if a.reverse:
        total = chosen['forward'] + chosen['reverse']
        rows += [('sum', total), ('max_mergeable_insert', total - 12)]
    rows.append(('status', state))
    with open(a.out, 'w', newline='') as fh:
        csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(rows)

    if status:
        print(message, file=sys.stderr)
        return status
    print(' '.join(str(chosen[n]) for n in reads))
    return 0


if __name__ == '__main__':
    sys.exit(main())
