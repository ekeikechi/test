"""Write synthetic forward-/reverse-seven-number-summaries.tsv with known answers.

Usage: python3 tests/make_seven_number.py OUTDIR PROFILE
Profiles (positions are after primer removal; 10,000 sampled reads):
  2x300        forward 283 nt, median < 25 from 283  -> quality cut 282, length cut 283
               reverse 281 nt (57% reach 281), median < 25 from 251
                                                     -> quality cut 250, length cut 280
  early_drop   as 2x300, but reverse median < 25 from 201 -> quality cut 200
  2x250        forward 233 nt, reverse 230 nt: F+R can never reach 520
  bad_start    as 2x300, but reverse median below 25 from position 1
"""
import os
import sys

QUANTILES = ['2%', '9%', '25%', '50%', '75%', '91%', '98%']
OFFSETS = [-20, -8, -2, 0, 1, 2, 2]   # other quantiles relative to the median


def write(path, length, frac_at_end, median, n=10000, tail=12):
    positions = list(range(1, length + tail + 1))
    count = [n if p < length else (round(n * frac_at_end) if p == length else 1)
             for p in positions]
    rows = [[''] + [str(p) for p in positions], ['count'] + [f'{c:.1f}' for c in count]]
    for q, off in zip(QUANTILES, OFFSETS):
        rows.append([q] + [f'{max(2, min(41, median(p) + off)):.1f}' for p in positions])
    with open(path, 'w') as fh:
        fh.write('\n'.join('\t'.join(r) for r in rows) + '\n')


def fwd(p):          # 38 to 260, then down to 30 at 282, 20 from 283
    return 38 if p <= 260 else (38 - 8 * (p - 260) / 22 if p <= 282 else 20)


def rev(p, drop=250):   # 38 to 150, then down to 25 at `drop`, 22 after
    if p <= 150:
        return 38
    if p <= drop:
        return 38 - 13 * (p - 150) / (drop - 150)
    return 22


PROFILES = {
    '2x300':      [(283, 0.91, fwd), (281, 0.57, rev)],
    'early_drop': [(283, 0.91, fwd), (281, 0.57, lambda p: rev(p, 200))],
    '2x250':      [(233, 0.90, lambda p: 38 if p <= 232 else 20),
                   (230, 0.90, lambda p: 36 if p <= 229 else 20)],
    'bad_start':  [(283, 0.91, fwd), (281, 0.57, lambda p: 20)],
}

if __name__ == '__main__':
    out, profile = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for name, (length, frac, median) in zip(['forward', 'reverse'], PROFILES[profile]):
        write(os.path.join(out, f'{name}-seven-number-summaries.tsv'), length, frac, median)
