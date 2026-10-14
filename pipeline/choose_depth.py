"""Choose rarefaction depths from per-sample read counts (stage 08, and the same rule in R).

  main depth        = the smallest sample left after MIN_SAMPLE_DEPTH (stage 06), so every
                      remaining sample is kept
  sensitivity depth = the deepest depth at which every group (combination of the
                      --groups metadata columns) keeps at least --min-keep of its samples
                      (default 0.8); omitted when no groups are given

Prints "MAIN SENSITIVITY" (SENSITIVITY empty without groups) and writes:
  --out         setting/value table (depths, rule, samples kept)
  --groups-out  one row per group: samples, kept at each depth, samples dropped at the
                sensitivity depth
"""
import argparse
import csv
import math
import sys


def read_depths(path):
    with open(path, newline='') as fh:
        rows = [r for r in csv.reader(fh, delimiter='\t') if r and not r[0].startswith('#')]
    return {r[0]: float(r[1]) for r in rows[1:]}


def read_metadata(path, cols):
    with open(path, newline='') as fh:
        rows = [r for r in csv.reader(fh, delimiter='\t') if r]
    header, body = rows[0], [r for r in rows[1:] if not r[0].startswith('#q2:')]
    missing = [c for c in cols if c not in header]
    if missing:
        sys.exit(f'{path}: no column(s) {", ".join(missing)} for the sensitivity groups')
    idx = [header.index(c) for c in cols]
    return {r[0]: ' / '.join(r[i] if i < len(r) else '' for i in idx) for r in body}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--depths', required=True, help='TSV: sample, reads (header line)')
    ap.add_argument('--metadata')
    ap.add_argument('--groups', default='', help='space-separated metadata columns')
    ap.add_argument('--min-keep', type=float, default=0.8)
    ap.add_argument('--out', required=True)
    ap.add_argument('--groups-out')
    a = ap.parse_args()

    depth = read_depths(a.depths)
    if not depth:
        sys.exit(f'{a.depths}: no samples')
    main_depth = int(min(depth.values()))
    rows = [('setting', 'value'), ('samples', len(depth)), ('main_depth', main_depth),
            ('main_rule', 'smallest remaining sample (all samples kept)')]
    sens = ''
    cols = a.groups.split()
    if cols:
        if not a.metadata:
            sys.exit('--groups needs --metadata')
        meta = read_metadata(a.metadata, cols)
        absent = [s for s in depth if s not in meta]
        if absent:
            sys.exit(f'samples missing from the metadata: {", ".join(sorted(absent)[:10])}')
        groups = {}
        for s, d in depth.items():
            groups.setdefault(meta[s], []).append(d)
        # per group: the depth at which ceil(min_keep * n) samples still have >= that depth
        allowed = {}
        for g, ds in groups.items():
            ds = sorted(ds, reverse=True)
            need = max(1, math.ceil(a.min_keep * len(ds)))
            allowed[g] = ds[need - 1]
        sens = int(min(allowed.values()))
        limiting = sorted(g for g, v in allowed.items() if int(v) == sens)
        rows += [('sensitivity_depth', sens),
                 ('sensitivity_rule', f'deepest depth keeping >= {a.min_keep:g} of samples in every '
                                      f'group of: {" x ".join(cols)}'),
                 ('limiting_groups', '; '.join(limiting)),
                 ('samples_kept_at_sensitivity', sum(d >= sens for d in depth.values()))]
        if a.groups_out:
            with open(a.groups_out, 'w', newline='') as fh:
                w = csv.writer(fh, delimiter='\t', lineterminator='\n')
                w.writerow(['group', 'samples', 'kept_at_main', 'kept_at_sensitivity',
                            'dropped_at_sensitivity'])
                for g in sorted(groups):
                    dropped = sorted(s for s, d in depth.items() if meta[s] == g and d < sens)
                    w.writerow([g, len(groups[g]), sum(d >= main_depth for d in groups[g]),
                                sum(d >= sens for d in groups[g]), ','.join(dropped)])
    with open(a.out, 'w', newline='') as fh:
        csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(rows)
    print(f'{main_depth} {sens}'.strip() if sens == '' else f'{main_depth} {sens}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
