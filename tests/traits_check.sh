#!/usr/bin/env bash
# Run R/pr2_traits.R on a synthetic PR2 database with known answers
# (tests/make_synthetic_pr2_traits.R). Needs Rscript (base R; ggplot2 optional).
# Usage: bash tests/traits_check.sh WORKDIR     (WORKDIR must not exist yet)
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/traits_check.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
mkdir -p "$WORK"
status=0
Rscript "${REPO}/tests/make_synthetic_pr2_traits.R" "$WORK/a"
Rscript "${REPO}/R/pr2_traits.R" "$WORK/a/export" "$WORK/a/out" "$WORK/a/pr2.rds" > "$WORK/a.log" 2>&1 \
    || { echo "FAIL  pr2_traits.R failed: $(tail -3 "$WORK/a.log")"; exit 1; }
T="$WORK/a/out/asv_traits.tsv"

# field ASV COLUMN -> value of COLUMN for that ASV in asv_traits.tsv
field() {
    awk -F'\t' -v id="$1" -v col="$2" \
        'NR == 1 { for (i = 1; i <= NF; i++) if ($i == col) c = i; next } $1 == id { print $c }' "$T"
}
check() {  # check ASV COLUMN EXPECTED DESCRIPTION   (EXPECTED "~text" = contains text)
    local got; got="$(field "$1" "$2")"
    if [[ "$3" == "~"* && "$got" == *"${3#\~}"* ]] || [[ "$got" == "$3" ]]; then
        echo "PASS  $1 $4: $2 = $got"
    else
        echo "FAIL  $1 $4: $2 = '$got', want '$3'"; status=1
    fi
}
check asv01 photo_category phototroph                "genus, all species phototrophs"
check asv01 trait_rank Genus                         "trait taken at the assigned rank"
check asv02 photo_category non_photosynthetic        "override: Paraphysomonas is colourless"
check asv02 overrides_applied "~genus Paraphysomonas" "override recorded"
check asv03 ecological_function mixed                "class-level ASV, species disagree (4/6)"
check asv03 photo_category mixed_or_unknown          "mixed -> mixed_or_unknown"
check asv04 photo_category plant                     "Embryophyceae reported as plant"
check asv05 photo_category phototroph                "dinoflagellate genus, all species CM"
check asv05 HAB "taxon contains 2 HAB species"       "HAB at genus level"
check asv06 photo_category non_photosynthetic       "Protoperidinium: Adl Table 2.1 codes C only"
check asv06 photo_category_source "Adl et al. 2019 Table 2.1" "decision attributed to Adl"
check asv15 photo_category phototroph                "Gymnodinium: Adl codes L M"
check asv15 photosynthetic_capable yes               "in the manuscript-1 set"
check asv16 photo_category dinoflagellate_unresolved "Amphidinium: species-level variation (Adl legend)"
check asv16 adl2019_photosynthesis variable          "variable genus"
check asv16 photosynthetic_capable doubtful          "only in the sensitivity set"
check asv17 adl2019_photosynthesis not_in_table      "genus absent from Adl Table 2.1"
check asv17 photo_category dinoflagellate_unresolved "unlisted dinoflagellate stays unresolved"
check asv18 adl2019_codes O                          "ciliate genus row (Strombidium, Table 2.2)"
check asv18 adl2019_diet omnivore                    "diet words"
check asv19 adl2019_codes B                          "ciliate not listed: Table 2.2 default B"
check asv19 adl2019_match "~Ciliophora default"      "default recorded"
check asv07 photosynthetic_capable yes               "acquired phototrophy counts as capable"
check asv04 photosynthetic_capable "no (plant)"      "plants excluded"
check asv07 photo_category acquired_phototrophy      "pSNCM (Mesodinium rubrum)"
check asv08 habitat_flag "~no marine WoRMS record among 3 species" "freshwater-only genus flagged"
check asv09 habitat_flag ""                          "1 species with habitat data: not flagged"
check asv10 photo_category unassigned                "unclassified ASV"
check asv10 photosynthetic_capable unassigned        "unclassified ASV is not counted as capable"
check asv11 trait_rank Genus                         "species not in PR2: genus used"
check asv12 trait_rank Supergroup                    "division not in PR2: goes up to supergroup"
check asv12 photo_category mixed_or_unknown          "supergroup is mixed"
check asv13 HAB "HAB species"                        "species-level HAB"
check asv13 worms_id 12345                           "WoRMS id at species level"
check asv14 ev_culture_share 0.667                   "2 of 3 reference sequences from cultures"
check asv14 ev_habitat_keywords "~sediment (1)"      "keyword class from isolation source"
check asv14 ev_hosts "Caretta caretta (1)"           "host evidence"
check asv14 silva_consensus "Eukaryota;SAR;Diatomea;Navicula" "majority-rule SILVA consensus stops where names disagree"
check asv14 silva_consensus_share 1.00                "all 3 sequences agree to genus"
check asv01 worms_marine "2 of 2"                    "WoRMS habitat as counts"

# manuscript-1 exports: capable ASVs; sensitivity adds the doubtful (unresolved) dinoflagellates
exp_ids() { awk -F'\t' 'NR > 2 { printf "%s ", $1 }' "$WORK/a/out/$1/feature-table.tsv" | sed 's/ $//'; }
want="asv01 asv05 asv07 asv08 asv09 asv11 asv13 asv14 asv15"
got="$(exp_ids export_photosynthetic)"
if [[ "$got" == "$want" ]]; then echo "PASS  export_photosynthetic has the 9 capable ASVs"
else echo "FAIL  export_photosynthetic has '$got', want '$want'"; status=1; fi
want="asv01 asv05 asv07 asv08 asv09 asv11 asv13 asv14 asv15 asv16 asv17"
got="$(exp_ids export_photosynthetic_sensitivity)"
if [[ "$got" == "$want" ]]; then echo "PASS  sensitivity export adds the 2 unresolved dinoflagellates"
else echo "FAIL  sensitivity export has '$got', want '$want'"; status=1; fi
# collapsed table keeps every read
tot_ft="$(awk -F'\t' 'NR > 2 { for (i = 2; i <= NF; i++) s += $i } END { print s }' "$WORK/a/export/feature-table.tsv")"
tot_cat="$(awk -F'\t' 'NR > 1 { for (i = 2; i <= NF; i++) s += $i } END { print s }' "$WORK/a/out/table_by_photo_category.tsv")"
if [[ "$tot_ft" == "$tot_cat" ]]; then echo "PASS  table_by_photo_category keeps all $tot_ft reads"
else echo "FAIL  reads $tot_ft in, $tot_cat in categories"; status=1; fi

# version mismatch between classifier and trait data is reported
Rscript "${REPO}/tests/make_synthetic_pr2_traits.R" "$WORK/b" 5.0.0
Rscript "${REPO}/R/pr2_traits.R" "$WORK/b/export" "$WORK/b/out" "$WORK/b/pr2.rds" > "$WORK/b.log" 2>&1
if grep -q "classifier was built from PR2 5.0.0 but the trait data are PR2 5.1.1" "$WORK/b.log"; then
    echo "PASS  PR2 version mismatch reported"
else echo "FAIL  version mismatch not reported"; status=1; fi

exit "$status"
