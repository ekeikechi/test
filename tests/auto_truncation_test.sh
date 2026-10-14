#!/usr/bin/env bash
# Unit test for pipeline/auto_truncation.py on synthetic quality profiles with known
# answers (tests/make_seven_number.py). Needs python3 only.
# Usage: bash tests/auto_truncation_test.sh WORKDIR     (WORKDIR must not exist yet)
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/auto_truncation_test.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
mkdir -p "$WORK"
status=0

# case NAME PROFILE "EXTRA ARGS" EXPECTED_STDOUT EXPECTED_EXIT [MESSAGE_FRAGMENT]
case_() {
    local name="$1" profile="$2" extra="$3" want="$4" want_exit="$5" fragment="${6:-}"
    local dir="$WORK/$profile" got rc=0
    [[ -d "$dir" ]] || python3 "${REPO}/tests/make_seven_number.py" "$dir" "$profile"
    # shellcheck disable=SC2086  # $extra holds several arguments
    got="$(python3 "${REPO}/pipeline/auto_truncation.py" \
        --forward "$dir/forward-seven-number-summaries.tsv" \
        --reverse "$dir/reverse-seven-number-summaries.tsv" \
        $extra --out "$WORK/$name.tsv" 2> "$WORK/$name.err")" || rc=$?
    if [[ "$got" == "$want" && "$rc" == "$want_exit" ]] \
       && { [[ -z "$fragment" ]] || grep -q "$fragment" "$WORK/$name.err"; }; then
        echo "PASS  $name -> ${got:-exit $rc}"
    else
        echo "FAIL  $name: got '$got' (exit $rc), want '$want' (exit $want_exit)"
        cat "$WORK/$name.err"; status=1
    fi
}

case_ quality_rule_only          2x300      ""                 "282 250" 0
case_ target_520_already_met     2x300      "--sum-target 520" "282 250" 0
case_ target_540_lengthens_R     2x300      "--sum-target 540" "282 258" 0
case_ early_drop_no_target       early_drop ""                 "282 200" 0
case_ early_drop_lengthened      early_drop "--sum-target 520" "282 238" 0
case_ 2x250_unreachable          2x250      "--sum-target 520" ""        2 "Overlap target not reachable: forward + reverse can reach at most 463 nt"
case_ reverse_bad_from_start     bad_start  "--sum-target 520" ""        2 "reverse reads: no truncation length qualifies"
case_ stricter_quality_q30       2x300      "--min-median-q 30" "282 211" 0

# the decision table records inputs, both cuts, the reason and the status
d="$WORK/early_drop_lengthened.tsv"
for kv in $'reverse_quality_cut\t200' $'reverse_length_cut\t280' $'reverse_chosen\t238' \
          $'reverse_reason\tlengthened for overlap' $'forward_reason\tquality cut' \
          $'sum\t520' $'max_mergeable_insert\t508' $'status\tok'; do
    grep -qx "$kv" "$d" || { echo "FAIL  decision table lacks '${kv/$'\t'/ = }'"; status=1; }
done
[[ $status == 0 ]] && echo "PASS  decision table complete"
if grep -qx $'status\toverlap target not reachable' "$WORK/2x250_unreachable.tsv"; then
    echo "PASS  unreachable target recorded in the decision table"
else
    echo "FAIL  unreachable status not recorded"; status=1
fi

exit "$status"
