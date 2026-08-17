# shellcheck shell=bash
# Progress + failure localisation for the Gutbusters pipeline stages.  SOURCED, not executed.
#
# WHY THIS EXISTS. The entrypoint runs three stages, but 11 of the 15 tool invocations live inside
# 5_eval.sh alone. So "Step failed: 5_eval.sh" points at an 80 KB log and says nothing else -- the
# stage is the wrong unit for a pipeline whose cost and failure surface are both per-TOOL. Every
# tool announces itself through gb_step; the entrypoint reads the state file on failure and names
# the tool that was running.
#
# It also makes FAIL-OPEN tools visible. The four tools in 5_eval.sh's parallel block each end in
# `|| echo "WARN: ... failed"` and let the run continue with median-imputed columns -- by design
# (see the block's own comment), because one missing feature must not lose the whole assembly. But
# a WARN 600 lines into a stage log is not a report. gb_tool_result puts the outcome in the
# top-level output and records it, so a degraded run cannot look like a clean one.
#
# State lives in ONE file because each stage is a separate process: a shell variable would reset
# at every stage boundary and the counter would restart from zero three times.

_gb_state="${GB_TMP:-${GB_OUTDIR:-/tmp}}/.gb_progress"
_gb_degraded="${_gb_state}.degraded"

# Total tool count for THIS run, computed from the same gates the stages branch on.
#
# ⚠️ THIS IS THE DENOMINATOR'S ONLY SOURCE OF TRUTH, and it is deliberately NOT derived by
# scanning the stages -- a plan that reads the future cannot exist. The drift it invites (someone
# adds a gb_step and forgets this function) is therefore made LOUD instead of silent: gb_step
# appends "!PLAN" once the counter passes the total, and gb_progress_summary reports a run whose
# final count did not land on the plan. A wrong percentage that announces itself is recoverable;
# a wrong percentage that reads plausibly is the bug this whole file exists to prevent.
gb_plan_total() {
	local n=0
	# ---- 1_prophage.sh ----
	n=$((n + 1))                                          # seqkit length filter
	if [ "${GB_KRAKEN:-1}" = "1" ]; then n=$((n + 2)); fi  # host pass + host parse
	if [ "${GB_BLAST:-0}" = "1" ]; then n=$((n + 1)); fi   # UHGV claimer (retired; gate kept)
	n=$((n + 1))                                          # CheckV end_to_end (excision)
	if [ "${GB_PHAGEBOOST:-0}" = "1" ]; then n=$((n + 1)); fi   # retired 2026-08-10
	if [ "${GB_GENOMAD_CLAIM:-0}" = "1" ]; then n=$((n + 1)); fi # claimer (OFF; stage 5 still runs geNomad)
	# ---- 4_clean.sh runs no software: file assembly only, hence no step. ----
	# ---- 5_eval.sh ----
	n=$((n + 1))   # DeepVirFinder
	n=$((n + 1))   # MetaPhaPred
	# The geNomad + DVF-embed + composition + km4 block counts as ONE step, not four: they run
	# CONCURRENTLY and finish out of order, so numbering them 7,8,9,10 would make the percentage
	# advance in an order unrelated to wall clock and imply a sequence that does not exist. Their
	# individual outcomes are reported by gb_tool_result instead, which is the honest shape.
	n=$((n + 1))   # parallel feature block
	n=$((n + 1))   # CheckV annotation
	n=$((n + 1))   # MMseqs2 viral-vs-cellular contrast
	n=$((n + 1))   # pyhmmer / PHROGs
	if [ "${GB_KRAKEN:-1}" = "1" ]; then n=$((n + 1)); fi  # tagging + host merge
	n=$((n + 1))   # admission scoring (score_stack4.py)
	n=$((n + 1))   # final FASTA + table extraction
	echo "$n"
}

# gb_progress_init -- called ONCE by the entrypoint, before the first stage.
gb_progress_init() {
	local total
	total="$(gb_plan_total)"
	printf '0\t%s\tstarting\n' "$total" > "$_gb_state" 2>/dev/null || true
	rm -f "$_gb_degraded" 2>/dev/null || true
}

# gb_step "<label>" -- announce the tool ABOUT TO RUN. Announcing on entry rather than on
# completion is the whole point: a tool that dies must already have been named.
gb_step() {
	local label="${1:-?}" k total pct warn
	k=0; total=0
	if [ -s "$_gb_state" ]; then read -r k total _ < "$_gb_state" || true; fi
	case "$k" in ''|*[!0-9]*) k=0 ;; esac
	case "$total" in ''|*[!0-9]*) total=0 ;; esac
	if [ "$total" -le 0 ]; then total="$(gb_plan_total)"; fi
	k=$((k + 1))
	printf '%s\t%s\t%s\n' "$k" "$total" "$label" > "$_gb_state" 2>/dev/null || true
	pct=0
	if [ "$total" -gt 0 ]; then pct=$(( k * 100 / total )); fi
	warn=""
	if [ "$k" -gt "$total" ]; then warn="  !PLAN (step $k exceeds the planned $total -- gb_plan_total is stale)"; fi
	printf '[Gutbusters] [%2d/%2d %3d%%] %s%s\n' "$k" "$total" "$pct" "$label" "$warn"
}

# gb_tool_result <name> <rc> [seconds] [what degrades] -- outcome of a FAIL-OPEN tool.
# Never call this for a tool whose failure aborts the stage; those are reported by the entrypoint.
gb_tool_result() {
	local name="${1:-?}" rc="${2:-1}" secs="${3:--}" note="${4:-}"
	if [ "$rc" = "0" ]; then
		printf '[Gutbusters]            + %-22s OK      %6ss\n' "$name" "$secs"
	else
		printf '[Gutbusters]            + %-22s FAILED  rc=%s -- %s\n' "$name" "$rc" "$note"
		printf '%s\t%s\n' "$name" "$note" >> "$_gb_degraded" 2>/dev/null || true
	fi
}

# gb_progress_where -- what was running when things stopped. Read by the entrypoint on failure.
gb_progress_where() {
	local k total label
	if [ -s "$_gb_state" ]; then
		IFS=$'\t' read -r k total label < "$_gb_state" || true
		printf '%s (step %s/%s' "${label:-unknown}" "${k:-?}" "${total:-?}"
		if [ -n "${total:-}" ] && [ "${total:-0}" -gt 0 ] 2>/dev/null; then
			printf ', %d%%' $(( ${k:-0} * 100 / total ))
		fi
		printf ')\n'
	else
		printf 'unknown (no progress state at %s)\n' "$_gb_state"
	fi
}

# gb_progress_summary -- end-of-run report. Prints ONLY when there is something to say, so a
# clean run stays quiet, and a degraded one cannot be mistaken for a clean one.
gb_progress_summary() {
	local k total _l
	if [ -s "$_gb_degraded" ]; then
		echo "[Gutbusters] DEGRADED: the run completed, but these tools failed and their columns"
		echo "[Gutbusters]           were filled from the model's frozen training medians:"
		while IFS=$'\t' read -r _n _note; do
			printf '[Gutbusters]           - %-22s %s\n' "$_n" "$_note"
		done < "$_gb_degraded"
	fi
	if [ -s "$_gb_state" ]; then
		IFS=$'\t' read -r k total _l < "$_gb_state" || true
		case "$k$total" in ''|*[!0-9]*) return 0 ;; esac
		# FEWER than planned is normal and benign: stages skip a tool when its input is empty
		# (CheckV on an empty FASTA, the Kraken2 name resolver with nothing classified). The plan
		# is computed from CONFIG gates and cannot see data-dependent skips, so this is reported
		# as information, not as a fault.
		if [ "$k" -lt "$total" ]; then
			echo "[Gutbusters] $((total - k)) of $total tool steps were skipped (nothing for them to process)."
		fi
		# MORE than planned is a real defect: a gb_step exists that gb_plan_total does not know
		# about, so every percentage this run printed was against a wrong denominator.
		if [ "$k" -gt "$total" ]; then
			echo "[Gutbusters] BUG: ran $k tool steps but planned $total -- gb_plan_total in" >&2
			echo "[Gutbusters]      scripts/gb_progress.sh is out of date with the stages, and" >&2
			echo "[Gutbusters]      every percentage above was computed against the wrong total." >&2
		fi
	fi
}
