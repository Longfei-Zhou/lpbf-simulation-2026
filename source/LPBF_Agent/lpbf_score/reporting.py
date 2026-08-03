"""Plain-language Markdown and self-contained HTML reports."""

from __future__ import annotations

import csv
import html
from pathlib import Path
from typing import Any


COMPONENTS: dict[str, dict[str, Any]] = {
    "coverage": {
        "label": "Top-surface coverage",
        "metric": "top_surface_cumulative_coverage_fraction",
        "format": "percent",
        "rule": ">=99% release; 95–99% review; <95% fail",
    },
    "fusion": {
        "label": "Interface fusion",
        "metric": "interface_cumulative_fusion_fraction",
        "format": "percent",
        "rule": ">=99% release; 95–99% review; <95% fail",
    },
    "fusion_margin": {
        "label": "LOF geometry margin",
        "metric": "lof_index_midpoint",
        "format": "number",
        "rule": "<=1 overlap; >1 LOF risk",
    },
    "keyhole_margin": {
        "label": "Keyhole proxy",
        "metric": "keyhole_aspect_ratio_p90_over_p10",
        "format": "number",
        "rule": "<0.5 lower; 0.5–0.8 review; >0.8 high risk",
    },
    "pool_consistency": {
        "label": "Pool stability",
        "metric": "melt_pool_dimension_mean_cv",
        "format": "number",
        "rule": "Lower CV is more stable",
    },
    "thermal_uniformity": {
        "label": "Solidification-field consistency",
        "metric": "thermal_log10_mad_mean",
        "format": "number",
        "rule": "Lower robust dispersion is more consistent",
    },
    "remelt": {
        "label": "Repeated remelting",
        "metric": "excess_remelt_mean_above_2",
        "format": "number",
        "rule": "No penalty for numMelt<=2; fewer repeats are preferred",
    },
}


def _number(value: Any, digits: int = 2) -> str:
    if value is None:
        return "—"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def _metric_text(value: Any, kind: str) -> str:
    if value is None:
        return "No data"
    if kind == "percent":
        return f"{100.0 * float(value):.2f}%"
    return _number(value, 3)


def _diagnostic_number(value: Any) -> str:
    if value is None:
        return "—"
    number = float(value)
    if number != 0.0 and (abs(number) >= 1e4 or abs(number) < 1e-3):
        return f"{number:.3e}"
    return f"{number:.3f}"


def _scaled_text(value: Any, scale: float, unit: str, digits: int = 1) -> str:
    if value is None:
        return "—"
    return f"{float(value) * scale:.{digits}f} {unit}"


def _decision_text(decision: str) -> str:
    return {
        "PASS": "Passed model screening",
        "FAIL": "Failed model screening",
        "REVIEW": "Engineering review required",
    }.get(decision, decision)


def _primary_problem(data: dict[str, Any]) -> dict[str, Any] | None:
    problems = data.get("problem_diagnosis", [])
    for severity in ("CRITICAL", "WARNING", "ADVISORY"):
        found = next((item for item in problems if item.get("severity") == severity), None)
        if found:
            return found
    return problems[0] if problems else None


def render_summary(data: dict[str, Any]) -> str:
    """Render a short report meant to be read before the technical appendix."""
    decision = str(data.get("decision", "REVIEW"))
    metrics = data.get("metrics", {})
    calibration = data.get("model_calibration", {})
    problem = _primary_problem(data)
    lines = [
        f"# Start Here: {data.get('layer_id', 'Layer')} Assessment",
        "",
        f"## {_decision_text(decision)} ({decision})",
        "",
        (
            f"**Primary reason: {problem.get('where', 'unknown location')} — "
            f"{problem.get('finding', 'inspect the technical report')}.**"
            if problem
            else "**No explicit release blocker was found at the current thresholds.**"
        ),
        "",
        "| Primary measure | Result | Interpretation |",
        "|---|---:|---|",
        f"| Continuous process score | **{_number(data.get('quality_score'), 1)}/100** | Margin across available continuous metrics; cannot override a hard failure |",
        f"| Top-surface fusion | **{_metric_text(metrics.get('top_surface_cumulative_coverage_fraction'), 'percent')}** | Checks scan completeness |",
        f"| Interface fusion | **{_metric_text(metrics.get('interface_cumulative_fusion_fraction'), 'percent')}** | Checks thermal connection to the previous layer |",
        f"| LOF index | **{_metric_text(metrics.get('lof_index_midpoint'), 'number')}** | <=1 is safer; >1 indicates LOF risk |",
        f"| Data completeness | **{_number(data.get('evidence_completeness_score'), 1)}/100** | Required files and fields present |",
        f"| Evidence adequacy | **{_number(data.get('evidence_adequacy_score'), 1)}/100** | Grid, sample range, and provenance support the conclusion |",
        f"| Model calibration | **{calibration.get('status', 'UNKNOWN')}** | {calibration.get('meaning', 'Not provided')} |",
        "",
    ]
    if decision == "FAIL" and float(data.get("quality_score", 0.0)) >= 75.0:
        lines.extend([
            "> A high continuous score can still fail because interface coverage and LOF are hard gates that cannot be averaged away.",
            "",
        ])

    lines.extend(["## Problems and next actions", ""])
    problems = data.get("problem_diagnosis", [])
    if problems:
        for index, item in enumerate(problems[:5], 1):
            actions = item.get("recommended_actions", [])
            action = actions[0] if actions else "Inspect the evidence in Assessment_Report.md."
            lines.append(
                f"{index}. **[{item.get('severity', 'INFO')}] "
                f"{item.get('where', 'whole layer')}:** {item.get('finding', '')}"
            )
            lines.append(f"   - Next: {action}")
    else:
        lines.append("1. No issue exceeds the reporting threshold; confirm calibration and model scope.")

    lines.extend([
        "",
        "## Seven component scores",
        "",
        "| Component | Score | Measured value | Rule |",
        "|---|---:|---:|---|",
    ])
    scores = data.get("component_scores", {})
    for key, spec in COMPONENTS.items():
        score = scores.get(key)
        value = metrics.get(spec["metric"])
        status = "Not scored" if score is None else f"{float(score):.1f}"
        lines.append(
            f"| {spec['label']} | {status} | {_metric_text(value, spec['format'])} | {spec['rule']} |"
        )

    diagnostics = data.get("process_diagnostics", {}).get("solidification", {})
    lines.extend(["", "## Solidification diagnostics (not release criteria)", ""])
    diagnostic_rows = (
        ("Cooling rate abs(dT/dt)", "cooling_rate_abs_k_s", "K/s"),
        ("Temperature gradient G", "temperature_gradient_k_m", "K/m"),
        ("Solidification velocity V", "solidification_velocity_m_s", "m/s"),
        ("G/V", "g_over_v_k_s_m2", "K·s/m²"),
        ("G×V", "g_times_v_k_s", "K/s"),
        ("Solidification event time tSol", "solidification_event_time_s", "s"),
    )
    lines.extend(["| Diagnostic | P10 | Median | P90 |", "|---|---:|---:|---:|"])
    for label, key, unit in diagnostic_rows:
        summary = diagnostics.get(key)
        if summary:
            lines.append(
                f"| {label} | {_diagnostic_number(summary.get('p10'))} {unit} | "
                f"{_diagnostic_number(summary.get('median'))} {unit} | "
                f"{_diagnostic_number(summary.get('p90'))} {unit} |"
            )
    lines.extend([
        "",
        "> These quantities have no universal larger-is-better direction. Use them for layer comparison, anomaly localization, and future material calibration.",
        "",
        "## File guide",
        "",
        "- `Dashboard.html`: most direct browser dashboard.",
        "- `00_READ_ME_FIRST.md`: this summary.",
        "- `Action_Plan.csv`: problems, evidence, and recommended actions.",
        "- `Assessment_Report.md`: formulas, sources, and full technical detail.",
        "- `assessment.json`: complete machine-readable result.",
        "",
        "> PASS/FAIL is a 3DThesis conduction-model screening result, not CT, metallography, or mechanical certification.",
    ])
    return "\n".join(lines) + "\n"


def _score_bar(label: str, score: Any) -> str:
    if score is None:
        return f'<div class="metric"><div><b>{html.escape(label)}</b><span>Not scored</span></div><div class="bar unavailable"></div></div>'
    value = max(0.0, min(100.0, float(score)))
    tone = "good" if value >= 80 else "warn" if value >= 60 else "bad"
    return (
        f'<div class="metric"><div><b>{html.escape(label)}</b><span>{value:.1f}</span></div>'
        f'<div class="bar"><i class="{tone}" style="width:{value:.1f}%"></i></div></div>'
    )


def render_dashboard(data: dict[str, Any]) -> str:
    """Render a dependency-free visual dashboard."""
    decision = str(data.get("decision", "REVIEW"))
    tone = {"PASS": "pass", "FAIL": "fail", "REVIEW": "review"}.get(decision, "review")
    problem = _primary_problem(data)
    metrics = data.get("metrics", {})
    scores = data.get("component_scores", {})
    calibration = data.get("model_calibration", {})
    component_html = "".join(
        _score_bar(spec["label"], scores.get(key)) for key, spec in COMPONENTS.items()
    )
    issue_rows = []
    for item in data.get("problem_diagnosis", [])[:8]:
        actions = item.get("recommended_actions", [])
        issue_rows.append(
            "<tr>"
            f"<td><span class='pill'>{html.escape(str(item.get('severity', 'INFO')))}</span></td>"
            f"<td><b>{html.escape(str(item.get('where', 'whole layer')))}</b><br>"
            f"{html.escape(str(item.get('finding', '')))}</td>"
            f"<td>{html.escape(str(item.get('evidence', '')))}</td>"
            f"<td>{html.escape('; '.join(str(value) for value in actions))}</td>"
            "</tr>"
        )
    if not issue_rows:
        issue_rows.append("<tr><td>INFO</td><td colspan='3'>No explicit anomaly at current thresholds.</td></tr>")

    where = problem.get("where", "current layer") if problem else "current layer"
    finding = problem.get("finding", "No explicit release blocker") if problem else "No explicit release blocker"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>LPBF Assessment Dashboard</title>
<style>
:root{{--ink:#172033;--muted:#667085;--line:#e5e7eb;--panel:#fff;--bg:#f5f7fb;--pass:#15803d;--fail:#c62828;--review:#b26a00}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}}
.wrap{{max-width:1120px;margin:auto;padding:28px}} h1,h2{{margin:0 0 14px}} h1{{font-size:26px}} h2{{font-size:18px}}
.hero,.panel{{background:var(--panel);border:1px solid var(--line);border-radius:16px;padding:22px;margin-bottom:18px;box-shadow:0 5px 18px #1822300d}}
.decision{{display:inline-block;padding:7px 13px;border-radius:999px;color:white;font-weight:800;letter-spacing:.04em}} .decision.pass{{background:var(--pass)}} .decision.fail{{background:var(--fail)}} .decision.review{{background:var(--review)}}
.lead{{font-size:18px;font-weight:700;margin:14px 0 4px}} .muted{{color:var(--muted)}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-top:18px}} .card{{background:#f8fafc;border:1px solid var(--line);border-radius:12px;padding:14px}} .card strong{{display:block;font-size:26px;margin-top:5px}} .card small{{color:var(--muted)}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}} .metric{{margin:13px 0}} .metric>div:first-child{{display:flex;justify-content:space-between}} .bar{{height:10px;background:#edf0f4;border-radius:20px;overflow:hidden;margin-top:5px}} .bar i{{display:block;height:100%}} .bar .good{{background:#22a06b}} .bar .warn{{background:#e4a11b}} .bar .bad{{background:#d64545}} .unavailable{{background:repeating-linear-gradient(45deg,#eee,#eee 5px,#ddd 5px,#ddd 10px)}}
table{{width:100%;border-collapse:collapse}} th,td{{text-align:left;vertical-align:top;border-bottom:1px solid var(--line);padding:10px 8px}} th{{color:var(--muted);font-size:13px}} .pill{{font-size:12px;background:#fee2e2;color:#991b1b;padding:3px 7px;border-radius:999px;font-weight:700}}
.warning{{border-left:5px solid #e4a11b;background:#fffaf0;padding:12px 14px;border-radius:8px}} code{{background:#eef2f6;padding:2px 5px;border-radius:4px}} @media(max-width:800px){{.cards,.grid{{grid-template-columns:1fr 1fr}}}} @media(max-width:520px){{.cards,.grid{{grid-template-columns:1fr}}.wrap{{padding:14px}}}}
</style></head><body><main class="wrap">
<section class="hero"><span class="decision {tone}">{html.escape(decision)} · {html.escape(_decision_text(decision))}</span>
<h1>{html.escape(str(data.get('layer_id','Layer')))}</h1>
<div class="lead">{html.escape(str(where))}: {html.escape(str(finding))}</div>
<div class="muted">A high continuous score cannot override a hard failure. Check interface fusion and LOF first.</div>
<div class="cards">
<div class="card"><small>Continuous process score</small><strong>{_number(data.get('quality_score'),1)}</strong><small>/100 across available metrics</small></div>
<div class="card"><small>Top / interface</small><strong>{_metric_text(metrics.get('top_surface_cumulative_coverage_fraction'),'percent')} / {_metric_text(metrics.get('interface_cumulative_fusion_fraction'),'percent')}</strong><small>authoritative numMelt coverage</small></div>
<div class="card"><small>Completeness / adequacy</small><strong>{_number(data.get('evidence_completeness_score'),0)} / {_number(data.get('evidence_adequacy_score'),0)}</strong><small>not quality scores</small></div>
<div class="card"><small>Model calibration</small><strong style="font-size:18px">{html.escape(str(calibration.get('status','UNKNOWN')))}</strong><small>{html.escape(str(calibration.get('meaning','Not provided')))}</small></div>
</div></section>
<div class="grid"><section class="panel"><h2>Seven component scores</h2>{component_html}</section>
<section class="panel"><h2>Primary physical measures</h2><table>
<tr><th>Measure</th><th>Value</th><th>Interpretation</th></tr>
<tr><td>LOF index</td><td><b>{_metric_text(metrics.get('lof_index_midpoint'),'number')}</b></td><td><=1 overlap; >1 risk</td></tr>
<tr><td>Melt-pool width</td><td>{_scaled_text(metrics.get('melt_pool_width_median_m'),1e6,'µm')}</td><td>Compare with hatch spacing</td></tr>
<tr><td>Melt-pool depth</td><td>{_scaled_text(metrics.get('melt_pool_depth_median_m'),1e6,'µm')}</td><td>Compare with layer thickness and interface coverage</td></tr>
<tr><td>Snapshot path span</td><td>{_metric_text(metrics.get('snapshot_scan_fraction_span'),'percent')}</td><td>Larger span is more representative</td></tr>
<tr><td>Peak temperature</td><td>{_number(metrics.get('snapshot_peak_temperature_max_k_diagnostic_only'),0)} K</td><td>Keyhole score is withheld outside model scope</td></tr>
</table><div class="warning"><b>Model boundary:</b> a conduction-only model cannot directly establish porosity, spatter, balling, cracks, or a real keyhole. PASS is not physical certification.</div></section></div>
<section class="panel"><h2>Problems and actions</h2><table><tr><th>Severity</th><th>Location and finding</th><th>Evidence</th><th>Next action</th></tr>{''.join(issue_rows)}</table></section>
<section class="panel"><h2>Report navigation</h2><p><code>00_READ_ME_FIRST.md</code> summary · <code>Action_Plan.csv</code> actions · <code>Assessment_Report.md</code> technical appendix · <code>assessment.json</code> machine data</p></section>
</main></body></html>"""


def write_intuitive_outputs(output_dir: Path, data: dict[str, Any]) -> None:
    """Write the files intended for humans before the technical appendix."""
    (output_dir / "00_READ_ME_FIRST.md").write_text(
        render_summary(data), encoding="utf-8", newline="\n"
    )
    (output_dir / "Dashboard.html").write_text(
        render_dashboard(data), encoding="utf-8", newline="\n"
    )
    with (output_dir / "Action_Plan.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["priority", "severity", "where", "finding", "evidence", "action", "confidence"],
        )
        writer.writeheader()
        for priority, item in enumerate(data.get("problem_diagnosis", []), 1):
            actions = item.get("recommended_actions", []) or [""]
            writer.writerow({
                "priority": priority,
                "severity": item.get("severity", "INFO"),
                "where": item.get("where", ""),
                "finding": item.get("finding", ""),
                "evidence": item.get("evidence", ""),
                "action": "; ".join(str(value) for value in actions),
                "confidence": item.get("confidence", ""),
            })
