import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from lpbf_score.scorer import (
    LayerAssessment,
    ProcessParameters,
    analyse_rdf,
    interpolated_pool_depth,
    parse_cli_layer_thickness,
    parse_path_file,
    powered_segments,
)


ROOT = Path(__file__).resolve().parent.parent


class PhysicsScorerTests(unittest.TestCase):
    def setUp(self):
        self.assessment = LayerAssessment(
            solidification_path=ROOT / "tests" / "fixtures",
            snapshots_path=None,
            output_dir=ROOT / "output" / "_test_unused",
            config_path=ROOT / "config" / "scoring.yaml",
            layer_id="test",
        )
        self.parameters = ProcessParameters(
            liquidus_k=1708.0,
            liquidus_source="test Material.txt",
            layer_thickness_m=25e-6,
            layer_thickness_source="test layer thickness",
            hatch_spacing_m=100e-6,
            hatch_spacing_source="test hatch spacing",
            grid_x_m=50e-6,
            grid_y_m=50e-6,
            grid_z_m=25e-6,
            power_w=200.0,
            efficiency=0.35,
            scan_speed_m_s=0.8,
        )
        count = 30
        gradient = np.linspace(8e6, 12e6, count)
        velocity = np.linspace(0.08, 0.12, count)
        self.solid = pd.DataFrame(
            {
                "x": np.arange(count) * 50e-6,
                "y": np.zeros(count),
                "z": np.zeros(count),
                "G": gradient,
                "V": velocity,
                "dTdt": gradient * velocity,
                "numMelt": np.where(np.arange(count) % 4 == 0, 2, 1),
            }
        )

    def evaluate(
        self,
        snapshots: pd.DataFrame,
        snapshot_count: int = 3,
        authoritative_coverage: float | None = 1.0,
        scan_fractions: list[float] | None = None,
    ):
        coverage = {
            "status": (
                "COMPLETE_EVIDENCE"
                if authoritative_coverage is not None
                else "INSUFFICIENT_EVIDENCE"
            ),
            "authoritative": authoritative_coverage is not None,
            "method": (
                "complete target-grid numMelt"
                if authoritative_coverage is not None
                else None
            ),
            "coverage_fraction": authoritative_coverage,
            "minimum_required_fraction": 0.99,
            "snapshot_cumulative_temperature": {
                "available": False,
                "scan_fraction_metadata": {
                    "values": scan_fractions or [],
                    "count": len(scan_fractions or []),
                },
            },
            "num_melt_target_grid": {
                "available": True,
                "tracking_mode": "None",
            },
            "provenance": {},
        }
        return self.assessment.score(
            self.solid,
            snapshots,
            self.parameters,
            solid_files=[Path("solid.csv")],
            snapshot_files=[Path(f"snapshot-{i}.csv") for i in range(snapshot_count)],
            summaries={},
            histories=None,
            coverage=coverage,
        )

    def test_good_geometry_passes(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        result = self.evaluate(snapshots)
        self.assertEqual(result["decision"], "PASS")
        self.assertTrue(result["predicted_success"])
        self.assertLess(result["metrics"]["lof_index_conservative"], 1.0)

    def test_lack_of_fusion_geometry_fails(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [115e-6, 120e-6, 125e-6],
                "depth_m": [38e-6, 40e-6, 42e-6],
                "length_m": [300e-6, 310e-6, 305e-6],
            }
        )
        result = self.evaluate(snapshots)
        self.assertEqual(result["decision"], "FAIL")
        self.assertFalse(result["predicted_success"])
        self.assertGreater(result["metrics"]["lof_index_conservative"], 1.0)

    def test_missing_geometry_requires_review(self):
        result = self.evaluate(pd.DataFrame(), snapshot_count=0)
        self.assertEqual(result["decision"], "REVIEW")
        self.assertIsNone(result["predicted_success"])
        self.assertIsNone(result["component_scores"]["fusion_margin"])

    def test_missing_cumulative_coverage_requires_review(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        result = self.evaluate(snapshots, authoritative_coverage=None)
        self.assertEqual(result["decision"], "REVIEW")
        self.assertIsNone(result["predicted_success"])
        self.assertIsNone(result["component_scores"]["coverage"])

    def test_incomplete_authoritative_coverage_fails(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        result = self.evaluate(snapshots, authoritative_coverage=0.90)
        self.assertEqual(result["decision"], "FAIL")
        self.assertFalse(result["predicted_success"])

    def test_near_release_coverage_requires_review_not_brittle_fail(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        result = self.evaluate(snapshots, authoritative_coverage=0.98)
        self.assertEqual(result["decision"], "REVIEW")
        self.assertFalse(
            any(flag.startswith("FAIL:") for flag in result["flags"])
        )

    def test_lof_grid_uncertainty_is_review_not_hard_fail(self):
        self.parameters.layer_thickness_m = 50e-6
        snapshots = pd.DataFrame(
            {
                "width_m": [150e-6, 151e-6, 160e-6, 170e-6],
                "depth_m": [75e-6, 75e-6, 75e-6, 75e-6],
                "length_m": [500e-6, 510e-6, 520e-6, 530e-6],
            }
        )
        result = self.evaluate(snapshots, snapshot_count=4)
        self.assertEqual(result["decision"], "REVIEW")
        self.assertTrue(result["metrics"]["lof_resolution_sensitive"])
        self.assertLess(result["metrics"]["lof_index_nominal"], 1.0)
        self.assertGreater(result["metrics"]["lof_index_conservative"], 1.0)

    def test_evidence_completeness_and_adequacy_are_separate(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        result = self.evaluate(snapshots, snapshot_count=3)
        self.assertGreater(
            result["evidence_completeness_score"],
            result["evidence_adequacy_score"],
        )
        self.assertIn("fusion_margin", result["component_scores"])

    def test_intuitive_diagnostic_outputs_are_written(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        result = self.evaluate(snapshots)
        with tempfile.TemporaryDirectory() as temp_dir:
            self.assessment.output_dir = Path(temp_dir)
            self.assessment.write_outputs(result, snapshots, summaries={})
            for filename in (
                "00_READ_ME_FIRST.md",
                "Dashboard.html",
                "Action_Plan.csv",
                "Assessment_Report.md",
                "Problem_Diagnosis.csv",
                "Coverage_Diagnostics.csv",
                "Score_Breakdown.csv",
                "assessment.json",
            ):
                self.assertTrue((Path(temp_dir) / filename).is_file())

    def test_pysimle_diagnostics_are_integrated_without_extra_score(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        result = self.evaluate(snapshots)
        diagnostics = result["process_diagnostics"]
        self.assertEqual(diagnostics["status"], "DIAGNOSTIC_ONLY_NOT_SCORED")
        self.assertGreater(
            diagnostics["solidification"]["g_over_v_k_s_m2"]["median"],
            0.0,
        )
        self.assertIn("self-normalisation using the same evaluated dataset", diagnostics["legacy_methods_rejected"])

    def test_snapshot_path_span_reduces_representativeness(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6] * 20,
                "depth_m": [110e-6] * 20,
                "length_m": [800e-6] * 20,
            }
        )
        result = self.evaluate(
            snapshots,
            snapshot_count=20,
            scan_fractions=list(np.linspace(0.15, 0.36, 20)),
        )
        self.assertAlmostEqual(
            result["metrics"]["snapshot_scan_fraction_span"], 0.21
        )
        self.assertLess(
            result["evidence_adequacy_parts"]["geometry_path_representativeness"],
            4.0,
        )
        self.assertTrue(
            any(
                item.get("category") == "Evidence representativeness"
                for item in result["problem_diagnosis"]
            )
        )

    def test_keyhole_score_is_withheld_outside_thermal_model_scope(self):
        snapshots = pd.DataFrame(
            {
                "width_m": [150e-6, 155e-6, 160e-6],
                "depth_m": [60e-6, 65e-6, 70e-6],
                "length_m": [500e-6, 510e-6, 520e-6],
                "max_temperature_k_diagnostic_only": [6000.0, 6200.0, 6100.0],
            }
        )
        result = self.evaluate(snapshots)
        self.assertIsNone(result["component_scores"]["keyhole_margin"])
        self.assertEqual(
            result["metrics"]["keyhole_score_status"],
            "NOT_SCORED_MODEL_OUT_OF_SCOPE",
        )

    def test_cli_layer_thickness_is_not_domain_resolution(self):
        parsed = parse_cli_layer_thickness(ROOT / "tests" / "fixtures" / "sample.cli")
        self.assertIsNotNone(parsed)
        self.assertAlmostEqual(parsed["layer_thickness_m"], 50e-6)

    def test_concatenated_extreme_path_mirrors_first_record_per_line(self):
        frame = parse_path_file(
            ROOT / "tests" / "fixtures" / "concatenated_path.txt"
        )
        self.assertIsNotNone(frame)
        self.assertEqual(len(frame), 2)
        self.assertEqual(len(powered_segments(frame)), 1)
        diagnostics = frame.attrs["path_parse_diagnostics"]
        self.assertTrue(diagnostics["extreme_concatenated_line_format_detected"])
        self.assertEqual(diagnostics["trailing_token_line_count"], 1)

    def test_empty_rdf_is_json_safe(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            rdf_path = Path(temp_dir) / "empty.RDF.Final.csv"
            rdf_path.write_text("x,y,z,tm,tl,cr\n", encoding="utf-8")
            result = analyse_rdf([rdf_path])
            self.assertIsNotNone(result)
            self.assertEqual(result["row_count"], 0)
            self.assertIsNone(result["negative_liquid_duration_fraction"])

    def test_pool_cell_count_uses_snapshot_grid_not_domain_grid(self):
        """Pool cell counts must use the grid that produced the geometry."""
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
                "grid_x_m": [10e-6, 10e-6, 10e-6],
                "grid_z_m": [5e-6, 5e-6, 5e-6],
            }
        )
        metrics = self.evaluate(snapshots)["metrics"]
        self.assertAlmostEqual(metrics["geometry_grid_x_m"], 10e-6)
        self.assertAlmostEqual(metrics["geometry_grid_z_m"], 5e-6)
        # The snapshot grid gives min(29.6, 21.2), not the Domain-grid 4.24.
        cells = metrics["minimum_cells_across_p10_pool_dimension"]
        self.assertGreater(cells, 15.0)
        self.assertAlmostEqual(
            cells,
            min(
                float(np.quantile(snapshots["width_m"], 0.10)) / 10e-6,
                float(np.quantile(snapshots["depth_m"], 0.10)) / 5e-6,
            ),
        )

    def test_pool_cell_count_falls_back_to_domain_grid(self):
        """Fall back to the Domain grid when snapshots lack grid metadata."""
        snapshots = pd.DataFrame(
            {
                "width_m": [300e-6, 310e-6, 295e-6],
                "depth_m": [110e-6, 105e-6, 108e-6],
                "length_m": [800e-6, 820e-6, 810e-6],
            }
        )
        metrics = self.evaluate(snapshots)["metrics"]
        self.assertAlmostEqual(metrics["geometry_grid_x_m"], self.parameters.grid_x_m)
        self.assertAlmostEqual(metrics["geometry_grid_z_m"], self.parameters.grid_z_m)


class InterpolatedPoolDepthTests(unittest.TestCase):
    """Exercise sub-grid pool-depth interpolation."""

    LIQUIDUS = 1708.0
    DX = DY = 10e-6
    DZ = 5e-6

    @staticmethod
    def column(z_values, temperatures):
        return pd.DataFrame(
            {
                "x": [0.0] * len(z_values),
                "y": [0.0] * len(z_values),
                "z": list(z_values),
                "T": list(temperatures),
            }
        )

    def test_depth_is_not_quantised_to_the_z_grid(self):
        # The liquidus crossing between -10 and -15 um lies at -11.4 um.
        pool = self.column([0.0, -5e-6, -10e-6], [2200.0, 1900.0, 1750.0])
        shoulder = self.column([-15e-6], [1600.0])
        depth, bottom_limited, fraction = interpolated_pool_depth(
            pool, shoulder, self.LIQUIDUS, self.DX, self.DY, self.DZ,
            domain_z_min=-40e-6,
        )
        self.assertAlmostEqual(depth, 11.4e-6, places=9)
        self.assertFalse(bottom_limited)
        self.assertAlmostEqual(fraction, 1.0)
        self.assertNotAlmostEqual(depth % self.DZ, 0.0, places=9)

    def test_pool_touching_domain_floor_is_flagged_as_lower_bound(self):
        # A molten cell on the domain floor makes depth a lower bound.
        pool = self.column([0.0, -5e-6, -10e-6], [2200.0, 1900.0, 1750.0])
        shoulder = self.column([], [])
        depth, bottom_limited, _ = interpolated_pool_depth(
            pool, shoulder, self.LIQUIDUS, self.DX, self.DY, self.DZ,
            domain_z_min=-10e-6,
        )
        self.assertTrue(bottom_limited)
        self.assertAlmostEqual(depth, 10e-6, places=9)

    def test_missing_shoulder_falls_back_to_half_cell(self):
        # Without a shoulder cell, use a half-cell without claiming interpolation.
        pool = self.column([0.0, -5e-6, -10e-6], [2200.0, 1900.0, 1750.0])
        shoulder = self.column([], [])
        depth, bottom_limited, fraction = interpolated_pool_depth(
            pool, shoulder, self.LIQUIDUS, self.DX, self.DY, self.DZ,
            domain_z_min=-40e-6,
        )
        self.assertAlmostEqual(depth, 12.5e-6, places=9)
        self.assertFalse(bottom_limited)
        self.assertAlmostEqual(fraction, 0.0)

    def test_deepest_column_wins_across_the_pool(self):
        # Melt depth is set by the deepest XY column.
        shallow = pd.DataFrame(
            {"x": [0.0, 0.0], "y": [0.0, 0.0],
             "z": [0.0, -5e-6], "T": [2200.0, 1750.0]}
        )
        deep = pd.DataFrame(
            {"x": [10e-6] * 3, "y": [0.0] * 3,
             "z": [0.0, -5e-6, -10e-6], "T": [2300.0, 2000.0, 1750.0]}
        )
        pool = pd.concat([shallow, deep], ignore_index=True)
        shoulder = pd.DataFrame(
            {"x": [0.0, 10e-6], "y": [0.0, 0.0],
             "z": [-10e-6, -15e-6], "T": [1600.0, 1600.0]}
        )
        depth, _, fraction = interpolated_pool_depth(
            pool, shoulder, self.LIQUIDUS, self.DX, self.DY, self.DZ,
            domain_z_min=-40e-6,
        )
        self.assertAlmostEqual(depth, 11.4e-6, places=9)
        self.assertAlmostEqual(fraction, 1.0)

    def test_empty_pool_is_safe(self):
        empty = pd.DataFrame(columns=["x", "y", "z", "T"])
        self.assertEqual(
            interpolated_pool_depth(
                empty, empty, self.LIQUIDUS, self.DX, self.DY, self.DZ
            ),
            (0.0, False, 0.0),
        )


if __name__ == "__main__":
    unittest.main()
