import tempfile
import unittest
from pathlib import Path

import numpy as np

from lpbf_score.coverage_generator import generate_coverage_points
from lpbf_score.scorer import (
    cli_solid_region_mask,
    match_cli_solid_geometry,
    target_region_mask,
)


ROOT = Path(__file__).resolve().parent.parent


class CoverageGeneratorTests(unittest.TestCase):
    def _case(self, root: Path) -> Path:
        case = root / "case"
        case.mkdir()
        (case / "Domain.txt").write_text(
            "X\n{\n Min 0\n Max 0.0002\n Res 0.0001\n}\n"
            "Y\n{\n Min 0\n Max 0.0002\n Res 0.0001\n}\n"
            "Z\n{\n Min -0.0001\n Max 0\n Res 0.00005\n}\n",
            encoding="utf-8",
        )
        (case / "Path.txt").write_text(
            "Mode X(mm) Y(mm) Z(mm) Pmod Vel(m/s)/Time(s)\n"
            "1 0 0 0 0 1e-7\n"
            "0 0.2 0 0 1 0.8\n",
            encoding="utf-8",
        )
        return case

    def test_generates_top_and_interface_points_without_final_csv(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            case = self._case(Path(temp_dir))
            points, summary = generate_coverage_points(
                case_dir=case,
                source_cli=ROOT / "tests" / "fixtures" / "sample.cli",
                hatch_spacing_um=100,
            )
            self.assertEqual(len(points), 6)
            self.assertEqual(summary["powered_segment_count"], 1)
            self.assertEqual(set(np.round(points["z_mm"], 9)), {0.0, -0.05})
            self.assertEqual(
                points.groupby("z_mm").size().sort_index().tolist(),
                [3, 3],
            )

    def test_custom_domain_requires_regular_domain_override(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            case = self._case(Path(temp_dir))
            (case / "Domain.txt").write_text(
                "Custom\n{\n File coverage_target_points.txt\n}\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "--domain-file"):
                generate_coverage_points(
                    case_dir=case,
                    layer_thickness_um=50,
                    hatch_spacing_um=100,
                )

    def test_optimized_corridor_mask_matches_brute_force(self):
        rng = np.random.default_rng(7)
        xy = rng.uniform(-0.001, 0.001, size=(500, 2))
        segments = rng.uniform(-0.0008, 0.0008, size=(12, 4))
        radius = 75e-6
        expected = np.zeros(len(xy), dtype=bool)
        for index, point in enumerate(xy):
            for x0, y0, x1, y1 in segments:
                start = np.array([x0, y0])
                vector = np.array([x1 - x0, y1 - y0])
                fraction = np.clip(
                    ((point - start) @ vector) / (vector @ vector),
                    0.0,
                    1.0,
                )
                if np.linalg.norm(point - (start + fraction * vector)) <= radius:
                    expected[index] = True
                    break
        np.testing.assert_array_equal(
            target_region_mask(xy, segments, radius),
            expected,
        )

    def test_cli_hole_and_exterior_voids_are_excluded(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cli = Path(temp_dir) / "with_hole.cli"
            cli.write_text(
                "$$HEADERSTART\n$$ASCII\n$$UNITS/1\n$$HEADEREND\n"
                "$$GEOMETRYSTART\n$$LAYER/0\n"
                "$$POLYLINE/1,1,5,0,0,0.3,0,0.3,0.3,0,0.3,0,0\n"
                "$$POLYLINE/1,0,5,0.1,0.1,0.2,0.1,0.2,0.2,0.1,0.2,0.1,0.1\n"
                "$$HATCHES/1,2,0,0.05,0.3,0.05,0,0.25,0.3,0.25\n"
                "$$GEOMETRYEND\n",
                encoding="utf-8",
            )
            segments = np.asarray(
                [
                    [0.0, 0.00005, 0.0003, 0.00005],
                    [0.0, 0.00025, 0.0003, 0.00025],
                ]
            )
            geometry, metadata = match_cli_solid_geometry(cli, segments)
            self.assertTrue(metadata["applied"])
            points = np.asarray(
                [
                    [0.00005, 0.00015],
                    [0.00015, 0.00015],
                    [0.00040, 0.00015],
                ]
            )
            np.testing.assert_array_equal(
                cli_solid_region_mask(points, geometry),
                np.asarray([True, False, False]),
            )

    def test_complete_cli_solid_is_denominator_not_path_intersection(self):
        """A missing scan corridor must not erase valid solid from the target."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case = root / "case"
            case.mkdir()
            (case / "Domain.txt").write_text(
                "X\n{\n Min 0\n Max 0.0004\n Res 0.0001\n}\n"
                "Y\n{\n Min 0\n Max 0.0004\n Res 0.0001\n}\n"
                "Z\n{\n Min -0.0001\n Max 0\n Res 0.00005\n}\n",
                encoding="utf-8",
            )
            (case / "Path.txt").write_text(
                "Mode X(mm) Y(mm) Z(mm) Pmod Vel(m/s)/Time(s)\n"
                "1 0 0.05 0 0 1e-7\n"
                "0 0.3 0.05 0 1 0.8\n",
                encoding="utf-8",
            )
            cli = root / "solid.cli"
            cli.write_text(
                "$$HEADERSTART\n$$ASCII\n$$UNITS/1\n$$HEADEREND\n"
                "$$GEOMETRYSTART\n$$LAYER/0\n"
                "$$POLYLINE/1,1,5,0,0,0.3,0,0.3,0.3,0,0.3,0,0\n"
                "$$HATCHES/1,1,0,0.05,0.3,0.05\n$$GEOMETRYEND\n",
                encoding="utf-8",
            )
            points, summary = generate_coverage_points(
                case_dir=case,
                source_cli=cli,
                layer_thickness_um=50,
                hatch_spacing_um=100,
            )
            mask_info = summary["cli_solid_mask"]
            self.assertTrue(mask_info["applied"])
            self.assertGreater(mask_info["solid_points_outside_path_corridor"], 0)
            self.assertEqual(
                len(points) // 2,
                mask_info["solid_target_point_count"],
            )


if __name__ == "__main__":
    unittest.main()
