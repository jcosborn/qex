"""Regression checks for the numeric comparator used by the HMC tests."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


DIFFNUM = Path(__file__).with_name("diffnum")
KINETIC_ATOL = 64 * sys.float_info.epsilon * 16 * 8**4


class DiffnumTests(unittest.TestCase):
    def compare(self, reference, actual, *options):
        with tempfile.TemporaryDirectory() as directory:
            reference_path = Path(directory) / "reference"
            actual_path = Path(directory) / "actual"
            reference_path.write_text(reference + "\n")
            actual_path.write_text(actual + "\n")
            return subprocess.run(
                [sys.executable, str(DIFFNUM), str(reference_path),
                 str(actual_path), *options],
                capture_output=True, text=True,
            )

    def kinetic_options(self, atol=KINETIC_ATOL):
        return ("1e-12", "--atol-after", "T:", str(atol))

    def test_observed_roundoff_needs_absolute_allowance(self):
        reference = "Beginning H: 73737.68198514577 T: -103.12184584984789 V: 73840.80383099562"
        actual = "Beginning H: 73737.68198514599 T: -103.12184584964416 V: 73840.80383099563"
        self.assertNotEqual(self.compare(reference, actual, "1e-12").returncode, 0)
        result = self.compare(reference, actual, *self.kinetic_options())
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_kinetic_allowance_scales_with_volume(self):
        reference, actual = "T: 0", "T: 2e-10"
        self.assertEqual(self.compare(reference, actual, *self.kinetic_options()).returncode, 0)
        small_atol = 64 * sys.float_info.epsilon * 16 * 4**4
        self.assertNotEqual(self.compare(reference, actual, *self.kinetic_options(small_atol)).returncode, 0)

    def test_absolute_and_relative_allowances_are_added(self):
        # Neither allowance alone suffices: 1.5e-9 < 1e-9 + 1e-12 * 1000.
        result = self.compare("T: 1000", "T: 1000.0000000015", *self.kinetic_options(1e-9))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_real_kinetic_error_fails(self):
        result = self.compare("T: -103.12184584984789", "T: -103.12184574984789", *self.kinetic_options())
        self.assertNotEqual(result.returncode, 0)

    def test_other_fields_keep_relative_tolerance(self):
        for label in ("H:", "V:", "even:", "odd:"):
            with self.subTest(label=label):
                result = self.compare(f"{label} 0.8", f"{label} 0.8000000002", *self.kinetic_options())
                self.assertNotEqual(result.returncode, 0)

    def test_legacy_relative_and_separator_arguments(self):
        self.assertEqual(self.compare("x: 1", "x: 1.0000000000001", "1e-12").returncode, 0)
        self.assertEqual(self.compare("x;1", "x;1.0000000000001", "1e-12", ";").returncode, 0)
        self.assertNotEqual(self.compare("x: 1", "x: 1.0001", "1e-12").returncode, 0)

    def test_malformed_or_nonfinite_output_fails(self):
        for actual in ("T:", "V: 0", "T: nan", "T: inf", "T: 0\nT: 0"):
            with self.subTest(actual=actual):
                self.assertNotEqual(self.compare("T: 0", actual, *self.kinetic_options()).returncode, 0)
        for value in ("nan", "inf", "-inf"):
            self.assertNotEqual(self.compare(f"T: {value}", f"T: {value}", *self.kinetic_options()).returncode, 0)

    def test_invalid_absolute_tolerance_fails(self):
        for atol in (-1, math.nan, math.inf):
            self.assertNotEqual(self.compare("T: 0", "T: 0", *self.kinetic_options(atol)).returncode, 0)


if __name__ == "__main__":
    unittest.main()
