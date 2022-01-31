# Core imports
import KratosMultiphysics
from KratosMultiphysics import KratosUnittest

# HDF5 imports
from KratosMultiphysics.HDF5Application import ModelPartPattern

# STL imports
import pathlib
import shutil


class TestModelPartPattern(KratosUnittest.TestCase):
    def setUp(self):
        """
        Create the following directory structure relative to the working directory:
        .
        └── test_path_pattern
            ├── valid_root
                └── mdpa name_step_1_time_0.0.h5
                    └── mdpa name_step_1_time_10.h5
                └── mdpa name_step_1_time_1.0.h5
                └── mdpa name_step_1_time_1.h5
                └── mdpa name_step_1_time_1.1.h5
                └── mdpa name_step_2_time_1.0.h5
                └── mdpa name_step_3_time_0.5.h5
                └── mdpa name_step_3_time_2.h5
                └── mdpa name_step_3_time_03.h5
                └── mdpa name_step_3_time_004.0.h5
                └── mdpa_name_step_1_time_1.h5
                └── m_dpa_name_step_1_time_1.h5
                └── mdpa name_step_04_time_1.h5
                └── mdpa name_step_04_time_2.hdf5
                └── ir.relevant
            └── invalid_root
                └── mdpa name_step_1_time_0.0.h5
                    └── mdpa name_step_1_time_10.h5
                └── mdpa name_step_1_time_1.0.h5
                └── mdpa name_step_1_time_1.h5
                └── mdpa name_step_1_time_1.1.h5
                └── mdpa name_step_2_time_1.0.h5
                └── mdpa name_step_3_time_0.5.h5
                └── mdpa name_step_3_time_2.h5
                └── mdpa name_step_3_time_03.h5
                └── mdpa name_step_3_time_004.0.h5
                └── mdpa_name_step_1_time_1.h5
                └── m_dpa_name_step_1_time_1.h5
                └── mdpa name_step_04_time_1.h5
                └── mdpa name_step_04_time_2.hdf5
                └── ir.relevant
            └── step_1
                └── mdpa name_time_1.h5
            └── step_2
                └── mdpa name_time_0.h5
        """
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting("test_path_pattern")
        mkdir_arguments = {"exist_ok" : False, "parents" : True}
        self.test_root = pathlib.Path().cwd() / "test_path_pattern"
        self.valid_root = self.test_root / "valid_root"
        self.invalid_root = self.test_root / "invalid_root"

        file_names = [
            "mdpa name_step_1_time_1.0.h5",
            "mdpa name_step_1_time_1.h5",
            "mdpa name_step_1_time_1.1.h5",
            "mdpa name_step_2_time_1.0.h5",
            "mdpa name_step_3_time_0.5.h5",
            "mdpa name_step_3_time_2.h5",
            "mdpa name_step_3_time_03.h5",
            "mdpa name_step_3_time_004.0.h5",
            "mdpa_name_step_1_time_1.h5",
            "m_dpa_name_step_1_time_1.h5",
            "mdpa name_step_04_time_1.h5",
            "mdpa name_step_04_time_2.hdf5",
            "ir.relevant"
        ]

        for directory in (self.test_root, self.valid_root, self.invalid_root):
            (directory / "mdpa name_step_1_time_0.0.h5").mkdir(**mkdir_arguments)
            (directory / "mdpa name_step_1_time_0.0.h5" / "mdpa name_step_1_time_10.h5").touch(exist_ok=False)

        for file_name in file_names:
            (self.valid_root / file_name).touch(exist_ok=False)
            (self.invalid_root / file_name).touch(exist_ok=False)

        (self.test_root / "step_1").mkdir(**mkdir_arguments)
        (self.test_root / "step_1" / "mdpa name_time_1.h5").touch(exist_ok=False)
        (self.test_root / "step_2").mkdir(**mkdir_arguments)
        (self.test_root / "step_2" / "mdpa name_time_0.h5").touch(exist_ok=False)

    def tearDown(self):
        """Remove ./test_path_pattern"""
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting("test_path_pattern")

    def test_IsAMatch(self):
        pattern_string = "mdpa_name_<model_part_name>_step_<step>_time_<time>_suffix"
        pattern = ModelPartPattern(pattern_string)

        valid_strings = [
            "mdpa_name_NAME_step_1_time_1.0_suffix",
            "mdpa_name_N_step_1_time_1.0_suffix",
            "mdpa_name_NAME_step_01_time_1.0_suffix",
            "mdpa_name_NAME_step_123456789_time_1.0_suffix",
            "mdpa_name_NAME_step_0101_time_1.0_suffix",
            "mdpa_name_NAME_step_1_time_1_suffix",
            "mdpa_name_NAME_step_1_time_-1.0_suffix",
            "mdpa_name_NAME_step_1_time_-1_suffix",
            "mdpa_name_NAME_step_1_time_0.123456789_suffix"
        ]

        invalid_strings = [
            "mdpa_name__step_1_time_1.0_suffix",
            "mdpa_name_NAME_step__time_1.0_suffix",
            "mdpa_name_NAME_step_1O_time_1.0_suffix",
            "mdpa_name_NAME_step_1a1_time_1.0_suffix",
            "mdpa_name_NAME_step_1_time__suffix",
            "mdpa_name_NAME_step_1_time_1.O_suffix",
            "mdpa_name_NAME_step_1_time_1.0_suffiy",
            "Mdpa_name_NAME_step_1_time_1.0_suffix"
        ]

        for string in valid_strings:
            self.assertTrue(pattern.IsAMatch(string))

        for string in invalid_strings:
            self.assertFalse(pattern.IsAMatch(string))

    def test_Match(self):
        pattern_string = "<model_part_name>/name_<model_part_name>_step_<step>_time_<time>suffix"
        pattern = ModelPartPattern(pattern_string)

        placeholders = ("<model_part_name>", "<step>", "<time>")
        values = [
            dict(zip(placeholders, ("mdpa name", "1", "1.0"))),
            dict(zip(placeholders, ("mdpa name", "01", "1"))),
            dict(zip(placeholders, (r'()[]\\"-+_*' + "'", "090", "0.010")))
        ]

        for dictionary in values:
            string = pattern_string
            for key, value in dictionary.items():
                string = string.replace(key, value)
            matches = pattern.Match(string)

            for key, value in dictionary.items():
                self.assertIn(key, matches)
                self.assertTrue(matches[key])
                for match_value in matches[key]:
                    self.assertEqual(match_value, value)

            self.assertEqual(len(matches["<model_part_name>"]), 2)

    def test_Glob(self):
        valid_paths = {
            self.valid_root / "mdpa name_step_1_time_0.0.h5",
            self.valid_root / "mdpa name_step_1_time_1.h5",
            self.valid_root / "mdpa name_step_1_time_1.0.h5",
            self.valid_root / "mdpa name_step_1_time_1.1.h5",
            self.valid_root / "mdpa name_step_2_time_1.0.h5",
            self.valid_root / "mdpa name_step_3_time_0.5.h5",
            self.valid_root / "mdpa name_step_3_time_2.h5",
            self.valid_root / "mdpa name_step_3_time_03.h5",
            self.valid_root / "mdpa name_step_3_time_004.0.h5",
            self.valid_root / "mdpa name_step_04_time_1.h5",
        }

        invalid_paths = {
            self.valid_root / "mdpa_name_step_1_time_1.h5",
            self.valid_root / "m_dpa_name_step_1_time_1.h5",
            self.valid_root / "mdpa name_step_04_time_2.hdf5",
            self.valid_root / "ir.relevant",
            self.invalid_root / "mdpa name_step_1_time_1.0.h5",
            self.invalid_root / "mdpa name_step_1_time_1.h5",
            self.invalid_root / "mdpa name_step_1_time_1.1.h5",
            self.invalid_root / "mdpa name_step_2_time_1.0.h5",
            self.invalid_root / "mdpa name_step_3_time_0.5.h5",
            self.invalid_root / "mdpa name_step_3_time_2.h5",
            self.invalid_root / "mdpa name_step_3_time_03.h5",
            self.invalid_root / "mdpa name_step_3_time_004.0.h5",
            self.invalid_root / "ir.relevant",
            self.test_root / "mdpa name_step_1_time_0.0.h5" / "mdpa name_step_1_time_0.0.h5"
        }

        pattern_string = "test_path_pattern/valid_root/<model_part_name>_step_<step>_time_<time>.h5"

        # Base pattern (all placeholders)
        pattern = ModelPartPattern(pattern_string)
        paths = [pathlib.Path(path) for path in pattern.Glob()]
        for path in valid_paths:
            self.assertIn(path, paths)

        extra = {
            self.valid_root / "mdpa_name_step_1_time_1.h5",
            self.valid_root / "m_dpa_name_step_1_time_1.h5"
        }
        self.assertEqual(len(paths), len(valid_paths) + len(extra))
        for path in valid_paths | extra:
            self.assertIn(path, paths)
        for path in invalid_paths - extra:
            self.assertNotIn(path, paths)

        # Partial pattern (<model_part_name> replaced)
        pattern = ModelPartPattern(pattern_string.replace("<model_part_name>", "mdpa name"))
        paths = [pathlib.Path(path) for path in pattern.Glob()]
        self.assertEqual(len(paths), len(valid_paths))
        for path in valid_paths:
            self.assertIn(path, paths)
        for path in invalid_paths:
            self.assertNotIn(path, paths)


if __name__ == "__main__":
    KratosUnittest.main()
