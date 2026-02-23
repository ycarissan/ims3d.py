#create unitttest for ims3d.py
import unittest
import os
import tempfile
import sys
import logging
import numpy as np

#import ims3d and ims3d_harv from .. directory
#  we assume we are in IMS3DPATH/tests directory
IMS3DPATH = os.path.dirname(os.path.abspath(__file__)) + '/../'
sys.path.append(IMS3DPATH)
import ims3d
import ims3d_harv
import geometry.geometry as geo_module
import grids.geode as geode_module


class TestGridSymmetryReduction(unittest.TestCase):
    """Tests pour la génération de grilles exploitant la symétrie moléculaire."""

    # Eau : groupe C2v  (O unique + un H unique sur 3 atomes)
    WATER_XYZ = (
        "3\nwater\n"
        "O  0.000000  0.000000  0.117790\n"
        "H  0.000000  0.757440 -0.471160\n"
        "H  0.000000 -0.757440 -0.471160\n"
    )

    NAPHTALENE_XYZ = os.path.join(IMS3DPATH, 'tutorial', 'naphtalene.xyz')

    @staticmethod
    def _write_xyz(tmpdir, content, name="mol.xyz"):
        path = os.path.join(tmpdir, name)
        with open(path, 'w') as f:
            f.write(content)
        return path

    @staticmethod
    def _make_geom(xyz_path):
        return geo_module.Geometry(xyz_path)

    @staticmethod
    def _make_grid(geom, depth=2, restrict=None):
        geo_grid = geode_module.geodesic_grid(depth=depth)
        return geode_module.generate_geodesic_grid(
            geom, geo_grid, None, restrict_to_atom_indices=restrict)

    def _check_no_vdw_overlap(self, geom, grid, depth=2, tol=1e-3):
        """Vérifie qu'aucun point de la grille n'est à l'intérieur d'une sphère VdW."""
        geo_grid  = geode_module.geodesic_grid(depth=depth)
        all_atoms = geom.getAllcenters()
        positions = np.array([[a['x'], a['y'], a['z']] for a in all_atoms])
        radii     = np.array([geo_grid.vdw_radii.get(a['label'], 1.5) for a in all_atoms])

        for pt in np.array(grid):
            dists = np.linalg.norm(positions - pt, axis=1)
            min_clearance = np.min(dists - radii)
            self.assertGreaterEqual(
                min_clearance, -tol,
                f"Point {pt} est à l'intérieur d'une sphère VdW "
                f"(clearance={min_clearance:.4f} Å)")

    # ── Groupe ponctuel ────────────────────────────────────────────────────

    def test_point_group_water(self):
        """Le groupe ponctuel de l'eau doit être C2v."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write_xyz(tmpdir, self.WATER_XYZ)
            self.assertEqual(self._make_geom(path).getPGA().sch_symbol, 'C2v')

    def test_unique_atoms_water(self):
        """L'eau (C2v) a 2 atomes uniques : O et un H."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write_xyz(tmpdir, self.WATER_XYZ)
            unique = self._make_geom(path).getPGA().get_equivalent_atoms()["eq_sets"].keys()
            self.assertEqual(len(unique), 2)

    # ── Filtrage VdW ───────────────────────────────────────────────────────

    def test_vdw_no_overlap_full_grid(self):
        """Aucun point de la grille complète (eau) ne doit être dans une sphère VdW."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write_xyz(tmpdir, self.WATER_XYZ)
            geom = self._make_geom(path)
            self._check_no_vdw_overlap(geom, self._make_grid(geom))

    def test_vdw_no_overlap_sym_grid(self):
        """Aucun point de la grille réduite (eau) ne doit être dans une sphère VdW."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write_xyz(tmpdir, self.WATER_XYZ)
            geom   = self._make_geom(path)
            unique = list(geom.getPGA().get_equivalent_atoms()["eq_sets"].keys())
            self._check_no_vdw_overlap(geom, self._make_grid(geom, restrict=unique))

    def test_vdw_no_overlap_sym_grid_naphthalene(self):
        """Aucun point de la grille réduite (naphtalène) ne doit être dans une sphère VdW."""
        geom   = self._make_geom(self.NAPHTALENE_XYZ)
        unique = list(geom.getPGA().get_equivalent_atoms()["eq_sets"].keys())
        self._check_no_vdw_overlap(geom, self._make_grid(geom, restrict=unique))

    # ── Réduction par symétrie ─────────────────────────────────────────────

    def test_sym_grid_fewer_points_water(self):
        """La grille réduite (eau, C2v) doit avoir moins de points que la grille complète."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path  = self._write_xyz(tmpdir, self.WATER_XYZ)
            geom  = self._make_geom(path)
            unique = list(geom.getPGA().get_equivalent_atoms()["eq_sets"].keys())
            self.assertLess(len(self._make_grid(geom, restrict=unique)),
                            len(self._make_grid(geom)))

    def test_sym_grid_reduction_factor_naphthalene(self):
        """La réduction de grille (naphtalène, D2h) doit être cohérente avec le ratio d'atomes uniques."""
        geom   = self._make_geom(self.NAPHTALENE_XYZ)
        pga    = geom.getPGA()
        unique = list(pga.get_equivalent_atoms()["eq_sets"].keys())
        full_grid = self._make_grid(geom)
        sym_grid  = self._make_grid(geom, restrict=unique)
        ratio_pts      = len(sym_grid) / len(full_grid)
        ratio_expected = len(unique) / len(geom.getAllcenters())
        self.assertLess(ratio_pts, 1.0,
                        "La grille symétrique doit avoir moins de points que la grille complète")
        self.assertAlmostEqual(ratio_pts, ratio_expected, delta=0.15,
                               msg=f"Ratio observé {ratio_pts:.2f} trop éloigné "
                                   f"du ratio attendu {ratio_expected:.2f}")


class TestIMS3D(unittest.TestCase):
    def test_IMS3D_exists(self):
        result = ims3d.main
        self.assertIsNotNone(result)

    def test_IMS3D_tutorial_com_file(self):
        reference_com_file = os.path.join(IMS3DPATH, 'tests', 'references', 'input_batch_00000_ref.com')
        generated_com_file = 'input_batch_00000.com'
        # append the current directory to a filename
        sys.argv = ['ims3d.py', '-r', '1', os.path.join(IMS3DPATH, 'tutorial', 'naphtalene.xyz')]
        #create a temporary directory and move into it
        with tempfile.TemporaryDirectory() as tempdir:
            os.chdir(tempdir)
            # Call the main function of ims3d
            result = ims3d.main()
            # Check if the generated file is similar to the reference file
            generated_content, reference_content = generate_and_compare_contents(reference_com_file, generated_com_file)
            self.assertEqual(generated_content, reference_content)

def generate_and_compare_contents(reference_com_file, generated_com_file):
    with open(generated_com_file, 'r') as generated_file:
        generated_content = generated_file.read()
    with open(reference_com_file, 'r') as reference_file:
        reference_content = reference_file.read()
        # if contents differ, print the differences
    if generated_content != reference_content:
        print("Generated content does not match reference content.")
        print("Generated content:")
        print(generated_content)
        print("Reference content:")
        print(reference_content)
    return generated_content,reference_content

class TestIMS3DHarv(unittest.TestCase):
    def test_IMS3D_harv_exists(self):
        # Example test case, replace with actual test logic
        result = ims3d_harv.main
        self.assertIsNotNone(result)  # Check if the result is not None
    
    def test_IMS3D_tutorial_log_file(self):
        reference_log_file = os.path.join(IMS3DPATH, 'tests', 'references', 'input_batch_00000_ref.log')
        generated_dat_file = 'ims.dat'
        reference_dat_file = os.path.join(IMS3DPATH, 'tests', 'references', 'ims_ref.dat')
        sys.argv = ['ims3d_harv.py', reference_log_file]
        #create a temporary directory and move into it
        with tempfile.TemporaryDirectory() as tempdir:
            os.chdir(tempdir)
            # Call the main function of ims3d
            result = ims3d_harv.main()
            # Check if the generated file is similar to the reference file
            generated_content, reference_content = generate_and_compare_contents(reference_dat_file, generated_dat_file)
            self.assertEqual(generated_content, reference_content)

if __name__ == '__main__':
    unittest.main()
