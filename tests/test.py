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
    """Tests pour la génération de grilles exploitant la symétrie moléculaire (naphtalène, D2h)."""

    NAPHTALENE_XYZ = os.path.join(IMS3DPATH, 'tutorial', 'naphtalene.xyz')

    @classmethod
    def setUpClass(cls):
        """Charge la géométrie et pré-calcule les grandeurs communes une seule fois."""
        cls.geom      = geo_module.Geometry(cls.NAPHTALENE_XYZ)
        cls.pga       = cls.geom.getPGA()
        cls.unique    = list(cls.pga.get_equivalent_atoms()["eq_sets"].keys())
        geo_grid      = geode_module.geodesic_grid(depth=2)
        cls.full_grid = geode_module.generate_geodesic_grid(cls.geom, geo_grid, None)
        cls.sym_grid  = geode_module.generate_geodesic_grid(
            cls.geom, geo_grid, None, restrict_to_atom_indices=cls.unique)

    @staticmethod
    def _expand_grid(sym_grid, sym_ops):
        """Étend sym_grid en appliquant les opérations de symétrie, avec déduplication."""
        result = [list(pt) for pt in sym_grid]
        seen = set()
        for pt in sym_grid:
            seen.add((round(pt[0], 5), round(pt[1], 5), round(pt[2], 5)))
        for op in sym_ops:
            for pt in sym_grid:
                coords    = np.array(pt)
                newcoords = op.operate(coords)
                if np.linalg.norm(newcoords - coords) < 1e-6:
                    continue
                key = (round(newcoords[0], 5), round(newcoords[1], 5), round(newcoords[2], 5))
                if key not in seen:
                    seen.add(key)
                    result.append(newcoords.tolist())
        return result

    def _check_no_vdw_overlap(self, grid, tol=1e-3):
        """Vérifie qu'aucun point de la grille n'est à l'intérieur d'une sphère VdW."""
        geo_grid  = geode_module.geodesic_grid(depth=2)
        all_atoms = self.geom.getAllcenters()
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

    def test_point_group(self):
        """Le groupe ponctuel du naphtalène doit être D2h."""
        self.assertEqual(self.pga.sch_symbol, 'D2h')

    def test_unique_atoms(self):
        """Le naphtalène (D2h) a 5 atomes uniques sur 18."""
        self.assertEqual(len(self.unique), 5)

    # ── Filtrage VdW ───────────────────────────────────────────────────────

    def test_vdw_no_overlap_full_grid(self):
        """Aucun point de la grille complète ne doit être dans une sphère VdW."""
        self._check_no_vdw_overlap(self.full_grid)

    def test_vdw_no_overlap_sym_grid(self):
        """Aucun point de la grille réduite ne doit être dans une sphère VdW."""
        self._check_no_vdw_overlap(self.sym_grid)

    # ── Réduction par symétrie ─────────────────────────────────────────────

    def test_sym_grid_fewer_points(self):
        """La grille réduite doit avoir moins de points que la grille complète."""
        self.assertLess(len(self.sym_grid), len(self.full_grid))

    def test_sym_grid_reduction_factor(self):
        """Le ratio de réduction doit être cohérent avec le ratio d'atomes uniques (±15 %)."""
        ratio_pts      = len(self.sym_grid) / len(self.full_grid)
        ratio_expected = len(self.unique) / len(self.geom.getAllcenters())
        self.assertAlmostEqual(ratio_pts, ratio_expected, delta=0.15,
                               msg=f"Ratio observé {ratio_pts:.2f} trop éloigné "
                                   f"du ratio attendu {ratio_expected:.2f}")

    # ── Expansion par symétrie ─────────────────────────────────────────────

    def test_expanded_grid_equals_full(self):
        """Après expansion par symétrie, la grille doit retrouver exactement la grille complète."""
        expanded = self._expand_grid(self.sym_grid, self.pga.get_symmetry_operations())
        self.assertEqual(
            len(expanded), len(self.full_grid),
            f"Grille étendue ({len(expanded)} pts) ≠ grille complète ({len(self.full_grid)} pts)"
        )


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
