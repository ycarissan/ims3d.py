import sys
import random
import pickle
from tqdm import tqdm

try :
    import numpy as np
except ModuleNotFoundError as error:
    numpy = None
    print("numpy module not found")

try :
    import scipy.spatial.transform
except ModuleNotFoundError as error:
    scipy = None
    print("scipy module not found")

try :
    import pymatgen
    from pymatgen.io.xyz import XYZ
    import pymatgen.transformations.standard_transformations
except ModuleNotFoundError as error:
    pymatgen = None
    print("pymatgen module not found")

if pymatgen == None or scipy == None:
    print("This should not occur within the proper conda environment and a call with the appropriate python3 interpreter.")
    sys.exit()

dummyElementLabel="Be"

class Geometry:
    def __init__(self, filename, orient = False):
        lines = open(filename, "r").readlines()
        self.header=(lines[1])
        self.atoms = []
        self.pseudoatoms = []
        self.spherecenters = []
        self.pga = None

        for l in lines[2:]:
            a = l.split()
            lbl = a[0].strip().upper()
            position = [float(a[1]), float(a[2]), float(a[3])]
            if lbl=="BQ" or lbl=="X" or lbl=="XX":
                print("BQ found")
                self.spherecenters.append( { 'label': dummyElementLabel, 'x': position[0], 'y': position[1], 'z': position[2] } )
            else:
                self.atoms.append( { 'label': lbl, 'x': position[0], 'y': position[1], 'z': position[2] } )

        if orient:
            filename_atoms_only = self.getgeomfilename_Atomsonly()
            self.geom_original_orientation = pymatgen.Molecule.from_file(filename_atoms_only)
            self.rotvec = get_rotation_vector_to_align_along_z(self.geom_original_orientation)
            
            orientation_rotation = scipy.spatial.transform.Rotation.from_rotvec(self.rotvec)

            for i in range(len(self.spherecenters)):
                el = self.spherecenters[i]
                position = [ el['x'], el['y'], el['z'] ]
                position = orientation_rotation.apply(position)
                self.spherecenters[i]['x'] = position[0]
                self.spherecenters[i]['y'] = position[1]
                self.spherecenters[i]['z'] = position[2]
            for i in range(len(self.atoms)):
                el = self.atoms[i]
                position = [ el['x'], el['y'], el['z'] ]
                position = orientation_rotation.apply(position)
                self.atoms[i]['x'] = position[0]
                self.atoms[i]['y'] = position[1]
                self.atoms[i]['z'] = position[2]

    def getPGA(self):
#        coords = [[0.000000, 0.000000, 0.000000],
#          [0.000000, 0.000000, 1.089000],
#          [1.026719, 0.000000, -0.363000],
#          [-0.513360, -0.889165, -0.363000],
#          [-0.513360, 0.889165, -0.363000]]
#         methane = Molecule(["C", "H", "H", "H", "H"], coords)
        return pymatgen.symmetry.analyzer.PointGroupAnalyzer(self.getPymatgenMolecule())

    def getUniqueElementsByLabels(self):
        unique_indices = self.getPGA().get_equivalent_atoms()["eq_sets"].keys()
        dict_atoms={}
        for lbl in self.getUniqueLabels():
            list_atoms=[]
            for idx in range(len(self.getAtoms())):
                if idx in unique_indices and self.getAtoms()[idx]["label"] == lbl:
                    list_atoms.append(self.getAtoms()[idx])
            dict_atoms[lbl] = list_atoms
        return dict_atoms

    def getElementsByLabels(self):
        dict_atoms={}
        for lbl in self.getUniqueLabels():
            list_atoms=[]
            for at in self.getAtoms():
                if at["label"] == lbl:
                    list_atoms.append(at)
            dict_atoms[lbl] = list_atoms
        return dict_atoms

    def getUniqueLabels(self):
        ulbls = set()
        for at in self.getAtoms():
            ulbls.add(at["label"])
        return ulbls

    def getPymatgenMolecule(self):
        return atoms2Molecule(self.getAllcenters())

    def writePymatgenMolecule(self, filename):
        molecule = self.getPymatgenMolecule()
        molecule.to(filename)

    def writePymatgenMoleculeSymmetryUnique(self, filename):
        equivalent_atoms = self.getPGA().get_equivalent_atoms()
        atoms=[]
        for k in equivalent_atoms["eq_sets"].keys():
            atoms.append(self.getAllcenters()[k])
        molecule = atoms2Molecule(atoms)
        molecule.to(filename=filename)
        return filename

    def getAtoms(self):
        return self.atoms

    def getPseudoAtoms(self):
        return self.pseudoatoms

    def getSphereCenters(self):
        return self.spherecenters

    def getAllcenters(self):
        return self.getAtoms() + self.getPseudoAtoms() + self.getSphereCenters()

    def getAtom(self, index):
        return self.atoms[index]

    def getXYZ(self, index):
        at = self.getAtom(index)
        return [ at['x'], at['y'], at['z'] ]

    def getcoords(self, atomlist):
        """ Return the position of the atoms which determine a cycle """
        coords = []
        for at in atomlist:
            pos = np.asarray(self.getXYZ(at), dtype=np.float64)
            coords.append(pos)
        return coords

    def getBarycenter(self, atomlist):
        coords = self.getcoords(atomlist)
        nat = len(atomlist)
        barycenter = np.asarray([0,0,0])
        for at in coords:
            barycenter = barycenter + at/nat
        return barycenter

    def addPseudoAtom(self, coords):
        self.pseudoatoms.append( { 'label': dummyElementLabel, 'x': coords[0], 'y': coords[1], 'z': coords[2] } )

    def getgeomfilename_Atomsonly(self):
        xyztmp_filename = "tmpfile_{:05d}.xyz".format(int(random.uniform(0, 99999)))
        fio = open(xyztmp_filename, "w+")
        fio.write("{}\n\n".format(len(self.atoms)))
        for atom in self.atoms:
            fio.write("{} {} {} {}\n".format(atom['label'], atom['x'], atom['y'], atom['z']))
        fio.close()
        return xyztmp_filename

def atoms2Molecule(atoms):
    coords=[]
    labels=[]
    for pt in atoms:
        labels.append(pt["label"])
        coords.append([pt["x"], pt["y"], pt["z"]])
#    print(labels)
#    print(coords)
    molecule = pymatgen.core.Molecule(labels, coords)
    return molecule

def get_angle_and_axis(op):
    """Return angle and rotation axis from an symmetry operation"""
    matQ = op.rotation_matrix
    Qxx = matQ[0, 0]
    Qyy = matQ[1, 1]
    Qzz = matQ[2, 2]
    Qzy = matQ[2, 1]
    Qyz = matQ[1, 2]
    Qxz = matQ[0, 2]
    Qzx = matQ[2, 0]
    Qyx = matQ[1, 0]
    Qxy = matQ[0, 1]
    x = Qzy-Qyz
    y = Qxz-Qzx
    z = Qyx-Qxy
    r = np.hypot(x,np.hypot(y,z))
    t = Qxx+Qyy+Qzz
    theta = np.arctan2(r,t-1)
    return theta, np.asarray([x/r, y/r, z/r])


def get_principal_axis(pga):
    theta_min = 2 * np.pi
    axis_min = np.asarray([0, 0, 1])
    for op in pga.get_symmetry_operations():
        theta, axis = get_angle_and_axis(op)
        if theta > np.pi/100 and theta < theta_min:
            theta_min = theta
            axis_min = axis
    return theta_min, axis_min

def get_rotation_vector_to_align_along_z(geom_sym):
    pga = pymatgen.symmetry.analyzer.PointGroupAnalyzer(geom_sym)
    theta, axis = get_principal_axis(pga)
#    print("Principal axis found {0[0]} {0[1]} {0[2]} angle: {1}".format(axis, theta))
    rotation_vector = np.cross([0, 0, 1], axis)
    rotation_vector = rotation_vector / np.linalg.norm(rotation_vector)
    rotation_angle = np.arcsin(np.linalg.norm(rotation_vector)/np.linalg.norm(axis))
    rotation_vector = rotation_vector * rotation_angle
#    rot = pymatgen.transformations.standard_transformations.RotationTransformation(rotation_axis, rotation_angle, angle_in_radians=True)
    return rotation_vector

def applySymmOps(sym_ops, points):
    for op in sym_ops:
        pts = op.operate_multi(points[:])
        points = np.concatenate([points, pts], axis=0)
    return points

def applySymmOps_onGrid(sym_ops, grid):
    toadd=[]
    for op in sym_ops:
        generated=np.array([])
        for i in tqdm(range(len(grid))):
            pt = grid[i]
            coords = np.array([pt["x"], pt["y"], pt["z"]])
            val=pt['ims']
            newcoords = op.operate(coords)
            if np.linalg.norm(newcoords-coords)>0.0001:
                generated = np.append(generated, {'label': 'Bq', 'x': newcoords[0], 'y': newcoords[1], 'z': newcoords[2], 'ims': val})
        if len(generated)>0:
            toadd.append(generated)
    for generated in toadd:
        grid = np.append(grid, generated)
    return grid

def readSymmOps():
    try:
        with open("symmetry_operations.bin","rb") as fio:
            sym_ops = pickle.load(fio)
            fio.close()
            return sym_ops
    except:
        return None

def main():
    geom = Geometry("methane.xyz")
    geom.writePymatgenMoleculeSymmetryUnique("meth.xyz")

if __name__=="__main__":
    main()
