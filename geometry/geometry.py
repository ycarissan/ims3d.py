import sys
import random

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
    import pymatgen.transformations.standard_transformations
except ModuleNotFoundError as error:
    pymatgen = None
    print("pymatgen module not found")

if pymatgen == None or scipy == None:
    print("This should not occur within the proper conda environment and a call with the appropriate python3 interpreter.")
    sys.exit()

class Geometry:
    def __init__(self, filename, orient = False):
        lines = open(filename, "r").readlines()
        self.header=(lines[1])
        self.atoms = []
        self.pseudoatoms = []
        self.spherecenters = []

        for l in lines[2:]:
            a = l.split()
            lbl = a[0].strip().upper()
            position = [float(a[1]), float(a[2]), float(a[3])]
            if lbl=="BQ" or lbl=="X" or lbl=="XX":
                print("BQ found")
                self.spherecenters.append( { 'label': "E", 'x': position[0], 'y': position[1], 'z': position[2] } )
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

    def getAllAtomCoords(self):
        coords=[]
        for iat in range(len(self.atoms)):
            coords.append(self.getXYZ(iat))
        return np.array(coords)

    def getAtom(self, index):
        return self.atoms[index]

    def getXYZ(self, index):
        at = self.getAtom(index)
        return [ at['x'], at['y'], at['z'] ]

    def getDistance(self, i,j):
        x0 = self.getXYZ(i)[0]
        y0 = self.getXYZ(i)[1]
        z0 = self.getXYZ(i)[2]
        x1 = self.getXYZ(j)[0]
        y1 = self.getXYZ(j)[1]
        z1 = self.getXYZ(j)[2]
        return np.linalg.norm([x0-x1, y0-y1, z0-z1])

    def getcoords(self, atomlist):
        """ Return the position of the atoms in the list given as argument """
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
        self.pseudoatoms.append( { 'label': 'E', 'x': coords[0], 'y': coords[1], 'z': coords[2] } )

    def getgeomfilename_Atomsonly(self):
        xyztmp_filename = "tmpfile_{:05d}.xyz".format(int(random.uniform(0, 99999)))
        fio = open(xyztmp_filename, "w+")
        fio.write("{}\n\n".format(len(self.atoms)))
        for atom in self.atoms:
            fio.write("{} {} {} {}\n".format(atom['label'], atom['x'], atom['y'], atom['z']))
        fio.close()
        return xyztmp_filename

#naive implementation    def getThreeNearestNeighbourgsIndices(self, index):
#naive implementation        dist=[]
#naive implementation        for i in range(len(self.atoms)):
#naive implementation            dist.append(self.getDistance(index, i))
#naive implementation        dist = np.array(dist)
#naive implementation        order = dist.argsort()
#naive implementation        ranks = order.argsort()
#naive implementation        return np.where((ranks > 0) & (ranks < 4))[0]

    def getThreeNearestNeighbourgsIndices(self, index):
        """ see https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates """
        xyz = self.getXYZ(index)
        coords = self.getAllAtomCoords()
        return scipy.spatial.KDTree(coords).query(xyz, [2, 3, 4])[1]

    def getNormalToThreeNearestNeighbougsPlane(self, index):
        i, j, k = self.getThreeNearestNeighbourgsIndices(index)
        A = np.array(self.getXYZ(i))
        B = np.array(self.getXYZ(j))
        C = np.array(self.getXYZ(k))
        AB=B-A
        AC=C-A
        normal_to_plane = np.cross(AB, AC)
        if normal_to_plane[2] < 0:
            normal_to_plane = -1 * normal_to_plane
        return normal_to_plane

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

def main():

    pga = pymatgen.symmetry.analyzer.PointGroupAnalyzer(newgeom)
    theta, axis = get_principal_axis(pga)
    print("Principal axis found {0[0]} {0[1]} {0[2]} angle: {1}".format(axis, theta))

if __name__=="__main__":
    main()
