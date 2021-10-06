""" From
https://github.com/eugene-eeo/spheres-from-triangles
Thanks to Eugene Eeo
"""
import numpy as np
import matplotlib.tri as mtri
from geometry.geometry import dummyElementLabel

from tqdm import tqdm
import simplejson as json

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import sys
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D

from collections import namedtuple

class geodesic_grid:
    vdw_radii_standard = {
            dummyElementLabel: 1.0,
            dummyElementLabel.capitalize(): 1.0,
            "H": 1.1,
            "C": 1.7,
            "O": 1.52,
            "N": 1.55,
            "S": 1.55,
            "B": 1.55,
            "BE": 1.55,
            "Be": 1.55,
            "CL": 1.75
            }
    def __init__(self, depth, vdw_radii=vdw_radii_standard, ignoreH = False, radius_all=None):
        self.depth = depth
        if not(radius_all == None):
            self.vdw_radii = {
                    dummyElementLabel:  radius_all,
                    dummyElementLabel.capitalize():  radius_all,
                    "H":  radius_all,
                    "C":  radius_all,
                    "O":  radius_all,
                    "N":  radius_all,
                    "S":  radius_all,
                    "B":  radius_all,
                    "Be":  radius_all,
                    "BE":  radius_all,
                    "CL": radius_all
                    }
        else:
            self.vdw_radii = vdw_radii
        self.ignoreH = ignoreH

Triangle = namedtuple("Triangle", "a,b,c")
Point = namedtuple("Point", "x,y,z")

def get_dict_classifier(grid):
    dict_grid={}
    for pt in tqdm(grid):
        norm = np.linalg.norm(pt)
        norm = "{:.5f}".format(norm)
        if not (norm in dict_grid.keys()):
            dict_grid[norm]=[]
        dict_grid[norm].append(pt)
    return dict_grid

def normalize(p):
    s = sum(u*u for u in p) ** 0.5
    return Point(*(u/s for u in p))


def midpoint(u, v):
    return Point(*((a+b)/2 for a, b in zip(u, v)))


def subdivide_hybrid3(tri, depth):
    def triangle(tri, depth):
        if depth == 0:
            yield tri
            return
        for t in subdivide_centroid(tri, 1):
            yield from edge(t, depth - 1)

    def centroid(tri, depth):
        if depth == 0:
            yield tri
            return
        for t in subdivide_midpoint(tri, 2):
            yield from triangle(t, depth - 1)

    def edge(tri, depth):
        if depth == 0:
            yield tri
            return
        for t in subdivide_edge(tri, 1):
            yield from centroid(t, depth - 1)

    return centroid(tri, depth)


def subdivide_hybrid2(tri, depth):
    def centroid(tri, depth):
        if depth == 0:
            yield tri
            return
        for t in subdivide_centroid(tri, 1):
            yield from edge(t, depth - 1)

    def edge(tri, depth):
        if depth == 0:
            yield tri
            return
        for t in subdivide_edge(tri, 1):
            yield from centroid(t, depth - 1)

    return centroid(tri, depth)


def subdivide_hybrid(tri, depth):
    def centroid(tri, depth):
        if depth == 0:
            yield tri
            return
        for t in subdivide_centroid(tri, 1):
            yield from edge(t, depth - 1)

    def edge(tri, depth):
        if depth == 0:
            yield tri
            return
        for t in subdivide_edge(tri, 1):
            yield from centroid(t, depth - 1)

    return edge(tri, depth)


def subdivide_midpoint2(tri, depth):
    if depth == 0:
        yield tri
        return
    #       p0
    #      /|\
    #     / | \
    #    /  |  \
    #   /___|___\
    # p1   m12   p2
    p0, p1, p2 = tri
    m12 = normalize(midpoint(p1, p2))
    # WRONG TRIANGULATION!
    yield from subdivide_midpoint2(Triangle(p0, m12, p1), depth-1)
    yield from subdivide_midpoint2(Triangle(p0, p2, m12), depth-1)


def subdivide_midpoint(tri, depth):
    if depth == 0:
        yield tri
        return
    #       p0
    #      /|\
    #     / | \
    #    /  |  \
    #   /___|___\
    # p1   m12   p2
    p0, p1, p2 = tri
    m12 = normalize(midpoint(p1, p2))
    yield from subdivide_midpoint(Triangle(m12, p0, p1), depth-1)
    yield from subdivide_midpoint(Triangle(m12, p2, p0), depth-1)


def subdivide_edge(tri, depth):
    if depth == 0:
        yield tri
        return
    #       p0
    #      /  \
    # m01 /....\ m02
    #    / \  / \
    #   /___\/___\
    # p1    m12   p2
    p0, p1, p2 = tri
    m01 = normalize(midpoint(p0, p1))
    m02 = normalize(midpoint(p0, p2))
    m12 = normalize(midpoint(p1, p2))
    triangles = [
        Triangle(p0,  m01, m02),
        Triangle(m01, p1,  m12),
        Triangle(m02, m12, p2),
        Triangle(m01, m02, m12),
    ]
    for t in triangles:
        yield from subdivide_edge(t, depth-1)


def subdivide_centroid(tri, depth):
    if depth == 0:
        yield tri
        return
    #       p0
    #       / \
    #      /   \
    #     /  c  \
    #    /_______\
    #  p1         p2
    p0, p1, p2 = tri
    centroid = normalize(Point(
        (p0.x + p1.x + p2.x) / 3,
        (p0.y + p1.y + p2.y) / 3,
        (p0.z + p1.z + p2.z) / 3,
    ))
    t1 = Triangle(p0, p1, centroid)
    t2 = Triangle(p2, centroid, p0)
    t3 = Triangle(centroid, p1, p2)

    yield from subdivide_centroid(t1, depth - 1)
    yield from subdivide_centroid(t2, depth - 1)
    yield from subdivide_centroid(t3, depth - 1)


def subdivide(faces, depth, method):
    for tri in faces:
        yield from method(tri, depth)

def get_geode(depth=2, method=None):
    if method==None:
        method = subdivide_edge
    # octahedron
    p = 2**0.5 / 2
    faces = [
        # top half
        Triangle(Point(0, 1, 0), Point(-p, 0, p), Point( p, 0, p)),
        Triangle(Point(0, 1, 0), Point( p, 0, p), Point( p, 0,-p)),
        Triangle(Point(0, 1, 0), Point( p, 0,-p), Point(-p, 0,-p)),
        Triangle(Point(0, 1, 0), Point(-p, 0,-p), Point(-p, 0, p)),

        # bottom half
        Triangle(Point(0,-1, 0), Point( p, 0, p), Point(-p, 0, p)),
        Triangle(Point(0,-1, 0), Point( p, 0,-p), Point( p, 0, p)),
        Triangle(Point(0,-1, 0), Point(-p, 0,-p), Point( p, 0,-p)),
        Triangle(Point(0,-1, 0), Point(-p, 0, p), Point(-p, 0,-p)),
    ]

    X = []
    Y = []
    Z = []
    T = []

    for i, tri in enumerate(subdivide(faces, depth, method)):
        X.extend([p.x for p in tri])
        Y.extend([p.y for p in tri])
        Z.extend([p.z for p in tri])
        T.append([3*i, 3*i+1, 3*i+2])

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    T = mtri.Triangulation(X, Y, np.array(T))

    return X, Y, Z, T

def get_geode_points(depth=2, method=None):
    X, Y, Z, T = get_geode(depth, method)
    points=[]
    for x, y, z in zip(X, Y, Z):
        points.append([x, y, z])
    return np.asarray(points)

def generate_geodesic_grid(geom, geodesic_grid, logger, symmetry = False):
    grid = []
    if symmetry:
        pga = geom.getPGA()
        symmetry_operations = pga.get_symmetry_operations()
    for atom in tqdm(geom.atoms+geom.spherecenters+geom.pseudoatoms):
        at    = np.array([ atom['x'], atom['y'], atom['z'] ])
        radius = geodesic_grid.vdw_radii[atom['label']]
        geodesic_points = get_geode_points(geodesic_grid.depth)
        geodesic_points = np.unique(geodesic_points.round(decimals=8),axis=0)
        for pt in geodesic_points:
            #
            # Compute the distance between the point and the current atom
            #
            point = at.copy()
            point[0] = point[0] + pt[0]*radius #translate point coordinates to the atomic repere
            point[1] = point[1] + pt[1]*radius #translate point coordinates to the atomic repere
            point[2] = point[2] + pt[2]*radius #translate point coordinates to the atomic repere
            addAtom = True
            for other_atom in geom.atoms+geom.spherecenters+geom.pseudoatoms:
                #
                # Check that we are looping over other atoms only (not the current one)
                #
                same_atom = False 
                other_at = np.array([ other_atom['x'], other_atom['y'], other_atom['z'] ])
                thrs=0.001
                if (np.abs(at[0]-other_at[0]) < thrs) and (np.abs(at[1]-other_at[1]) < thrs) and (np.abs(at[2]-other_at[2]) < thrs):
                    same_atom = True

                #dist_at_other_at = np.linalg.norm(np.array( [ at[i] - other_at[i] for i in range(3)] ) )
                #if dist_at_other_at < 1e-6:
                #    same_atom = True
                if not(same_atom):
                    #
                    # Compute the distance between the point and the other atom
                    #
                    dist_point_other_at = np.linalg.norm(np.array( [ point[i] - other_at[i] for i in range(3)] ) )
                    #
                    # Get the vdw radius of the other atom
                    #
                    other_radius = geodesic_grid.vdw_radii[other_atom['label']]
                    #
                    # if the point is within the vdw radius of the other atom, skip it
                    #
                    if (dist_point_other_at < other_radius):
                        addAtom = False

            if addAtom and symmetry and len(symmetry_operations)>1:
                for op in symmetry_operations:
                    for i in range(len(grid)):
                        if np.linalg.norm(grid[i]) != np.linalg.norm(point):
                            addAtom = False
                            break
                        if op.are_symmetrically_related(grid[i], point):
                            addAtom = False
                            break
                    if not addAtom:
                        break

            if addAtom:
                logger.debug(
                        "Bq     {0[0]:16.10f} {0[1]:16.10f} {0[2]:16.10f} {1}\n".format(point,point[0]+point[1]))
                grid.append(point)
    return grid

def writegrid(grid):
    np.savetxt("points_values.csv", grid, delimiter=",", header='#x,y,z,v')

if __name__ == '__main__':

    method = {
        "hybrid":   subdivide_hybrid,
        "hybrid2":  subdivide_hybrid2,
        "hybrid3":  subdivide_hybrid3,
        "midpoint": subdivide_midpoint,
        "midpoint2": subdivide_midpoint2,
        "centroid": subdivide_centroid,
        "edge":     subdivide_edge,
        }[sys.argv[1]]
    depth = int(sys.argv[2])
    color = getattr(cm, sys.argv[3] if len(sys.argv) >= 4 else 'coolwarm')

    X, Y, Z, T = get_geode(depth, method)

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(T, Z, lw=0.2, edgecolor="black", color="grey", alpha=1, cmap=color)
    plt.axis('off')
    plt.show()
