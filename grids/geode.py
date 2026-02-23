""" From
https://github.com/eugene-eeo/spheres-from-triangles
Thanks to Eugene Eeo
"""
import numpy as np
import matplotlib.tri as mtri
from geometry.geometry import dummyElementLabel

from tqdm import tqdm

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import sys
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D

from collections import namedtuple

try:
    from numba import jit
except ModuleNotFoundError as error:
    print("Module numba not found: the conda environment needs update.\nRunnig the following command will fix the problem:\n\nconda env update --file ${IMS3D_HOME}/conda/ims3d_conda_env.yml --prune\n")
    quit()


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
            "CL": 1.75,
            "F": 1.75
            }
    def __init__(self, depth, vdw_radii=vdw_radii_standard, ignoreH=False, radius_all=None):
        self.depth = depth
        if radius_all is not None:
            self.vdw_radii = {k: radius_all for k in self.vdw_radii_standard}
        else:
            self.vdw_radii = vdw_radii
        self.ignoreH = ignoreH


Triangle = namedtuple("Triangle", "a,b,c")
Point    = namedtuple("Point",    "x,y,z")


@jit(nopython=True, cache=True)
def optimized_norm(pt):
    return np.linalg.norm(pt)


def get_dict_classifier_key(pt):
    norm = optimized_norm(np.array(pt))
    return "{:.5f}".format(norm)


def get_dict_classifier(grid):
    dict_grid = {}
    for pt in grid:
        norm = get_dict_classifier_key(pt)
        if norm not in dict_grid:
            dict_grid[norm] = []
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
    p0, p1, p2 = tri
    m12 = normalize(midpoint(p1, p2))
    yield from subdivide_midpoint2(Triangle(p0, m12, p1), depth-1)
    yield from subdivide_midpoint2(Triangle(p0, p2, m12), depth-1)


def subdivide_midpoint(tri, depth):
    if depth == 0:
        yield tri
        return
    p0, p1, p2 = tri
    m12 = normalize(midpoint(p1, p2))
    yield from subdivide_midpoint(Triangle(m12, p0, p1), depth-1)
    yield from subdivide_midpoint(Triangle(m12, p2, p0), depth-1)


def subdivide_edge(tri, depth):
    if depth == 0:
        yield tri
        return
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
    p0, p1, p2 = tri
    centroid = normalize(Point(
        (p0.x + p1.x + p2.x) / 3,
        (p0.y + p1.y + p2.y) / 3,
        (p0.z + p1.z + p2.z) / 3,
    ))
    yield from subdivide_centroid(Triangle(p0, p1, centroid),   depth - 1)
    yield from subdivide_centroid(Triangle(p2, centroid, p0),   depth - 1)
    yield from subdivide_centroid(Triangle(centroid, p1, p2),   depth - 1)


def subdivide(faces, depth, method):
    for tri in faces:
        yield from method(tri, depth)


def get_geode(depth=2, method=None):
    if method is None:
        method = subdivide_edge
    p = 2**0.5 / 2
    faces = [
        Triangle(Point(0, 1, 0), Point(-p, 0, p), Point( p, 0, p)),
        Triangle(Point(0, 1, 0), Point( p, 0, p), Point( p, 0,-p)),
        Triangle(Point(0, 1, 0), Point( p, 0,-p), Point(-p, 0,-p)),
        Triangle(Point(0, 1, 0), Point(-p, 0,-p), Point(-p, 0, p)),
        Triangle(Point(0,-1, 0), Point( p, 0, p), Point(-p, 0, p)),
        Triangle(Point(0,-1, 0), Point( p, 0,-p), Point( p, 0, p)),
        Triangle(Point(0,-1, 0), Point(-p, 0,-p), Point( p, 0,-p)),
        Triangle(Point(0,-1, 0), Point(-p, 0, p), Point(-p, 0,-p)),
    ]

    X, Y, Z, T_list = [], [], [], []
    for i, tri in enumerate(subdivide(faces, depth, method)):
        X.extend([p.x for p in tri])
        Y.extend([p.y for p in tri])
        Z.extend([p.z for p in tri])
        T_list.append([3*i, 3*i+1, 3*i+2])

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    T = mtri.Triangulation(X, Y, np.array(T_list))

    return X, Y, Z, T


def get_geode_points(depth=2, method=None):
    X, Y, Z, T = get_geode(depth, method)
    points = np.stack([X, Y, Z], axis=1)
    return points


def generate_geodesic_grid(geom, geodesic_grid, logger, symmetry=False):
    grid = []
    all_atoms = geom.atoms + geom.spherecenters + geom.pseudoatoms

    if symmetry:
        pga = geom.getPGA()
        symmetry_operations = pga.get_symmetry_operations()

    # ── Calcul unique de la géodésique (avant la boucle) ──────────────────
    geodesic_points = get_geode_points(geodesic_grid.depth)
    geodesic_points = np.unique(geodesic_points.round(decimals=8), axis=0)
    # shape: (P, 3)

    # ── Pré-calcul des positions et rayons de tous les atomes ─────────────
    all_positions = np.array([[a['x'], a['y'], a['z']] for a in all_atoms])
    all_radii     = np.array([geodesic_grid.vdw_radii[a['label']] for a in all_atoms])

    for idx, atom in enumerate(tqdm(all_atoms)):
        at     = all_positions[idx]
        radius = all_radii[idx]

        # Translation + mise à l'échelle vectorisée — shape: (P, 3)
        tmpgrid = at + geodesic_points * radius

        # ── Filtrage vectorisé ────────────────────────────────────────────
        other_mask      = np.ones(len(all_atoms), dtype=bool)
        other_mask[idx] = False
        other_positions = all_positions[other_mask]  # (M-1, 3)
        other_radii     = all_radii[other_mask]       # (M-1,)

        if len(other_positions) > 0:
            # diff : (P, M-1, 3) → dists : (P, M-1)
            diff  = tmpgrid[:, np.newaxis, :] - other_positions[np.newaxis, :, :]
            dists = np.linalg.norm(diff, axis=2)
            valid = np.all(dists >= other_radii, axis=1)  # (P,)
        else:
            valid = np.ones(len(tmpgrid), dtype=bool)

        valid_points = tmpgrid[valid]

        # ── Symétrie (inchangée, rarement utilisée) ───────────────────────
        if symmetry:
            kept = []
            generated_by_symmetry = set()
            for i in range(len(valid_points)):
                if i in generated_by_symmetry:
                    continue
                kept.append(valid_points[i])
                k = get_dict_classifier_key(valid_points[i])
                for j in range(i + 1, len(valid_points)):
                    if j not in generated_by_symmetry:
                        if get_dict_classifier_key(valid_points[j]) == k:
                            for op in symmetry_operations:
                                if op.are_symmetrically_related(valid_points[i], valid_points[j]):
                                    generated_by_symmetry.add(j)
                                    break
            grid.extend(kept)
        else:
            grid.extend(valid_points.tolist())

    return grid


def writegrid(grid):
    np.savetxt("points_values.csv", grid, delimiter=",", header='#x,y,z,v')


if __name__ == '__main__':

    method = {
        "hybrid":    subdivide_hybrid,
        "hybrid2":   subdivide_hybrid2,
        "hybrid3":   subdivide_hybrid3,
        "midpoint":  subdivide_midpoint,
        "midpoint2": subdivide_midpoint2,
        "centroid":  subdivide_centroid,
        "edge":      subdivide_edge,
        }[sys.argv[1]]
    depth = int(sys.argv[2])
    color = getattr(cm, sys.argv[3] if len(sys.argv) >= 4 else 'coolwarm')

    X, Y, Z, T = get_geode(depth, method)

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(T, Z, lw=0.2, edgecolor="black", color="grey", alpha=1, cmap=color)
    plt.axis('off')
    plt.show()
