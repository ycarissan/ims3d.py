import numpy as np

class angular_grid:
    vdw_radii_standard = {
            "E": 1.0, # pseudo atom
            "H": 1.1,
            "C": 1.7,
            "O": 1.52,
            "N": 1.55,
            "Cl": 1.75,
            }
    def __init__(self, ntheta=12, vdw_radii=vdw_radii_standard, ignoreH=False, radius_all=None):
        self.ntheta = ntheta
        if radius_all is not None:
            self.vdw_radii = {k: radius_all for k in self.vdw_radii_standard}
        else:
            self.vdw_radii = vdw_radii
        self.ignoreH = ignoreH


def _build_unit_sphere_offsets(ntheta):
    """
    Pré-calcule les vecteurs unitaires de la grille angulaire (shape: (P, 3)).
    Appelé une seule fois, indépendant des atomes.
    """
    thetas = np.linspace(0,     np.pi,   ntheta,   endpoint=False)
    phis   = np.linspace(0, 2 * np.pi, 2 * ntheta, endpoint=False)

    # Grille 2D → tableau plat de (P, 3)
    TH, PH = np.meshgrid(thetas, phis, indexing='ij')
    TH = TH.ravel()
    PH = PH.ravel()

    offsets = np.stack([
        np.sin(TH) * np.cos(PH),
        np.sin(TH) * np.sin(PH),
        np.cos(TH),
    ], axis=1)  # shape: (P, 3)
    return offsets


def generate_angular_grid(geom, angular_grid, logger):
    """
    Génère la grille 3D de points fantômes (Bq) autour de chaque atome.
    Version vectorisée : les boucles Python sur les points/atomes sont
    remplacées par des opérations matricielles NumPy.
    """
    ntheta = angular_grid.ntheta

    # ── Filtrage des atomes à traiter ─────────────────────────────────────
    all_atoms = geom.atoms + geom.spherecenters + geom.pseudoatoms
    if angular_grid.ignoreH:
        active_atoms = [a for a in all_atoms if a['label'] != 'H']
    else:
        active_atoms = all_atoms

    if not active_atoms:
        return [], []

    # ── Pré-calcul des positions et rayons ────────────────────────────────
    all_positions = np.array([[a['x'], a['y'], a['z']] for a in active_atoms])
    all_radii     = np.array([angular_grid.vdw_radii[a['label']] for a in active_atoms])

    # ── Vecteurs unitaires de la sphère (calculés une seule fois) ─────────
    unit_offsets = _build_unit_sphere_offsets(ntheta)  # shape: (P, 3)
    P = len(unit_offsets)

    grid    = []
    normals = []

    for idx, atom in enumerate(active_atoms):
        at     = all_positions[idx]       # shape: (3,)
        radius = all_radii[idx]

        # Points candidats autour de cet atome : shape (P, 3)
        normal_vecs = unit_offsets * radius          # mise à l'échelle
        points      = at + normal_vecs               # translation

        # ── Filtrage vectorisé ────────────────────────────────────────────
        # Distances de tous les P points vers tous les M autres atomes
        # other_positions : (M-1, 3)  /  other_radii : (M-1,)
        other_mask      = np.ones(len(active_atoms), dtype=bool)
        other_mask[idx] = False
        other_positions = all_positions[other_mask]   # (M-1, 3)
        other_radii     = all_radii[other_mask]        # (M-1,)

        if len(other_positions) > 0:
            # diff : (P, M-1, 3)
            diff  = points[:, np.newaxis, :] - other_positions[np.newaxis, :, :]
            dists = np.linalg.norm(diff, axis=2)          # (P, M-1)
            valid = np.all(dists >= other_radii, axis=1)  # (P,)
        else:
            valid = np.ones(P, dtype=bool)

        # ── Collecte des points valides ───────────────────────────────────
        for i in np.where(valid)[0]:
            pt  = points[i]
            nv  = normal_vecs[i]
            logger.debug(
                "Bq     {0[0]:16.10f} {0[1]:16.10f} {0[2]:16.10f}\n".format(pt))
            grid.append(pt)
            normals.append(list(nv))

    return grid, normals


def writegrid(grid, normals=None):
    fio_pts  = open("points_values.csv", "w+")
    fio_norm = open("normals.csv", "w+")
    fio_pts.write("#x,y,z,v\n")
    fio_norm.write("#nx,ny,nz\n")
    for ipt in range(len(grid)):
        pt = grid[ipt]
        nv = normals[ipt]
        fio_pts.write(
            "{:12.8f}, {:12.8f}, {:12.8f}, {:12.8f}, {:12.8f}, {:12.8f}\n".format(
                pt[0], pt[1], pt[2], pt[0]+pt[1], pt[0]+pt[1], pt[0]+pt[1]))
        fio_norm.write(
            "{:12.8f}, {:12.8f}, {:12.8f}\n".format(nv[0], nv[1], nv[2]))


def readgrid():
    grid    = np.loadtxt("points_values.csv", delimiter=",", skiprows=1)
    normals = np.loadtxt("normals.csv",       delimiter=",", skiprows=1)
    return grid, normals


def addNormals(ims_grid, grid, normals):
    for i in range(len(ims_grid)):
        samepoint = np.all(np.isclose(
            [ims_grid[i]['x'], ims_grid[i]['y'], ims_grid[i]['z']],
            grid[i][:3], atol=1e-4))
        if samepoint:
            ims_grid[i]['nx'] = normals[i][0]
            ims_grid[i]['ny'] = normals[i][1]
            ims_grid[i]['nz'] = normals[i][2]
        else:
            print("{} != {}".format(ims_grid[i]['x'], grid[i][0]))
            print("{} != {}".format(ims_grid[i]['y'], grid[i][1]))
            print("{} != {}".format(ims_grid[i]['z'], grid[i][2]))
            exit(99)
