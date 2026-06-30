import numpy as np
from scipy.spatial import Delaunay

from .utils import dist_cal



def circum_radius_calculation(pta, ptb, ptc, ptd):
    vol_mat = np.vstack([pta, ptb, ptc, ptd]).transpose()
    vol_mat = np.vstack([vol_mat, np.array([1,1,1,1])])
    det_M = np.linalg.det(vol_mat)
    V = np.abs(det_M)/6

    a = dist_cal(ptd, pta)
    b = dist_cal(ptd, ptb)
    c = dist_cal(ptd, ptc)

    A = dist_cal(ptb, ptc)
    B = dist_cal(pta, ptc)
    C = dist_cal(pta, ptb)

    q_medium = (a * A + b * B + c * C) * (a * A + b * B - c * C) * (a * A - b * B + c * C) * (-a * A + b * B + c * C)
    radius = q_medium ** 0.5 / (24 * V)
    
    return radius




def alpha_shape2d(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    # def add_edge(edges, i, j):
    #     """
    #     Add an edge between the i-th and j-th points,
    #     if not in the list already
    #     """
    #     if (i, j) in edges or (j, i) in edges:
    #         # already added
    #         assert (j, i) in edges, "Can't go twice over same directed edge right?"
    #         if only_outer:
    #             # if both neighboring triangles are in shape, it's not a boundary edge
    #             edges.remove((j, i))
    #         return
    #     edges.add((i, j))

    def add_edge(edges, i, j):
        edge = tuple(sorted((i, j)))
        if edge in edges:
        
            edges.remove(edge)
            return
        edges.add(edge)

    
    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges



def para3_mesh_generator(paras_vars, emu_para, tf_masks):
    ## Might be useful for future development
    para3_meshes = {}
    for k, v in tqdm(paras_vars.items()):
        temp_pts = emu_para[tf_masks[v].all(axis = 1)]
        if temp_pts.shape[0] > 5000:
            temp_pts = temp_pts.sample(5000)
            
        temp_pts_np = temp_pts[list(k)].values
        temp_pts_index = alpha_shape_3D(temp_pts_np, 7)
    
        para3_meshes[k] = trimesh.Trimesh(vertices=temp_pts_np, faces=list(temp_pts_index))

    
    return para3_meshes




def alpha_shape_3D(points, alpha):
    
    # def add_edge(edges, i, j, k):
    #     #perms = (set(p) for p in itertools.permutations([i, j, k]))
    #     #print(perms)
    #     #common_sets = perms.intersection(edges)
    #     #exists = len(common_sets) > 0
        
    #     if (i, j, k) in edges or (i, k, j) in edges or (j, i, k) in edges or (j, k, i) in edges or (k, i, j) in edges or (k, j, i) in edges:
    #         edges.remove((i, k, j))
    #         edges.remove((j, i, k))
    #         edges.remove((j, k, i))
    #         edges.remove((k, i, j))
    #         edges.remove((k, j, i))
    #         return
            
        # edges.add((i, j, k))

    def add_edge(edges, i, j, k):
        edge = tuple(sorted((i, j, k)))
        if edge in edges:    
            edges.remove(edge)
            return
        edges.add(edge)

# #######

    tetra = Delaunay(points)
    # Extract tetrahedra
    edges = set()
    for ia, ib, ic, id in tetra.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        pd = points[id]
        # Compute radius of circumsphere
        circumsphere_radius = circum_radius_calculation(pa, pb, pc, pd)
        if circumsphere_radius < 1.0 / alpha:
            add_edge(edges, ia, ib, ic)
            add_edge(edges, ia, ib, id)
            add_edge(edges, ia, ic, id)
            add_edge(edges, ib, ic, id)
            
    return edges




def dist_cal(pt1, pt2):
    dist = (np.sum((pt1 - pt2) **2))**0.5
    return dist



def circum_center_calculation(pta, ptb, ptc, ptd):
    A = np.vstack([ptb - pta, ptc - pta, ptd - pta])
    B = np.array([np.sum(ptb ** 2) - np.sum(pta ** 2), 
              np.sum(ptc ** 2) - np.sum(pta ** 2), 
              np.sum(ptd ** 2) - np.sum(pta ** 2)]).reshape(3,1)
    B = B * 0.5
    circuncenter_xyz = np.linalg.inv(A) @ B
    return circuncenter_xyz.flatten()


