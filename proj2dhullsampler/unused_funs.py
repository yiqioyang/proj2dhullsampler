

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



