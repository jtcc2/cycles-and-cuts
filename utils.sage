def get_start_curve(Fp2):
    """Gives the j-invariant of a supersingular elliptic curve over Fp2"""
    p = Fp2.characteristic()
    q = next(q for q in Primes() if q%4 == 3 and kronecker_symbol(-q,p) == -1)
    K = QuadraticField(-q)
    H = K.hilbert_class_polynomial()
    j0 = H.change_ring(Fp2).any_root()
    return j0

def get_L_graph(p, L_list):
    """Finds the L-isogeny graph over Fpbar, for the given list of L's.
    Note that each edge is (start_vertex, end_vertex, label) where the we're using the label to store the degree of the isogeny."""
    Fp2 = GF(p**2)
    j = get_start_curve(Fp2)
    E = EllipticCurve_from_j(j)
    edges = []
    for ell in L_list:
        G_ell = E.isogeny_ell_graph(ell, directed=True, label_by_j=True)
        for edge in G_ell.edges():
            edges.append((edge[0], edge[1], ell))
    return DiGraph(edges,loops=True,multiedges=True)

def vertex_complement(G, C):
    """Return the vertex complement of a subset of vertices."""
    comp = G.vertices(vertex_property = (lambda v : v not in C))
    return comp

def volume(G, C):
    """
    Computes volume of subset of vertices.
    For d-regular graphs this is just d times the size of C.
    Although we include this to handle isogeny graphs with special curves (p not 1 mod 12), which are not regular at those special curves.
    """
    return sum([G.out_degree(v) for v in C])

def edge_expansion(G, C):
    """Compute edge expansion of cut C in graph G."""
    Ccomp = vertex_complement(G, C)
    num_edges_between = len(G.edge_boundary(C, Ccomp))
    return num_edges_between / min(volume(G, C), volume(G, Ccomp))

def get_spectral_ordering(G, use_numpy_over_sage=False):
    """Returns spectral ordering of vertices from Fiedler vector"""
    L = G.laplacian_matrix()
    # Get eigenvalues/vectors
    if use_numpy_over_sage:
        import numpy as np
        from numpy.linalg import eig
        a = np.array(L.rows())
        w,v=eig(a)
        es = list(zip(w, v))
    else:
        es_nonflat = L.eigenvectors_right()
        es = []
        for e in es_nonflat:
            for vec in e[1]:
                es.append((e[0], vec))
    # Find fiedler vector
    def first_entry(e):
        return e[0]
    es.sort(key=first_entry, reverse=True)
    lambda2 = es[1][0]
    fiedler_vector = list(es[1][1])
    # Order vertices by entry in fiedler vector
    vs_and_es = list(zip(fiedler_vector, G.vertices()))    # list of tuples [(eigenvector_entry, vertex_label), (eigenvector_entry, vertex_label), ...]
    vs_and_es.sort(key=first_entry, reverse=True)
    return [v[1] for v in vs_and_es]

def fiedler(G, use_numpy_over_sage=False, ordering=None):
    """Executes the min edge expansion, sweep cut, version of Fiedler's algorithm. Returns edge expansion and cut."""
    if ordering == None:
        ordering = get_spectral_ordering(G, use_numpy_over_sage)
    # Iterate over sweep cuts
    best_cut_so_far = None
    best_edge_expansion_so_far = None
    for i in range(1, len(G.vertices())):
        Ci = ordering[:i]
        Ci_comp = vertex_complement(G, Ci)
        # Find sweep cut with minimal edge expansion, or more precisely  min( h(C_i) ,  h(V not C_i) )
        Ci_expansion = edge_expansion(G, Ci)
        Ci_comp_expansion = edge_expansion(G, Ci_comp)
        exp_max = max(Ci_expansion, Ci_comp_expansion)
        if best_edge_expansion_so_far == None or exp_max < best_edge_expansion_so_far:
            best_edge_expansion_so_far = exp_max
            best_cut_so_far = Ci
    return best_edge_expansion_so_far, best_cut_so_far

def get_neighbour_ordering(G, start_set):
    """Given graph  `G` returns the list of vertices ordered with respect to the neighbour ordering staring from vertices in `start_set` list."""
    ordered_vertices = copy(start_set)
    left_to_add_neighbors = copy(start_set)
    while len(ordered_vertices) < len(G.vertices()):
        expand_vertex = left_to_add_neighbors[0]
        left_to_add_neighbors.remove(expand_vertex)
        for n in G.neighbors(expand_vertex):
            if n not in ordered_vertices:
                ordered_vertices.append(n)
                left_to_add_neighbors.append(n)
    return ordered_vertices

def get_greedy_ordering(G, start_set):
    """Given graph  `G` returns the list of vertices ordered with respect to the greedy neighbour ordering staring from vertices in `start_set` list."""

    def phi(G, C):
        C_comp = vertex_complement(G, C)
        C_expansion = edge_expansion(G, C)
        C_comp_expansion = edge_expansion(G, C_comp)
        return max(C_expansion, C_comp_expansion)

    ordered_vertices = copy(start_set)
    cut_neighbors = []
    for v in ordered_vertices:
        for n in G.neighbors(v):
            if n not in cut_neighbors:
                cut_neighbors.append(n)
    while len(ordered_vertices) < len(G.vertices()) - 1:
        # iterate over cut_neighbors and add one which gives smallest edge expansion
        minnum = 2
        next_neighbor = None
        for neighbor in cut_neighbors:
            cut = copy(ordered_vertices)
            cut.append(neighbor)
            expansion = phi(G, cut)
            if expansion < minnum:
                minnum = expansion
                next_neighbor = neighbor    
        # add it to ordered_vertices
        ordered_vertices.append(next_neighbor)
        # add its neighbors to cut_neighbors, and remove it from cut_neighbors
        cut_neighbors.remove(next_neighbor)
        for n in G.neighbors(next_neighbor):
            if n not in cut_neighbors:
                if n not in ordered_vertices:
                    cut_neighbors.append(n)
    # Add last vertex
    for v in G.vertices():
        if v not in ordered_vertices:
            ordered_vertices.append(v)
            break
    return ordered_vertices

def quadratic_ideals_above_prime(O, ell):
    """
    Returns list of quadratic ideals above a prime. Returns empty list if inert, list of 1 ideal if ramified, and list of 2 ideals if split.
    Assumes ell is coprime to conductor of O. Assumes discriminant of quadratic field containing O is squarefree.
    Note: Sage has 'factor' method which works on ideals or maximal orders, but doesn't support non-maximal, hence why we've implemented this ourselves.
    """
    f = O.conductor()
    z = O.number_field().gens()[0]
    N = ZZ(z**2)
    Zell = Integers(ell)
    if ell > 2:
        if Zell(N) == Zell(0): # ramified
            if mod(N, 4) == 1:
                return [O.ideal([ell, f*((ell + z)/2)])]
            return [O.ideal([ell, f*z])]
        if Zell(N).is_square(): # split
            x = ZZ(Zell(N).sqrt())
            if mod(x, 2) == 0:
                x = ell - x
            if mod(N, 4) == 1:
                return [O.ideal([ell, f*((x+z)/2)]), O.ideal([ell, f*((x-z)/2)])]
            return [O.ideal([ell, f*x + f*z]), O.ideal([ell, f*x - f*z])]
        return [] # inert
    else: # ell = 2
        if mod(N, 8) == 1: return [O.ideal([2, f*((1+z)/2)]), O.ideal([2, f*((1-z)/2)])]
        if mod(N, 2) == 0: return [O.ideal([2, f*z])]
        if mod(N, 4) == 3: return [O.ideal([2, 1+f*z])]
        return [] # N = 5 mod 8

def imaginary_quadratic_elt_to_string(elt):
    """Given an element of an imaginary quadratic field, returns a nicely formatted string representing that element."""
    outstr = ""
    z = parent(g).gens()[0]
    z1 = ZZ(z**2).squarefree_part()
    z2 = sqrt(z**2 / z1)
    firstpart = elt[0] != 0
    secondpart = z2*elt[1] != 0
    if elt[0] != 0:
        outstr += str(elt[0])
        if elt[1] > 0 and secondpart: outstr += " + "
    if not secondpart: return outstr
    if z2*elt[1] < 0: outstr += " - "
    if (z2*elt[1]).abs() == 1: return outstr + "sqrt(" + str(z1) + ")"
    return outstr + str((-1)**((z2*elt[1] < 0))*z2*elt[1]) + ".sqrt(" + str(z1) + ")"