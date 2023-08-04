def image( el, mats ):
    """el is an element of a group algebra FG and mats is a list of square matrices matrices over F. 
    Substitutes the elements in mats in the place of the generators of G.
    
    Can be used to calculate a homomorphism FG -> M(n,F) for fields where the standard homomorphism 
    machinery of Sage fails."""

    #initialize with the zero matrix
    im = 0*mats[0]
    R_mats = mats[0][0,0].parent()
    
    #take the dictionary with the terms and their coefficients
    terms = el.monomial_coefficients()
    
    for mon, coeff in terms.items():
        exps = mon.exponents()
        pr = prod( mats[m]**exps[m] for m in range(len(exps)))
        im += R_mats(coeff)*pr

    return im 


def Lambda_i_basis( idem ):
    
    R = idem.coefficients()[0].parent()
    p = R.prime()
    G = idem.parent().group()
    Ggens = G.gens()
    x, y = Ggens
    if x*idem == idem and y*idem == idem:
        return [ matrix( R, 1, 1, [1] ), matrix( R, 1, 1, [1] )], [G.one()]
    elif x*idem != idem:
        x0 = x
    else:
        x0 = y
    
    mat, els, g_els = [], [], []
    for k in range( p-1 ):
        el = x0**k*idem 
        mat.append( el.coefficients())
        els.append( el )
        g_els.append( x0**k )

    mat = Matrix( mat )
    mats = [[] for _ in range( len( Ggens ))]

    for i in range( len( Ggens )):
        g = Ggens[i]
        for el in els:
            coeffs = mat.solve_left( vector((g*el).coefficients()))
            mats[i].append( coeffs )

    return  [ matrix( mat ) for mat in mats ], g_els
    
def relations_eiZpG( ei ):

    A = ei.parent()
    G = A.group()
    p = ZZ(prime_divisors( G.order())[0])
    F = A.base_ring()
    mat = []

    for g in G:
        pr = g*ei
        moncoef = pr.monomial_coefficients()
        vec = vector( [ moncoef[g] for g in G ])
        mat.append( vec )

    mat = matrix( mat, sparse = False )
    mat_kernel = mat.left_kernel()
    return matrix( mat_kernel.basis())

def check_lifting_action_module( ei, mats ):
    
    A = ei.parent()
    G = A.group()
    p = ZZ(prime_divisors( G.order())[0])
    F = A.base_ring()
    mat = relations_eiZpG( ei )
    els = [ g for g in G ]
    
    for row in mat:
        r = sum( row[k]*A(els[k]) for k in range( G.order()))
        im = image( r, mats )
        if im != 0*im:
            return False
        
    return True

def check_lifting_action_diagram( d ):
    ids = d.idempotents
    mods = d.action_Vi

    for i in range( len( ids )):
        if not check_lifting_action_module( ids[i], mods[i] ):
            return False

    return True


# calculate the Jacobson radical
# also works for C_p x C_p only

def jacobson_radical( mats ):

    F = mats[0].parent().base_ring() # the underlying field
    # At the moment F must be a finite field. 
    assert F.is_field() and F.is_finite()

    if mats[0].dimensions() == (0,0):
        return matrix( F, [] )

    l = mats[0].dimensions()[0] # the dimension of the module
    V = F**l # the module

    jac_gens = [] # this will store the generators

    # the jacobson radical is generated by b - b*m where
    # b is an element of the basis of V and m is in mats 
    # THIS ONLY WORKS IF F IS A FIELD WITH p ELEMENTS!!!
    for m in mats:
        for b in V.basis():
            jac_gens.append( b-b*m )
    
    # We calculate the row-reduced form of the matrix whose rows are jac_gens
    jac_gens = matrix( jac_gens ).rref()
    # we delete the zero rows and return 
    return jac_gens[0:len(jac_gens.pivots())]


# V is a matrix generating a subspace which contains vec
# mat is the matrix of a group element acting on V in the basis given by mat
# returns the image of vec under mat
