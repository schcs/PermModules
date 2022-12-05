def depth( el ):
    return min( valuation(c) for c in el )

def depth_list( ellist ):
    return min( depth( el ) for el in ellist )

def depth_matrix( mat ):
    r, c = mat.nrows(), mat.ncols()
    return depth( vector([ mat[i,j] for i in range(r) for j in range( c )])) 

def image( el, mats ):

    im = 0*mats[0]

    terms = el.monomial_coefficients()
    for k in terms.keys():
        exps = k.exponents()
        pr = prod( mats[m]**exps[m] for m in range(len(exps)))
        im += terms[k]*pr

    return im 

def idempotents( FG ):
    """Determines the idempotents of the group algebra FG.
    At the moment G must be C_p x C_p and F the field of p-adic numbers."""

    G = FG.group()
    p = ZZ(prime_divisors( G.order())[0])
    n, c = G.gens()
    R = FG.base_ring()

    g_hat = sum( [ FG(x) for x in G ])

    idems = [ p**-2*g_hat ]
    idems.append( (p**-2)*(p*sum([ FG(n**k) for k in range(p)]) - g_hat ))

    for i in range(p): 
        idems.append( p**-2*(p*sum([ FG((c*n**i)**k) for k in range( p )]) - g_hat ))
    
    return idems

def e_iU_plus_U_over_U( e, mats ):

    R = e.parent().base_ring()
    im = image( e, mats )

    if im.is_zero():
        dt = 0
    else: 
        dt = depth_matrix( im ) 

    im *= ZZ(p)**-dt

    r = im.ncols()
    for i in range( r ):
        v = zero_vector( R, r )
        v[i] = ZZ(p)**-dt
        im = im.stack( v )
    
    R1 = IntegerModRing( ZZ(p)**-dt )
    #return im
    return matrix( R1, hnf( im )[0] )


def diagram( G, mats ):
    """Constructs the Butler diagram for the Z_pG-module U. G must be a p-group. The action of 
    the generators of G on U is determined by the matrices in mats."""

    p = ZZ(prime_divisors( G.order())[0])
    R, Q = Zp( p ), Qp( p )
    FG = GroupAlgebra( G, F )
    ids = idempotents( FG )

    return [ e_iU_plus_U_over_U( e, mats ) for e in ids ]

    








