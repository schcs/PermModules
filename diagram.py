class ButlerDiagram:
    """The data structure to store the information related to a Butler diagram.
    If d is a Butler diagral, then d has the following components:
    d.group:                    the acting p-group G
    d.p:                        the prime p
    d.padic_ring:               the underlying p-adic ring
    self.padic_field:           the underlying p-adic field F
    self.residue_class_ring:    the residue class ring for which the V and Vi are modules
    self.idempotents:           the idempotents of FG
    self.action:                the action of G on the original module U
    self.V:                     the space V in the diagram
    self.Vi:                    the list of spaces Vi in the diagram
    self.action_V:              the matrices defining the action of G on V
    self.action_Vi:             the sequence of matrices defining the action of G on the Vi

    """
    def __init__( self, group, p, R, F, R0, ids, mats, V, Vi, V_gens, Vi_gens ):
        self.group = group
        self.p = p
        self.padic_ring = R
        self.padic_field = F
        self.residue_class_ring = R0
        self.idempotents = ids
        self.action = mats
        self.V = V
        self.Vi = Vi
        self.action_V = V_gens
        self.action_Vi = Vi_gens 

    # the function  that defines how a Butler diagram is printed
    def __repr__( self ):
        return f'Butler diagram for the group {self.group} over {self.padic_ring}'

def leading_position( vec ):
    return vec.support()[0]
    
def depth( el ):
    """The depth of a vector with entries in a p-adic field or integer ring. 
       Returns k where p^k is the minimum valuation among the entries in el."""
    return min( valuation(c) for c in el )

def depth_list( ellist ):
    """Returns the minimum among the depths of ellinst."""
    return min( depth( el ) for el in ellist )

def depth_matrix( mat ):
    """mat is a matrix with entries in a p-adic field or integer ring. Returns the 
    minimum among the valuations of the entries of mat."""
    r, c = mat.nrows(), mat.ncols()
    return depth( vector([ mat[i,j] for i in range(r) for j in range( c )])) 

def depth_matrix_list( mlist ):
    """return the minimum among the depths of the matrices in mlist"""
    return min( depth_matrix( mat ) for mat in mlist )

def image( el, mats ):
    """el is an element of a group algebra FG and mats is a list of square matrices matrices over F. 
    Substitutes the elements in mats in the place of the generators of G.
    
    Can be used to calculate a homomorphism FG -> M(n,F) for fields where the standard homomorphism 
    machinery of Sage fails."""

    #initialize with the zero matrix
    im = 0*mats[0]
    
    #take the dictionary with the terms and their coefficients
    terms = el.monomial_coefficients()

    for mon, coeff in terms.items():
        exps = mon.exponents()
        pr = prod( mats[m]**exps[m] for m in range(len(exps)))
        im += coeff*pr

    return im 

def klein_subgroup( a, b ):
    return [a**0, a, b, a*b ]

def idempotents_222( FG ):
    
    G = FG.group()
    p = ZZ(prime_divisors( G.order())[0])
    n, b, c = G.gens()
    
    # the sum of the elements of G

    g_hat = sum( [ FG(x) for x in G ])
    idems = [ ZZ(p)**-3*g_hat ]
    max = [ klein_subgroup( n, b ), klein_subgroup( n, c ), klein_subgroup( b, c ),
            klein_subgroup( n*b, c ), klein_subgroup( n*c, b ), klein_subgroup( b*c, n ), 
            klein_subgroup( n*b, n*c )]

    for m in max:
        idems.append( ZZ(p)**-3*( 2*sum( FG( x ) for x in m ) - g_hat ))
    
    return idems


def idempotents( FG ):
    """Determines the idempotents of the group algebra FG.
    At the moment G must be C_p x C_p and F the field of p-adic numbers!"""


    G = FG.group()
    p = ZZ(prime_divisors( G.order())[0])
    
    if G.elementary_divisors() == (2,2,2):
        return idempotents_222( FG )

    n, c = G.gens()

    # the sum of the elements of G
    g_hat = sum( [ FG(x) for x in G ])

    # the idempotent e_0
    idems = [ p**-2*g_hat ]

    # the idempotent corresponding to the subgroup <n>
    idems.append( (p**-2)*(p*sum([ FG(n**k) for k in range(p)]) - g_hat ))

    for i in range(p): 
        # the idempotent corresponding to the subgroup <cn^i>
        idems.append( p**-2*(p*sum([ FG((c*n**i)**k) for k in range( p )]) - g_hat ))
    
    return idems

def e_iU_plus_U_over_U( e, mats, depth = -1, V = [] ):
    """Suppose U is an FG-module for a finite p-group G and F a p-adic field. 
    The action of the generators of G is given by the matrices in mats. 
    For an idempontent u, the function calculates the space (e_iU+U)/U as a subspace of 
    the free Z/p^kZ-module where k=depth."""

    #OBSOLATE WILL BE REMOVED

    # get some basic data about FG
    R = e.parent().base_ring()
    p = R.prime()
    R0 = IntegerModRing( p**-depth )

    # calculate the image of e over R0
    im = matrix( R0, ZZ(p)**-depth*image( e, mats ))

    # calculate HNF for standard basis
    im = matrix( R0, hnf( im )[0] )

    # remove zeros
    return im[[x for x in range(im.nrows()) if not im[x].is_zero() ]]

#
def diagram_V( G, mats ):
    """Suppose U is an FG-module for a finite p-group G and F a p-adic field. 
    The action of the generators of G is given by the matrices in mats. 
    Constructs the space V = (sum(V_i) + U)/U where the V_i are the subspaces 
    (e_iU+U)/U for the idempontents e_i of FG."""

    # OBSOLATE WILL BE REMOVED


    p = ZZ(prime_divisors( G.order())[0])
    F = mats[0].base_ring()
    FG = GroupAlgebra( G, F )
    ids = idempotents( FG )
    ims = [ image( id, mats ) for id in ids ]
    dt = depth_matrix_list( ims )
    R = IntegerModRing( p**-dt )    

    gens_V = []
    leads = []

    for im in ims:
        im0 = matrix(R, ZZ(p)**-dt*im) 

        for r in im0:
            if r.is_zero():
                continue
            rr = row_reduce( matrix( R, gens_V ), r )[0]
            #return rr
            if not rr.is_zero():
                le = rr.support()[0]
                q = rr[le]*ZZ(p)**-valuation( rr[le], p )
                rr /= q**-1
                _, pos = search( leads, le )
                leads.insert( pos, le )
                gens_V.insert( pos, rr ) 
    
    return gens_V
    return matrix( R, gens_V )

def butler_diagram( G, mats ):
    """Constructs the Butler diagram for the Z_pG-module U. G must be a p-group. The action of 
    the generators of G on U is determined by the matrices in mats."""

    # first collect data from G and mats

    p = ZZ(prime_divisors( G.order())[0])
    Q = mats[0][0,0].parent()
    FG = GroupAlgebra( G, Q )

    # calculate the idempotents 
    ids = idempotents( FG )
    nr_ids = len( ids )

    # calculate the images of the idempontents under the representation
    ims = [ image( e, mats ) for e in ids ]

    # we determine the largest k such that p^-k appears in the matrix entries 
    # -k is the depth
    depth = depth_matrix_list( ims )

    # we multiply with p**-depth so that we get matrices over Z_p
    ims = [ p**-depth*im for im in ims ]

    # Z/p**-depth Z will be the residue class ring over which the spaces V and Vi are defined
    R0 = IntegerModRing( ZZ(p)**-depth )

    # se set up the lists for the generators of V and the Vi
    gens_V = []
    gens_Vi = [ [] for _ in range( nr_ids )]

    # leads_V will contain the position of the leading entries of V
    leads_V = []

    # the same for each of the Vi
    leads_Vi = [[] for _ in range( nr_ids )]


    # Compute the generators for the subspaces Vi
    for k in range( nr_ids ):   
        for row in ims[k]:
            for z in range( -depth ):
                vec = vector( R0(x) for x in p**z*row )
                if vec.is_zero():
                    break
                gens_Vi[k].append( vec )
        
        gens_Vi[k] = hnf( matrix( gens_Vi[k] ))[0]
        gens_Vi[k] = matrix( [ x for x in gens_Vi[k] if not x.is_zero()])

    for i in range( nr_ids ):
        for r in gens_Vi[i]:
            gens_V.append( r )

    gens_V = hnf( matrix( gens_V ))[0]
    gens_V = matrix( [ x for x in gens_V if not x.is_zero()])

    # compute the matrices that define the G action on V
    mats_V = []
    for m in mats:
        mei = matrix( R0, m )
        mats_V.append( matrix( [ row_reduce( matrix( gens_V ), 
                            r*mei, is_member = True )[1] for r in gens_V ] ))


    # compute the matrices that define the G-action on the Vi 
    mats_Vi = [ [] for _ in range( nr_ids )]
    for i in range( nr_ids ):
        for m in mats:
            mei = matrix( R0, m )
            mats_Vi[i].append( matrix( [ row_reduce( matrix( gens_Vi[i] ), 
                    r*mei, is_member = True  )[1] for r in gens_Vi[i] ] ))
    
    return ButlerDiagram( G, p, Q, Q.integer_ring(), 
                IntegerModRing( p**-depth ), ids, mats, gens_V, gens_Vi, mats_V, mats_Vi )