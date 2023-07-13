class ButlerDiagram:
    """The data structure to store the information related to a Butler diagram.
    If d is a Butler diagral, then d has the following components:
    d.group:                    the acting p-group G
    d.p:                        the prime p
    d.padic_ring:               the underlying p-adic ring
    self.padic_field:           the underlying p-adic field F
    self.residue_class_ring:    the residue class ring for which the V and Vi are modules
    self.idempotents:           the idempotents of FG
#    self.action:                the action of G on the original module U  TO BE DELETED
    self.V:                     the space V in the diagram
    self.Vi:                    the list of spaces Vi in the diagram
    self.action_V:              the matrices defining the action of G on V
    self.action_Vi:             the sequence of matrices defining the action of G on the Vi

    """
    def __init__( self, group, p, R, F, R0, ids, V, Vi, V_gens, Vi_gens ):
        self.group = group
        self.p = p
        self.padic_ring = R
        self.padic_field = F
        self.residue_class_ring = R0
        self.idempotents = ids
#        self.action = mats
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
    R_mats = mats[0][0,0].parent()
    
    #take the dictionary with the terms and their coefficients
    terms = el.monomial_coefficients()
    
    for mon, coeff in terms.items():
        exps = mon.exponents()
        pr = prod( mats[m]**exps[m] for m in range(len(exps)))
        im += R_mats(coeff)*pr

    return im 

def klein_subgroup( a, b ):
    return [a**0, a, b, a*b ]

def idempotents_222( FG ):
    
    G = FG.group()
    p = ZZ(prime_divisors( G.order())[0])
    n, b, c = G.gens()
    
    # the sum of the elements of G

    g_hat = sum( [ FG(x) for x in G ])
    print( g_hat )
    idems = [ ZZ(p)**-3*g_hat ]
    max = [ klein_subgroup( n, b ), klein_subgroup( n, c ), klein_subgroup( b, c ),
            klein_subgroup( n*b, c ), klein_subgroup( n*c, b ), klein_subgroup( b*c, n ), 
            klein_subgroup( n*b, n*c )]

    for m in max:
        idems.append( ZZ(p)**-3*( 2*sum( FG( x ) for x in m ) - g_hat ))
    
    return idems

# calculate the support of an idempotent
# idempotents are of the form 1/p^2(sum(G)) or 1/p^2(sum(H)-sum(G))
# In the first case, it returns nothing, in the second case it returns a generator of H
# Only works for C_p x C_p
def idempotent_subgroup( idem ):
    
    p = idem.coefficients()[0].parent().prime()    
    mon_coeff = idem.monomial_coefficients()
    for k in mon_coeff.keys():
        if mon_coeff[k] == 1/p-1/p**2 and k != k**0:
            return k
    
    return False
 
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
    k = depth_matrix( mat )
    #return mat, k
    return hnf((p**-k)*mat, normalize = true )[0]

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

def butler_diagram( lat ):
    """Constructs the Butler diagram for the Z_pG-module U. G must be a p-group. The action of 
    the generators of G on U is determined by the matrices in mats."""

    # first collect data from lat
    G, mats = lat.group, lat.mat_gens
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
                vec = vector( R0(x) for x in p**z*row ); 
                if vec.is_zero():
                    continue
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
                IntegerModRing( p**-depth ), ids, gens_V, gens_Vi, mats_V, mats_Vi )