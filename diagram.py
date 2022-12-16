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

def idempotents( FG ):
    """Determines the idempotents of the group algebra FG.
    At the moment G must be C_p x C_p and F the field of p-adic numbers!"""

    G = FG.group()
    p = ZZ(prime_divisors( G.order())[0])
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


    # compute the generators of the subspace V and for the subspaces Vi at the 
    # same time
    for k in range( nr_ids ):   
        for row in ims[k]:

            # row is the image of e_j under the k-th idempotent considered over 
            # the residue class ring R0
            row = vector( R0(x) for x in row )
            
            # if row is zero, nothing to do
            if row.is_zero():
                continue
            
            # we reduce row against the already found generators of Vi[k]
            row_Vi = row_reduce( matrix( gens_Vi[k]), row )[0]

            # if the reduction is zero, then nothing more to do
            if row_Vi.is_zero():
                continue
            else:
                # we insert the reduction of row into the generators of Vi[k]
                le = row_Vi.support()[0] # the leading entry of row

                # we normalize the leading entry
                q = row_Vi[le]*ZZ(p)**-valuation( row_Vi[le], p ) 
                row_Vi /= q**-1

                # find the positin of row_Vi in the matrix Vi[k]
                _, pos = search( leads_Vi[k], le )

                # insert row_Vi into the right position and insert its leading entry also
                leads_Vi[k].insert( pos, le )
                gens_Vi[k].insert( pos, row_Vi )
            
            # now we further reduce row_Vi agains the generators of V
            row_V = row_reduce( matrix( gens_V ), row_Vi )[0]

            if not row_V.is_zero():
                # we repeat the same computation as we did after the last if statement
                le = row_V.support()[0]
                q = row_V[le]*ZZ(p)**-valuation( row_V[le], p )
                row_V /= q**-1
                _, pos = search( leads_V, le )
                leads_V.insert( pos, le )
                gens_V.insert( pos, row_V )    

    # compute the matrices that define the G action on V
    mats_V = []
    for m in mats:
        mats_V.append( matrix( [ row_reduce( matrix( gens_V ), 
                            r*matrix( R0, m ), is_member = true )[1] for r in gens_V ] ))


    # compute the matrices that define the G-action on the Vi 
    mats_Vi = [ [] for _ in range( nr_ids )]
    for i in range( nr_ids ):
        for m in mats:
            mei = matrix( R0, m )
            mats_Vi[i].append( matrix( [ row_reduce( matrix( gens_Vi[i] ), 
                    r*mei, is_member = True  )[1] for r in gens_Vi[i] ] ))

    # Build the Butler diagram and return
    return ButlerDiagram( G, p, Q, Q.integer_ring(), R0, ids, mats, gens_V, gens_Vi, mats_V, mats_Vi )
                 
    








