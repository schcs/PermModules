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