# In this file we collect functions that calculate a lattice from a Butler diagram.

class padicLattice:
    """The data structure to store information related to a GZ_p-lattice."""
    def __init__( self, G, mat_gens ):
        self.group = G                                 # the acting group
        self.mat_gens = mat_gens                       # the matrix generator
                                                       # that correspond to the generators of G
        self.p = prime_divisors( G.order())[0]         # the prime p
        self.padic_ring = mat_gens[0][0][0].parent()   # the p-adic ring

    def __repr__( self ):
        return f'{self.p}-adic lattice for {self.group}'

    
def lattice_from_diagram( diag ):

    # get some data from diag
    G, F, R = diag.group, diag.padic_field, diag.padic_ring
    p = F.prime()
    rankG = len( diag.action_V )

    # define the group algebra
    FG = GroupAlgebra( G, F )
    # calculate the idempotents
    ids = diag.idempotents
    no_ids = len( ids )
    r = len( diag.Vi )
    dimV = diag.V.dimensions()[1] 
    V = GF(p)**dimV 
    
    # calculating the G-action on the Lambda_i
    lambdas = [ Lambda_i_basis( i ) for i in ids ]
    # we calculate the Jacobson radicals for the modules in the diagram
    jacobson_radicals = [ jacobson_radical( mats ) for mats in diag.action_Vi ]
    
    # this will hold the multiplicities of the Lambda_i
    lambda_multiplicities = []

    # starting to build the list of images
    images = []

    for i in range( r ):
        Vi = V.span( diag.Vi[i] )
        #return jacobson_radicals[i], diag.Vi[i]
        gens_rad = jacobson_radicals[i]*diag.Vi[i]
        rad = V.span( gens_rad )
        comp = complement( Vi, rad )
        ri = comp.dimension()
        lambda_multiplicities.append( ri )
        # the list of elements that we need to define the images of the generators of the 
        # Lambda_i

        acting_elements =  lambdas[i][1]
        acting_elements_mat = [ image( FG( x ), diag.action_Vi[i] ) for x in acting_elements ]

        # this will hold the images of the map Lambda_i^r_i -> V
        images_i = []
    
        for j in range( ri ):
            for m in acting_elements_mat:
                images_i.append( image_vector_under_action( diag.Vi[i], m, comp.basis()[j]))

        images += images_i 
    
    images = matrix( GF( p ), images )
    ims_ker = matrix( ZZ, [ x for x in images.left_kernel().basis()])
    ims_ker = [ x for x in matrix( R, ims_ker )]
    
    W = R**len( ims_ker[0] )
    for b in W.basis():
        ims_ker.append( p*b )
    
    mat = hnf( matrix( R, ims_ker ))[0]
    
    mat = matrix( x for x in mat if not x.is_zero())

    dim_U = mat.dimensions()[0]
    gens_G_on_V  = []
    for i in range( rankG ):
        blocks = sum( [[ lambdas[j][0][i] for _ in range( lambda_multiplicities[j] )] for j in range( r )], [] )
        gens_G_on_V.append( block_diagonal_matrix( blocks ))
    
    gens_G_on_U = []
    for i in range( rankG ):
        gens_G_on_U.append( matrix( [ row_reduce( mat, mat[j]*gens_G_on_V[i] )[1] for j in range( dim_U )]))

    return padicLattice( G, gens_G_on_U )

