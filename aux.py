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


# compute the complement of U in W
def complement( W, U ):
    
    # check if subspace 
    assert U.is_subspace( W )

    # calculate basis for U
    vectsU = [ u for u in U.basis()] 

    # this will hold the generating set of the complement   
    vects = []
    for v in W.basis():
        if not (v in U):        
            vects.append( v )
            vectsU.append( v )
            U = W.subspace( vectsU )

    return W.span( vects )


def image_vector_under_action( V, mat, vec ):

    # get the coefficients of vec in the linear combination 
    # of the rows of mat
    return V.solve_left( vec )*mat*V 
