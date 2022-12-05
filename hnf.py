def hnf( mat0, normalize = True ):
    
    """Hermite normal form for matrices over p-adic integers
    Takes matrix mat0 as input. It must be a matrix with entries of p-adic integers. 
    Returns the Hermite normal form m0 of mat0 together with the transformation matrix
    trans such that trans*mat0 = m0.
    
    If the optional argument normalize is set true (default), then the rows of m0 are normalized 
    in the sense that the leading entry of every row will be 1*p^k+... with some k>=0.""" 

    mat = copy( mat0 )
    m, n = mat.nrows(), mat.ncols()

    R = parent( mat[0,0] )
    i, j = 0, 0
    p = R.prime()
    trans = identity_matrix( R, m )

    while i < m and j < n:
        if not False in [ mat[k,j] == 0 for k in range(i,m) ]:
            j += 1
        else:
            while True:
                stop_flag = true  
                for k in range( i, m ):
                    for l in range( k+1, m ):
                        if mat[k,j] != 0 and mat[l,j] != 0:
                            stop_flag = false
                            if valuation( mat[k,j]) <= valuation( mat[l,j] ):
                                q = mat[l,j]/mat[k,j]
                                mat.add_multiple_of_row( l, k, -q )
                                trans.add_multiple_of_row( l, k, -q )
                                #print( "add", -q,  "times row ", k, "to row", l )
                            else: 
                                q = mat[k,j]/mat[l,j]
                                mat.add_multiple_of_row( k, l, -q )
                                trans.add_multiple_of_row( k, l, -q )
                                #print( "add", -q,  "times row ", l, "to row", k )

                if stop_flag:
                    break 
            
            for k in range(i,m):
                if mat[k,j] != 0:
                    break 
            
            if k != i:
                mat.swap_rows( k, i )
                trans.swap_rows( k, i )
                #print( "swap rows ", k, " and ", i )
            
            #print( "before norm", mat[i,j] )
            if normalize:
                q = mat[i,j]*(p**-valuation( mat[i,j] ))
                mat.rescale_row( i, q**-1 )
                trans.rescale_row( i, q**-1 )
                #print( "rescale row ", i, "by ", q**-1 )

            
            #print( "after norm", mat[i,j] )

            for l in range(i):
                if valuation( mat[l,j] ) >= valuation( mat[i,j] ):
                    q = mat[l,j]/mat[i,j]
                    mat.add_multiple_of_row( l, i, -q )
                    trans.add_multiple_of_row( l, i, -q )
                    #print( "add", -q,  "times row ", i, "to row", l )

            
            i += 1; j += 1
    return mat, trans