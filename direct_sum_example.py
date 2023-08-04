def direct_sum_diagram( dim_list, p ):

    dim_sum = sum( dim_list )
    V = GF( p )**dim_sum 

    bas_V = V.basis()
    spaces = [ matrix( bas_V ) ]

    d_sum = 0 
    for d in dim_list:
        spaces.append( matrix( bas_V[d_sum:d_sum+d] ))
        d_sum += d

    R = pAdicRing( p, 10, print_mode = "digits" )
    F = pAdicField( 3, 10, print_mode = "digits" )
    G = AbelianGroup( [p,p] )
    ids = idempotents( GroupAlgebra( G, F ))

    act_V0 = [ identity_matrix( GF( p ), dim_sum ) for _ in range( 2 ) ]
    act_V = [act_V0] + [[ identity_matrix( GF( p ), d ) for _ in range( 2 ) ] for d in dim_list ]
    
    return ButlerDiagram( G, 3, R, F, GF(p), ids, spaces[0], spaces, act_V0, act_V )
