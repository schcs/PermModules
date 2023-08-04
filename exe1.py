G=AbelianGroup([3,3])
k=1
r=5
V_i_gens=[[[1,0,0,0,0],[0,0,1,1,1],[0,-1,-1,0,1]],[[1,0,0,0,0],[0,1,0,0,0]], [[1,0,0,0,0],[0,0,1,0,0]],[[1,0,0,0,0],[0,0,0,1,0]],[[1,0,0,0,0],[0,0,0,0,1]]]
mats=[[matrix(IntegerModRing(3),[[1,0,0],[0,1,0],[0,0,1]]),matrix(IntegerModRing(3),[[1,0,0],[0,1,0],[0,0,1]])], 
       [matrix(IntegerModRing(3),[[1,0],[0,1]]),matrix(IntegerModRing(3),[[1,0],[1,1]])],
       [matrix(IntegerModRing(3),[[1,0],[1,1]]),matrix(IntegerModRing(3),[[1,0],[0,1]])],
       [matrix(IntegerModRing(3),[[1,0],[1,1]]),matrix(IntegerModRing(3),[[1,0],[2,1]])],
       [matrix(IntegerModRing(3),[[1,0],[1,1]]),matrix(IntegerModRing(3),[[1,0],[1,1]])]]


mat_lambda_i=[[matrix(ZZ,[[1]]),matrix(ZZ,[[1,0],[0,1]]), matrix(ZZ,[[0,-1],[1,-1]]), matrix(ZZ,[[0,-1],[1,-1]]), matrix(ZZ,[[0,-1],[1,-1]])],[matrix(ZZ,[[1]]),matrix(ZZ,[[0,-1],[1,-1]]), matrix(ZZ,[[1,0],[0,1]]), matrix(ZZ,[[-1,1],[-1,0]]), matrix(ZZ,[[0,-1],[1,-1]])]]

V_gens = identity_matrix( 5 )
mats_V = [ matrix( GF(3), [[1,0,1,1,1], [0,1,0,0,0], [0,0,1,0,0], [ 0,0,0,1,0], [0,0,0,0,1]]).transpose(), 
           matrix( GF(3), [[1,1,0,2,1],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]]).transpose() ]

R = pAdicRing( 3, 10, print_mode = "digits" )
F = pAdicField( 3, 10, print_mode = "digits" )
ids = idempotents( GroupAlgebra( G, F ))
diag = ButlerDiagram( G, 3, R, F, GF(3), ids, matrix( V_gens ), 
                            [ matrix( x ) for x in V_i_gens ], mats_V, mats )
