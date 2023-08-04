G=AbelianGroup([2,2])
k=1
r=4
V_i_gens=[[[0,0]],[[1,0]],[[1,0],[0,1]],[[0,1]]]
mats=[[matrix(IntegerModRing(2),[[1]]),matrix(IntegerModRing(2),[[1]])],[matrix(IntegerModRing(2),[[-1]]),matrix(IntegerModRing(2),[[1]])],[matrix(IntegerModRing(2),[[1,0],[0,1]]),matrix(IntegerModRing(2),[[-1,0],[0,-1]])],[matrix(IntegerModRing(2),[[-1]]),matrix(IntegerModRing(2),[[-1]])]]
mat_lambda_i=[[matrix(ZZ,[[1]]),matrix(ZZ,[[1]]),matrix(ZZ,[[-1]]),matrix(ZZ,[[-1]])], [matrix(ZZ,[[1]]),matrix(ZZ,[[-1]]),matrix(ZZ,[[1]]),matrix(ZZ,[[-1]])]]


V_gens = identity_matrix( 2 )
mats_V = [ matrix( GF(3), [[1,0],[0,1]] ), 
           matrix( GF(3), [[1,0],[0,1]]) ]


R = pAdicRing( 2, 10, print_mode = "digits" )
F = pAdicField( 2, 10, print_mode = "digits" )
ids = idempotents( GroupAlgebra( G, F ))
diag = ButlerDiagram( G, 2, R, F, GF(2), ids, matrix( GF(2), V_gens ), 
                            [ matrix( GF(2), x ) for x in V_i_gens ], mats_V, mats )
