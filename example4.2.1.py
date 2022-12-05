# the Sage code to do the construction for Example 4.2.1 of Marlon dissertation

p = 2
F = QpER( p, print_mode = 'digits', secure = False )
R = ZpER( p, print_mode = 'digits', secure = False )
G = AbelianGroup( [2,2] )
FG = GroupAlgebra( G, F )
RG = GroupAlgebra( G, R )

mats = [ matrix( F, [[-1,1,0],[0,1,0],[0,0,-1]] ), matrix( F, [[-1,1,1],[0,1,0],[0,0,1]])]






