# the Sage code to do the construction for Example 5.1

p = 3
F = QpER( p, print_mode = 'digits', secure = False )
R = ZpER( p, print_mode = 'digits', secure = False )
G = AbelianGroup( [3,3] )
FG = GroupAlgebra( G, F )
RG = GroupAlgebra( G, R )

mats = [ Matrix( F, [[0,1,0],[0,0,1],[1,0,0]] ), Matrix( F, [[1,0,0],[0,1,0],[0,0,1]] )]






