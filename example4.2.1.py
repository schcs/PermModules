# the Sage code to do the construction for Example 4.2.1 of Marlon dissertation

prec = 3
FF = QpER( 2, prec,  print_mode = 'digits', secure = False )
G = AbelianGroup( [2,2] )

mats = [ matrix( FF, [[-1,1,0],[0,1,0],[0,0,-1]] ), 
         matrix( FF, [[-1,1,1],[0,1,0],[0,0,1]])]






