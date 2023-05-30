prec = 2
G = AbelianGroup( [2,2,2] )
FF = QpER( 2, prec,  print_mode = 'digits', secure = False )

mats = [ 
    matrix( FF, [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1],
        [0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1],
        [0, -1, 0, 0, 0, 0, 1, -1, -1, -1, -1],
        [0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1]] ),
    matrix( FF, [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1],
        [0, 1, -1, 0, 1, 0, 0, 0, 0, 1, 1],
        [1, 0, 1, -1, -1, -1, -1, 0, 0, -1, -1],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1]] ),
    matrix( FF,[[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 1, -1, 0, -1, 0],
        [1, 1, 0, 0, -1, -1, -1, -1, 0, -1, 0],
        [0, 0, 0, -1, 1, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1]])
]