def klein_subgroup( a, b ):
    return [a**0, a, b, a*b ]

def idempotents_222( FG ):
    
    G = FG.group()
    p = ZZ(prime_divisors( G.order())[0])
    n, b, c = G.gens()
    
    # the sum of the elements of G

    g_hat = sum( [ FG(x) for x in G ])
    print( g_hat )
    idems = [ ZZ(p)**-3*g_hat ]
    max = [ klein_subgroup( n, b ), klein_subgroup( n, c ), klein_subgroup( b, c ),
            klein_subgroup( n*b, c ), klein_subgroup( n*c, b ), klein_subgroup( b*c, n ), 
            klein_subgroup( n*b, n*c )]

    for m in max:
        idems.append( ZZ(p)**-3*( 2*sum( FG( x ) for x in m ) - g_hat ))
    
    return idems

# calculate the support of an idempotent
# idempotents are of the form 1/p^2(sum(G)) or 1/p^2(sum(H)-sum(G))
# In the first case, it returns nothing, in the second case it returns a generator of H
# Only works for C_p x C_p
def idempotent_subgroup( idem ):
    
    p = idem.coefficients()[0].parent().prime()    
    mon_coeff = idem.monomial_coefficients()
    for k in mon_coeff.keys():
        if mon_coeff[k] == 1/p-1/p**2 and k != k**0:
            return k
    
    return False

def idempotents( FG ):
    """Determines the idempotents of the group algebra FG.
    At the moment G must be C_p x C_p and F the field of p-adic numbers!"""


    G = FG.group()
    p = ZZ(prime_divisors( G.order())[0])
    
    if G.elementary_divisors() == (2,2,2):
        return idempotents_222( FG )

    n, c = G.gens()

    # the sum of the elements of G
    g_hat = sum( [ FG(x) for x in G ])

    # the idempotent e_0
    idems = [ p**-2*g_hat ]

    # the idempotent corresponding to the subgroup <n>
    idems.append( (p**-2)*(p*sum([ FG(n**k) for k in range(p)]) - g_hat ))

    for i in range(p): 
        # the idempotent corresponding to the subgroup <cn^i>
        idems.append( p**-2*(p*sum([ FG((c*n**i)**k) for k in range( p )]) - g_hat ))
    
    return idems

def right_regular_action( g ):

    g_elts = g.list()
    g_gens = g.gens()
    perms = []
    sym_gp = SymmetricGroup( g.order()) 
    for g in g_gens:
        perm = sym_gp( [ g_elts.index( x*g ) + 1 for x in g_elts ])
        perms.append( perm )
    
    return sym_gp.subgroup( perms, canonicalize = False )


def subgroups_with_cyclic_quotient( g ):
    gp = right_regular_action( g )

    subs_gp = [ x for x in gp.subgroups() if gp.quotient(x).is_cyclic() ]
    no_gens = len( g.gens())
    gen_dict = { str(gp.gens()[k]): g.gens()[k] for k in range( no_gens )}
    #return gen_dict
    
    subs = []
    for h in subs_gp:
        gens_h_perm = h.gens()
        gens_h = []
        for x in gens_h_perm: 
            if x == x**0:
                continue
            w = x.word_problem( gp.gens(), as_list = True, display = False )
            gens_h.append( prod( gen_dict[x[0]]**x[1] for x in w ))
        subs.append( g.subgroup( gens_h ))

    return subs 

# checks if m/n is minimal normal in g/n
# g >= m >= n

def is_min_normal_subgroup( g, n, m, p ):

    return n.is_subgroup( m ) and (( g.order()//n.order() > p and m.order()//n.order() == p ) or 
                                   ( g.order()//n.order() == p and m.order() == g.order() ))
    
def idempotents_of_group( g, F = False ):

    subs = subgroups_with_cyclic_quotient( g )
    p = prime_divisors( g.order())[0]
    
    if type( F ) == bool:
        F = pAdicField( p, 10, print_mode = "digits" )
    
    FG = GroupAlgebra( g, F )    
    H_hat = lambda H: H.order()**-1*sum( FG( x ) for x in g if x in H )
    ids = [ H_hat( g )]

    for n in subs:
        idem = FG.one()
        to_append = False

        for m in subs:
            if is_min_normal_subgroup( g, n, m, p ):
                #print( m, n )
                idem *= H_hat( n ) - H_hat( m )
                to_append = True 
        
        if to_append: 
            ids.append( idem )

    return ids

    

    return subgroups_cyclic_quot

def group_element_from_gap( g, el ):
    
    rep_el = el.ExtRepOfObj()
    if rep_el == []:
        return g.one()

    g_gens = g.gens()
    return prod( [ g.gens()[int(rep_el[2*i+1])-1]**int(rep_el[2*i+2]) for i in range( len( rep_el )//2 )])
