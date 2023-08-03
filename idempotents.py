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

def group_element_from_gap( g, el ):
    
    rep_el = el.ExtRepOfObj()
    if rep_el == []:
        return g.one()

    g_gens = g.gens()
    return prod( [ g.gens()[int(rep_el[2*i+1])-1]**int(rep_el[2*i+2]) for i in range( len( rep_el )//2 )])
    
def get_idempotents_from_gap( g ):

    p = prime_divisors( g.order())[0]

    gapA = gap.GroupRing( QQ, g )
    gapids = gap.CentralIdempotentsOfAlgebra( gapA )
    print( "Central idempotents computer by GAP" )

    F = pAdicField( p, 10, print_mode = "digits" )
    A = GroupAlgebra( g, F )

    ids = []
    for i in gapids:
        comp_list = i.CoefficientsAndMagmaElements()
        ids.append( sum( F(comp_list[2*i+2])*A(group_element_from_gap( g, comp_list[2*i+1])) 
                for i in range( len( comp_list )//2 )))

    return ids 
