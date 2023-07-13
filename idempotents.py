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
