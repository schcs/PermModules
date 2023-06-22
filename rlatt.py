## Código escrito por Marlon
## Modificado por Csaba


class Reticulado2:
    """The data structure to store information related to a GZ_p-lattice."""
    __init__( self, G, mat_gens ):

        self.group = G
        self.mat_gens = mat_gens 
        self.p = prime_divisors( G.order())[0]
        self.padic_ring = mat_gens[0][0][0].parent()
    
# mat is matrix
# gen_U is generating set of Z-module. 
# U invariant under mat
# returns the matrix of the action of mat on U?

def act_U(mat,gen_U):
    l=[mat*u for u in gen_U] 

    m = transpose(matrix(QQ,gen_U))
    
    coor= len(l)*[] 
    for j in range(len(l)):
        coor[j]=m.solve_right(vector(QQ,l[j]))
    mact = transpose(matrix(QQ,coor))
    return mact

# p é primo
# k >= 1
# gen1 é lista de vetores com entradas em Z/(p^k) de comprimento l
# gen2 é lista de vetores com entradas em Z/(p^k) de comprimento l
# os elementos de gen1 pertencem ao módulo gerado por gen2
# mat é matriz quadrada m x m onde m == len( gen2 ) com entradas em Z/(p^k)
#
# output: devolve uma lista das imagens dos vetores em gen1 pelo morfismo V -> V determinado pela matriz mat 

def mudança(p,k,gen1,gen2,mat):# gen1 é lista de elementos de um módulo que desejamos conhecer a imagem por mat e gen2 os geradores desse módulo que ja conhecemos a imagem por mat
    l=len(vector(gen1[0]))
    F=ZZ**l
    R=IntegerModRing(p**k)
    mgen=transpose(matrix(R,[x for x in gen2]))

    m=[[] for _ in range(len(gen1))]
    for j in range(len(gen1)):
        m[j]=mgen.solve_right(vector(R,gen1[j])) 

    coorde=[[] for _ in range(len(gen1))]
    for j in range(len(gen1)):
        coorde[j]=mat*m[j]

    im=[[] for _ in range(len(gen1))]
    for j in range(len(gen1)):
        im[j]=sum(coorde[j][i]*F(gen2[i]) for i in range(len(gen2)))
    return im

# create the lattice from the butler diagram

def lattice( diag ):
    G = diag.G
    mats = diag.action_Vi
    
    k,r,mats, V_i_gens, mat_lambda_i
    # primeiro calculamos JVi e os quocientes Vi/JVi

    p=ZZ(prime_divisors(G.order())[0])
    
    R=Zp(p)

    l=len(vector(V_i_gens[0][0]))

    F0=ZZ**l

    F1=F0.submodule(p**k*gens(F0)[j] for j in range(l))

    F2=F0/F1

    Vi=[[] for _ in range(r)]

    for j in range(r):
        Vi[j]=F2.submodule([F2(x) for x in V_i_gens[j] ]) 
    
    JV_gens=[[] for _ in range(r)]

    for j in range(r): 
        morf = (mats[j][i]).transpose()*matrix( V_i_gens[j] )
        JV_gens[j]=[[F2(morf[m])-F2(V_i_gens[j][m]) for m in range(len(V_i_gens[j]))]  for i in range(len(gens(G)))]+[[p*F2(x) for x in gens(Vi[j])]] # a ação de g deve ser dada por um morfismo de Vj
     

    JV=[[] for _ in range(r)]

    for j in range(r):
        JV[j]=F2.submodule([F2(x) for x in sum( JV_gens[j], [] )])
     
    W=[[] for _ in range(r)]         # lista dos quocientes V_i/JV_i

    D=[]               # listas das dimensões

    for j in range(r): 
        W[j]= Vi[j].V()/JV[j].V()
        d=len(gens(W[j]))
        D.append(d)
    
    s= D[0]+(p-1)*(sum( D[j] for j in range(1,r)))    

    W_gens=[[] for _ in range(r)]

    for j in range(r):
        W_gens[j]=[x for x in gens(W[j])]

    W_lift=[[] for _ in range(r)]

    for j in range(r):
        W_lift[j]=[W[j](x).lift() for x in W_gens[j]]
 
    V_pro=[[] for _ in range(r)] # lista das listas dos levantamentos dos geradores de V_i/JV_i
 
    for j in range(r):
        V_pro[j]=[F2(x) for x in W_lift[j]]
     
    Vi_f_im=[[] for _ in range(2,r)]
 
    for j in range(2,r):
        Vi_f_im[j-2]=[[mudança(p,k, [V_pro[j][m]], V_i_gens[j], mats[j][0]**i)  for i in range(p-1)] for m in range(D[j])]
     
    Vi_f_im=sum(sum(sum(Vi_f_im, [] ), []), [])   
    V0_f_im=[x for x in V_pro[0]]
 
    V1_f_im=[[mudança(p,k, [V_pro[1][m]], V_i_gens[1], mats[1][1]**i)  for i in range(p-1)] for m in range(D[1])]
    V1_f_im=sum(sum(V1_f_im, [] ), [] )
 
    V_f_im=V0_f_im + V1_f_im + Vi_f_im
              
    F3=ZZ**s
    F4=F3.submodule([p**k*gens(F3)[j] for j in range(s)])
    F5=F3/F4
    f=F5.hom([F2(x) for x in V_f_im])
    U0=f.kernel()
    U0_gens=[F5(x) for x in gens(U0)]
    U_gens=[F5(x).lift() for x in U0_gens]+[p**k*y for y in gens(F3)]
    U_mat=matrix(ZZ,U_gens).echelon_form(include_zero_rows=False)
    U=U_mat.rows()
    
    F_act=[[] for _ in range(len(gens(G)))]
    for j in range(len(gens(G))):
        blocks = sum( D[x]*mat_lambda_i[x] for x in range(len(D)), [] )
        F_act[j] = block_diagonal_matrix( blocks ) #mat_act(D,mat_lambda_i[j],r)
        
    U_acts=[[] for _ in range(len(gens(G)))]
    for j in range(len(gens(G))):
        U_acts[j]=act_U(F_act[j],U)
        
    return Reticulado(G,p,R,Vi,mats,k,l,D,U,f,F_act,U_acts,mat_lambda_i)
