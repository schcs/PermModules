class Reticulado:
    """The data structure to store the information related to a Butler diagram.
    If d is a Butler diagral, then d has the following components:
    d.group:                    the acting p-group G
    d.p:                        the prime p
    d.padic_ring:               the underlying p-adic ring
    self.Vi:                    the list of spaces Vi in the diagram
    self.action_Vi:             the sequence of matrices defining the action of G on the Vi
    self.k:                     o menor inteiro que podemos mergulhar V em um Z/p**kZ-modulo livre 
    self.ranksubmodule:         menor l em que podemos ver V como submodulo de (Z/p**kZ)**l
    self.D:                     lista das dimensões de Vi/JVi
    self.reticulado:            lista dos geradores de U
    self.morphism:              o homomorfismo phi tal que U é seu kernel
    self.F_act:                 a lista de matrizaes que definem a ação de G em F >= U
    self.U_acts:                a lista de matrizes que definem a ação de G em U
    self.mat_lambda_i:          a lista das listas de matrizes das ações em cada  ZpGe_i para cada gerador de G
    """
    def __init__(self, G, p, R, Vi,mats,k,l,D,U,f,F_act, U_acts, mat_lambda_i):
        self.group = G
        self.p = p
        self.padic_ring=R
        self.Vi=Vi
        self.action_Vi=mats
        self.k=k
        self.ranksubmodule=l
        self.D=D
        self.reticulado=U
        self.morphism=f
        self.F_act=F_act
        self.U_acts=U_acts
        self.mat_lambda_i=mat_lambda_i
def __repr__( self ):
    return f'lattice U associated with the Butler diagram for the group {self.group} over {self.padic_ring}'    


# a função que devolve a lista das imagens dos geradores de V_i pela multiplicação por um elemento de G
# isn't this function just matrix multiplication by g?

# g is a matrix, geradores é uma lista de elementos de ZZ^m. 
# returns [ ]
def morfis(geradores, g):

    # parent of geradores
    F=ZZ**len(geradores[0])

    # g is k x k matrix
    k = len(g.rows()[0])

    # m is an array of length k
    m=[[] for _ in range(k)]

    for j in range(k):
        # this looks like vector multiplication 
        m[j]=sum(g[i][j]*F(geradores[i]) for i in range(k))

    return m


# parece que morfis faz a mesma coisa que a seguinte função
def morfis2( geradores, g ):
    return g.transpose()*matrix( geradores )

# try sum( L, [] )
def sumlist(L):
    l=[]
    for j in range(len(L)):
        l=l+L[j]
    return l


def mat_act(D,act_mats,r):
    l=[[] for j in range(r)]
    for j in range(r):
        l[j]=[act_mats[j] for i in range(D[j])]
    l0=sumlist(l)
    m=block_diagonal_matrix(l0)
    return m


def act_U(mat,gen_U):
    l=[mat*gen_U[j] for j in range(len(gen_U))]
    m=transpose(matrix(QQ,gen_U))
    coor=[[] for _ in range(len(l))]
    for j in range(len(l)):
        coor[j]=m.solve_right(vector(QQ,l[j]))
    mact=transpose(matrix(QQ,coor))
    return mact
# p é primo
# k >= 1
# gen1 é lista de vetores com entradas em Z/(p^k) de comprimento l
# gen2 é lista de vetores com entradas em Z/(p^k) de comprimento l
# os elementos de gen1 pertencem ao módulo gerado por gen2
# mat é matriz quadrada m x m onde m == len( gen2 ) com entradas em Z/(p^k)
#
# output: devolve uma lista das imagens dos vetores em gen1 pelo morfismo V -> V determinado pela matriz mat 

def mudança(p,k,gen1,gen2,mat):# gen1 é lista de elementos de um módulo que desejamos conhcer a imagem por mat e gen2 os geradores desse módulo que ja conhecemos a imagem por mat
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

def lattice(G,k,r,mats, V_i_gens, mat_lambda_i):

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
        JV_gens[j]=[[F2(morfis(V_i_gens[j],mats[j][i])[m])-F2(V_i_gens[j][m]) for m in range(len(V_i_gens[j]))]  for i in range(len(gens(G)))]+[[p*F2(x) for x in gens(Vi[j])]] # a ação de g deve ser dada por um morfismo de Vj
     

    JV=[[] for _ in range(r)]

    for j in range(r):
        JV[j]=F2.submodule([F2(x) for x in sumlist(JV_gens[j])])
     
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
     
    Vi_f_im=sumlist(sumlist(sumlist(Vi_f_im)))   
    V0_f_im=[x for x in V_pro[0]]
 
    V1_f_im=[[mudança(p,k, [V_pro[1][m]], V_i_gens[1], mats[1][1]**i)  for i in range(p-1)] for m in range(D[1])]
    V1_f_im=sumlist(sumlist(V1_f_im))
 
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
        F_act[j]=mat_act(D,mat_lambda_i[j],r)
        
    U_acts=[[] for _ in range(len(gens(G)))]
    for j in range(len(gens(G))):
        U_acts[j]=act_U(F_act[j],U)
        
    return Reticulado(G,p,R,Vi,mats,k,l,D,U,f,F_act,U_acts,mat_lambda_i)
