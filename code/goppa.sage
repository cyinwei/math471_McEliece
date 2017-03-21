from coding import *

class McEliecePKS:
    def _init_(self, m=4, k=12):
        self.k = k
        self.m = m
        self.t = int((n-k))
        self.n = 2^m-1
        self.pub_key = None
        self.priv_key = None
        goppacode = BinGoppaCode(m,k)
        # Generating scramble matrix S
        data = [[random.randint(0,1) for i in range(len(goppacode.gen_mx))]for i in range(len(goppacode.gen_mx))]
        S = Matrix(data=data).to_Zmod(2)
        # Generating permutation matrix P
        data = [[random.randint(0,1) for i in range(len(goppacode.gen_mx[0]))]for i in range(len(goppacode.gen_mx[0]))]
        P = Matrix(data=data).to_Zmod(2)
        G_ =   S*P*goppacode.gen_mx
        self.pub_key = (G_, self.t)
        self.priv_key = (goppacode.gpoly, goppacode.gen_mx, S, P)
        
    def encrypt(self, m, pub_key):
        c = m*(pub_key[0])
        # random error vector 
        e = [0 for i in range(len(pub_key[0]))]
        for i in range(t):
            e[random.randint(0,len(pub_key[0]))] = 1
        c = c + e
        return c
    
    def decrypt(self, c):
        #Getting inverse of P
        m = priv_key[3].join_with(Matrix.get_identity(priv_key[3].rows))
        m = m.get_reduced_echelon()
        m = m.submatrix(priv_key[3].rows,priv_key[3].rows,m.rows,m.cols)
        m = m*c
        # Correcting errors
        m = goppacode.decode(m)
        # Final calculations
        m = (priv_key[0].transpose()).join_with(m.transpose())
        m = m.get_reduced_echelon()
        # finding inverse of S
        s_ = priv_key[2].join_with(Matrix.get_identity(priv_key[2].rows))
        s_ = m.get_reduced_echelon()
        s_ = m.submatrix(priv_key[2].rows,priv_key[2].rows,m.rows,m.cols)
        m = s_*m
        return m

class BinGoppaCode:
    def _init_(self, m=4, k=12):
        self.k = k
        self.m = m
        self.t = int((n-k))
        self.n = 2^m-1
        self.gen_mx = [[]]
        self.parity_mx = [[]]
        self.L = []
        self.poly = []
        
        #Constructing parity matrix
        #Generating galois field
        gf2.<a> = GF(2)
        gf2m.<a> = GF(2**m)
        self.gf2 = gf2
        self.gf2m = gf2m
        
        #Choosing n elements for L, choosing from power 2 to power 2+n 
        L = [a**i for i in range(2,n+2)]
        self.L = L
        #Generating goppa polynomial - fixed for now
        
        x = PolynomialRing(gf2m,repr(a)).gen()
        gpoly = x^t + x + L[4]
        
        #Checking that satisfies Goppa polynomial conditions e.g.Irreducible and elements of L not roots
        if not gpoly.is_irreducible():
            return 1    # try another polynomial check for each L not being a root
        #Generate parity matrix H
        for elem in L:
            if gpoly(elem) == gf2m(0):
                return 1
                
        #gl_ = [(g(elem))**-1 for elem in L] 
        #H = [[()*elem for elem in gl_] for i in range(1,t)]
        H_gRS = matrix([[L[j]^(i) for j in range(n)] for i in range(t)])
        H_gRS = H_gRS*diagonal_matrix([1/gpoly(L[i]) for i in range(n)])
        
        self.H_Goppa = matrix(gf2,m*H_gRS.nrows(),H_gRS.ncols())
        for i in range(H_gRS.nrows()):
            for j in range(H_gRS.ncols()):
                be = bin(eval(H_gRS[i,j].int_repr()))[2:]
                be = '0'*(m-len(be))+be; be = list(be)
                H_Goppa[m*i:m*(i+1),j]=vector(map(int,be))
                
        Krnl = H_Goppa.right_kernel()
        self.G_Goppa = Krnl.basis_matrix()
        
    def encode(self, u):
        return u*G_Goppa
    
    def decode(self, y):
        z = var('z')
        PR = PolynomialRing(gf2m,'z')
        # init integer bigN
        bigN = D-1
        # declare sigma's (sigma_{-1} ... sigma_bigN)
        # declare omega's (omega_{-1} ... omega_bigN) 
        sigma = vector(PR,bigN+2)
        omega = vector(PR,bigN+2) 
        delta = vector(gf2m,bigN+2)
        # init sigma_{-1} and sigma_0 as well as omega_{-1} and omega_0 
        sigma[-1+1] = PR(0)
        sigma[0+1] = PR(1)
        flag = 2*bigN # z^flag represents the rational function 1/z
        omega[-1+1] = z^(flag)
        omega[0+1] = PR(0)
        # init mu and delta
        mu = -1
        delta[-1+1] = 1
        
        
        s = H_gRS*y.transpose()
        if s==matrix(self.gf2m,H_gRS.nrows(),1):
            return y;
        
        b = PR([s[_,0] for _ in range(s.nrows())])
        # init sigma_{-1} and sigma_0 as well as omega_{-1} and omega_0 
        sigma[-1+1] = PR(0)
        sigma[0+1] = PR(1)
        flag = 2*bigN # z^flag represents the rational function 1/z
        omega[-1+1] = z^(flag)
        omega[0+1] = PR(0)
        # init mu and delta
        mu = -1
        delta[-1+1] = 1
        for i in range(bigN):
            delta[i+1] = (sigma[i+1]*b).coeffs()[i]
            sigma[i+1+1] = sigma[i+1](z)-(delta[i+1]/delta[mu+1])*z^(i-mu)*sigma[mu+1](z);
            if (omega[mu+1].degree()==flag):
                omega[i+1+1] = omega[i+1](z)-(delta[i+1]/delta[mu+1])*z^(i-mu-1);
            else:
                omega[i+1+1] = omega[i+1](z)-(delta[i+1]/delta[mu+1])*z^(i-mu)*omega[mu+1](z)
            rord = max(sigma[i+1].degree(),1+omega[i+1].degree()) # recurrence order
            if (delta[i+1]<>0)and(2*rord<=i):
                mu = i
                
        ELP = sigma[bigN+1] # the Error Locator Polynomial
        # compute zeroes of ELP, compute error positions, compute error vector ee
        ee = vector(self.gf2,[0 for _ in range(n)])
        for i in range(self.N):
            if (ELP(x^i)==self.gf2m(0)): #an error occurred 
                ee[mod(self.N-i,self.N)] += 1 # in position N-i
        cc = y+ee
        return cc
    
    
pks = McEliecePKS()
testC = BinGoppaCode()
F.<x> = GF(2)
u = vector(F,[randint(0,1) for i in range(12)])
print testC.encode(u)

        
