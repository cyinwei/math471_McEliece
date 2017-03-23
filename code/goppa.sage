#from coding import *

class McEliecePKS:
    def __init__(self, m, k):
        self.k = k
        self.m = m
        self.n = (2**self.m)-1
        self.t = (self.n-self.k)
        self.pub_key = None
        self.priv_key = None
        goppacode = BinGoppaCode(self.m,self.k)
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
    def __init__(self, m, k):
        self.k = 12
        self.m = 4
        self.n = (2**m)-1
        self.t = (self.n-self.k)
        self.D = self.n-self.k+1
        self.gen_mx = [[]]
        self.parity_mx = [[]]
        self.L = []
        self.poly = []
        
        #Constructing parity matrix
        #Generating galois field
        gf2.<a> = GF(2)
        gf2m.<a> = GF(2**self.m)
        self.a = a
        self.gf2 = gf2
        self.gf2m = gf2m
        
        #Choosing n elements for L, choosing from power 2 to power 2+n 
        L = [a**i for i in range(0,self.n)]
        self.L = L
        #Generating goppa polynomial - fixed for now
        
        x = PolynomialRing(gf2m,repr(a)).gen()
        gpoly = x^(3) + x + 1
        
        #Checking that satisfies Goppa polynomial conditions e.g.Irreducible and elements of L not roots
        if not gpoly.is_irreducible():
            print "holi irre1"
            return    # try another polynomial check for each L not being a root
        #Generate parity matrix H
        for elem in L:
            if gpoly(elem) == gf2m(0):
                print "holi irre2"
                return
                
        #gl_ = [(g(elem))**-1 for elem in L] 
        #H = [[()*elem for elem in gl_] for i in range(1,t)]
        self.H_gRS = matrix([[L[j]^(i) for j in range(self.n)] for i in range(3)])
        self.H_gRS = self.H_gRS*diagonal_matrix([1/gpoly(L[i]) for i in range(self.n)])
        
        self.H_Goppa = matrix(gf2,self.m*self.H_gRS.nrows(),self.H_gRS.ncols())
        for i in range(self.H_gRS.nrows()):
            for j in range(self.H_gRS.ncols()):
                be = bin(eval(self.H_gRS[i,j].int_repr()))[2:]
                be = '0'*(self.m-len(be))+be; be = list(be)
                self.H_Goppa[self.m*i:self.m*(i+1),j]=vector(map(int,be))
        
        self.nl = self.H_Goppa.ncols(); 
        Krnl = self.H_Goppa.right_kernel()
        self.G_Goppa = Krnl.basis_matrix()
        
        
    def encode(self, u):
        return u*self.G_Goppa
    
    def decode(self, y):
        z = var('z')
        PR = PolynomialRing(self.gf2m,'z')
        # init integer bigN
        bigN = self.D-1
        # declare sigma's (sigma_{-1} ... sigma_bigN)
        # declare omega's (omega_{-1} ... omega_bigN) 
        sigma = vector(PR,bigN+2)
        omega = vector(PR,bigN+2) 
        delta = vector(self.gf2m,bigN+2)
        # init sigma_{-1} and sigma_0 as well as omega_{-1} and omega_0 
        sigma[-1+1] = PR(0)
        sigma[0+1] = PR(1)
        flag = 2*bigN # z^flag represents the rational function 1/z
        omega[-1+1] = z^(flag)
        omega[0+1] = PR(0)
        # init mu and delta
        mu = -1
        delta[-1+1] = 1
        
        
        s = (self.H_gRS)*(y.transpose())
        if s==matrix(self.gf2m,self.H_gRS.nrows(),1):
            print "rati"
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
        ee = matrix(self.gf2,[0 for _ in range(self.nl)])
        for i in range(self.n):
            if (ELP(x**i))==self.gf2m(0): #an error occurred 
                print "sii"
                ee[mod(self.n-i,self.n)] += 1 # in position N-i
        cc = y+ee
        print "y:",y
        print "e:",ee
        return cc
    

testC = BinGoppaCode(4,12)
print  "G:",testC.G_Goppa,"\n"
print "H:",testC.H_Goppa, "\n"

gf2 = GF(2)
c = testC.encode(matrix(gf2, [0,1,0] ))
print "c:",c
ce = c + matrix(gf2, [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1])
print "ce:",ce
d = testC.decode(ce)

print "d:",d


        
