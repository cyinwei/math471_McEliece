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
    def __init__(self, m, t):
        self.m = m
        self.t = t
        self.n = 2**m
        self.k = self.n-m*t
        self.gen_mx = [[]]
        self.parity_mx = [[]]
        self.L = []
        self.D = self.n-self.k+1
        
        #self.m = 4
        #self.t = 3
        #self.n = (2**m)
        #self.k = 12
        
        #Generating galois field
        gf2.<a> = GF(2)
        gf2m.<a> = GF(2**self.m)
        self.a = a
        self.gf2 = gf2
        self.gf2m = gf2m
        
        #Choosing n elements for L, choosing from power 2 to power 2+n 
        L = [a**i for i in range(2,2+self.n)]
        self.L = L
        #Generating goppa polynomial - fixed for now
        
        x = PolynomialRing(gf2m,repr(a)).gen()
        gpoly = x^t + x + L[4]
        print gpoly
        self.gpoly = gpoly
        #Checking that satisfies Goppa polynomial conditions e.g.Irreducible and elements of L not roots
        if not gpoly.is_irreducible():
            print "holi irre1"
            return    # try another polynomial check for each L not being a root
        #Generate parity matrix H
        for elem in L:
            if gpoly(elem) == gf2m(0):
                print "holi irre2"
                return
                
        #Constructing parity matrix
                
        self.H_gRS = matrix([[L[j]^(i) for j in range(self.n)] for i in range(self.t)])
        self.H_gRS = self.H_gRS*diagonal_matrix([1/gpoly(L[i]) for i in range(self.n)])
        
        self.H_Goppa = matrix(gf2,self.m*self.H_gRS.nrows(),self.H_gRS.ncols())
        
        for i in range(self.H_gRS.nrows()):
            for j in range(self.H_gRS.ncols()):
                #print 'cnt=',i*n+j;
                #print '__init__: poly=',H_check_poly[i,j];
                be = bin(self.H_gRS[i,j].integer_representation())[2:];
                be = be[::-1];
                #print '__init__: be=',H_check_poly[i,j].integer_representation();
                #be = '0'*(m-len(be))+be; 
                be = be+'0'*(m-len(be));
                be = list(be);
                self.H_Goppa[m*i:m*(i+1),j] = vector(map(int,be));
                #print "H_Goppa ",H_Goppa[m*i:m*(i+1),j].transpose();
        
        self.G_Goppa = self.H_Goppa.transpose().kernel().basis_matrix();
        G_Goppa_poly = self.H_gRS.transpose().kernel().basis_matrix();		
                
        #Construct the syndrome calculator. This will be used
        #to simplify the calculation of syndromes for decoding.
        PR_F_2m = gpoly.parent()
        X = PR_F_2m.gen();
        self.SyndromeCalculator = matrix(PR_F_2m, 1, len(L));
        for i in range(len(L)):
            self.SyndromeCalculator[0,i] = (X - L[i]).inverse_mod(gpoly);
        
    def encode(self, u):
        return u*self.G_Goppa
    
    def _split(self,p):
        # split polynomial p over F into even part po
        # and odd part p1 such that p(z) = p2 (z) + z p2 (z)
        Phi = p.parent()
        p0 = Phi([sqrt(c) for c in p.list()[0::2]]);
        p1 = Phi([sqrt(c) for c in p.list()[1::2]]);
        return (p0,p1);
    
    def _lattice_basis_reduce(self, s):
        # a <- s   b <- v 
        g = self.gpoly;
        t = g.degree();
        a = []; a.append(0);
        b = []; b.append(0);
        (q,r) = g.quo_rem(s);
        (a[0],b[0]) = simplify((g - q*s, 0 - q))

        #If the norm is already small enough, we
        #are done. Otherwise, intialize the base
        #case of the recursive process.
        #print 'norm: ',self._norm(a[0],b[0]);
        
        #This is the way in which Bernstein indicates
    #the norm of a member of the lattice is
    #to be defined.
        X = g.parent().gen();
        norm = 2^((a[0]^2+X*b[0]^2).degree());
        
        if norm > 2^t:
            a.append(0); b.append(0);
            (q,r) = s.quo_rem(a[0]);
            (a[1],b[1]) = (r, 1 - q*b[0]);
            if a[1] == 0:
                return (s,1);	
        else:
            return (a[0], b[0]);
        i = 1;
        while (2^((a[i]^2+X*b[i]^2).degree())) > 2^t:   # while(norm(a[i],b[i]) > 2^t)
            a.append(0); b.append(0);
            (q,r) = a[i-1].quo_rem(a[i]);
            (a[i+1],b[i+1]) = (r, b[i-1] - q*b[i]);
            i+=1;
        return (a[i],b[i]);
    
    def decode(self, y_):
    
        y = copy(y_);
        X = self.gpoly.parent().gen();
        synd = self.SyndromeCalculator*y.transpose();		
        #print 'Decode: synd=', synd;
    
        syndrome_poly = 0;
        for i in range (synd.nrows()):
            syndrome_poly += synd[i,0]*X^i
        
        #We will decode codewords using Pattersons Algorithm, errorlocator polynomial returned.
        error = matrix(GF(2),1,self.H_Goppa.ncols());
        (g0,g1) = self._split(self.gpoly); 
        # returns the g-inverse of polynomial p
        (d,u,v) = xgcd(g1,self.gpoly)
        g1_inverse = u.mod(self.gpoly)
        sqrt_X = g0*g1_inverse;
        T = syndrome_poly.inverse_mod(self.gpoly);
        (T0,T1) = self._split(T - X);
        R = (T0+ sqrt_X*T1).mod(self.gpoly);
        #Perform lattice basis reduction.
        (alpha, beta) = self._lattice_basis_reduce(R);
        #Construct the error-locator polynomial.
        sigma = (alpha*alpha) + (beta*beta)*X;
        #print 'SyndromeDecode: sigma=', BinRepr(sigma);
        #Pre-test sigma
        if (X^(2^self.m)).mod(sigma) != X:	
            #print "sigma: Decodability Test Failed";
            return y+error;# return a zero vector 
        #For every root of the error polynomial,
        #correct the error induced at the corresponding index.
        for i in range(len(self.L)):
            if sigma(self.L[i]) == 0:
                error[0,i] = 1;
        
        
        return y+error;
    
    

testC = BinGoppaCode(4,2)
print  "G:",testC.G_Goppa,"\n"
print "H:",testC.H_Goppa, "\n"

gf2 = GF(2)
c = testC.encode(matrix(gf2, [0,0,0,0,0,0,0,1] ))
print "c: ",c
ce = c + matrix(gf2, [0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0])
print "ce:",ce
d = testC.decode(ce)

print "d :",d


        
