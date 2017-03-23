from random import randint

class McEliecePKS:
    def __init__(self, m, t):
        self.m = m
        self.t = t
        self.n = 2**m
        self.k = self.n-m*t
        self.pub_key = None
        self.priv_key = None
        self.goppacode = BinGoppaCode(m,t)
        
        gf2.<a> = GF(2)
        gf2m.<a> = GF(2**m)
        ## Generating scramble matrix S
        while True:
            S = matrix(gf2, [[randint(0,1) for i in range(self.goppacode.G_Goppa.nrows())]for i in range(self.goppacode.G_Goppa.nrows())])
            if S.is_invertible():
                break
            
        ## Generating permutation matrix P
        rng = range(n); P = matrix(GF(2),n);
        for i in range(n):
            p = floor(len(rng)*random());
            P[i,rng[p]] = 1; rng=rng[:p]+rng[p+1:];
            
        G_ =   S*self.goppacode.G_Goppa*P
        self.pub_key = (G_, self.t)
        self.priv_key = (self.goppacode.gpoly, self.goppacode.G_Goppa, S, P)
        
    def encrypt(self, m, pub_key):
        gf2 = GF(2)
        G_ = pub_key[0]
        t = pub_key[1]
        c = m*G_
        # random error vector 
        e = [0 for i in range(G_.ncols())]
        for i in range(t):
            e[randint(0,G_.ncols()-1)] = 1
        e = matrix(gf2, e)
        c = c + e
        return c
    
    def decrypt(self, c):
        priv_key = self.priv_key
        P = priv_key[3]
        S = priv_key[2]
        gpoly = priv_key[0]
        G_Goppa = priv_key[1]
        
        #Getting inverse of P
        #m = P.join_with(Matrix.get_identity(P.rows))
        #m = m.get_reduced_echelon()
        #m = m.submatrix(P.rows,P.rows,m.rows,m.cols)
        #m = m*c
        m = c*(~P)
        
        # Correcting errors - Solving system
        m = self.goppacode.decode(m)
        #m = (gpoly.transpose()).join_with(m.transpose())
        #m = m.get_reduced_echelon()
        ## finding inverse of S
        #s_ = S.join_with(Matrix.get_identity(S.rows))
        #s_ = m.get_reduced_echelon()
        #s_ = m.submatrix(S.rows,S.rows,m.rows,m.cols)
        #m = s_*m
        m = (S*G_Goppa).solve_left(m)
        
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
        
        #Generating galois field
        gf2.<a> = GF(2)
        gf2m.<a> = GF(2**self.m)
        self.a = a
        self.gf2 = gf2
        self.gf2m = gf2m
        
        #Choosing n elements for L, choosing from power 2 to power 2+n 
        L = [a**i for i in range(2,2+self.n)]
        self.L = L
        
        #Generating goppa polynomial - fixed form for now
        x = PolynomialRing(gf2m,repr(a)).gen()
        gpoly = x^t + x + L[4]
        self.gpoly = gpoly
        #Checking that satisfies Goppa polynomial conditions e.g.Irreducible and elements of L not roots
        if not gpoly.is_irreducible():
            print "not irreducible"
            return    # try another polynomial check for each L not being a root
        #Generate parity matrix H
        for elem in L:
            if gpoly(elem) == gf2m(0):
                print "an element from L is a root of g"
                return
                
        #Constructing parity matrix
                
        self.H_gRS = matrix([[L[j]^(i) for j in range(self.n)] for i in range(self.t)])
        self.H_gRS = self.H_gRS*diagonal_matrix([1/gpoly(L[i]) for i in range(self.n)])
        self.H_Goppa = matrix(gf2,self.m*self.H_gRS.nrows(),self.H_gRS.ncols())
        
        for i in range(self.H_gRS.nrows()):
            for j in range(self.H_gRS.ncols()):
                be = bin(self.H_gRS[i,j].integer_representation())[2:];
                be = be[::-1];
                be = be+'0'*(m-len(be));
                be = list(be);
                self.H_Goppa[m*i:m*(i+1),j] = vector(map(int,be));
                
        self.G_Goppa = self.H_Goppa.transpose().kernel().basis_matrix();
        G_Goppa_poly = self.H_gRS.transpose().kernel().basis_matrix();		
        
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
        while (2^((a[i]^2+X*b[i]^2).degree())) > 2^t: # while(norm(a[i],b[i]) > 2^t)
            a.append(0); b.append(0);
            (q,r) = a[i-1].quo_rem(a[i]);
            (a[i+1],b[i+1]) = (r, b[i-1] - q*b[i]);
            i+=1;
        return (a[i],b[i]);
    
    def decode(self, y_):
   # Decoding using Patterson's Algorithm
        y = copy(y_);
        X = self.gpoly.parent().gen();
        synd = self.SyndromeCalculator*y.transpose();		
    
        syndrome_poly = 0;
        for i in range (synd.nrows()):
            syndrome_poly += synd[i,0]*X^i
        
        error = matrix(GF(2),1,self.H_Goppa.ncols());
        (g0,g1) = self._split(self.gpoly); 
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
        #Pre-test sigma if fails, then zerro error vector is returned
        if (X^(2^self.m)).mod(sigma) != X:	
            return y+error; 
        #Generating the error correcting vector
        for i in range(len(self.L)):
            if sigma(self.L[i]) == 0:
                error[0,i] = 1;
        
        return y+error;


# Main

m=5    
t=4
n=2**m
k=n-m*t

## McEliece Testing
testE = McEliecePKS(m,t)
#print "priv:", testE.priv_key 
#print "pub:", testE.pub_key
gf2 = GF(2)
m = matrix(gf2,[randint(0,1) for i in range(k)])
print "Message:"
print "m: ",m
c = testE.encrypt(m, testE.pub_key)
print "Encrypting m"
print "c: ",c
d = testE.decrypt(c)
print "Decrypting c"
print "d :",d

## Goppa Code Testing

#testC = BinGoppaCode(m,t)
#print  "G:",testC.G_Goppa,"\n"
#print "H:",testC.H_Goppa, "\n"

#gf2 = GF(2)
#m = matrix(gf2,[randint(0,1) for i in range(k)])
#print "m: ",m
#c = testC.encode(m)
#print "c: ",c
#e = [0 for i in range(c.ncols())]
#for i in range(t):
    #e[randint(0,c.ncols()-1)] = 1
#e = matrix(gf2, e)
#print "e:",e
#ce = c + e
#print "ce:",ce
#d = testC.decode(ce)
#print "d :",d


        
