from galois import *
from coding import *

class McEliecePKS:
    def _init_(self, n, k , t):
        
        self.t = t
        self.k = k 
        self.m = (n-k)/t
        self.n = n
        self.pub_key = None
        self.priv_key = None
        goppacode = BinGoppaCode(n,k,t)
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
    
        
            
            
        
        
        
        
        
