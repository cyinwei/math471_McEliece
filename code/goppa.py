from galois import *
from coding import *

def evalPol(pol,a):
    i = 0
    res = 
    for coff in pol.i.coefficients:
        res += coff*(a ** i)
        
class BinGoppaCode:
    def _init_(self, n, k, t):
        self.t = t
        self.k = k
        self.m = (n-k)/t
        self.n = n
        self.gen_mx = [[]]
        self.parity_mx = [[]]
        self.L = []
        self.poly = []
        
        #Constructing parity matrix
        #Generating galois field
        gfm = GF(2**m)
        gf2 = GF(2)
        #Choosing n elements for L, choosing from power 2 to power 2+n 
        L = [gfm[i] for i in range(3,3+n)]
        #Generating goppa polynomial - fixed for now
        coffs [gfm[0] for i in range(t)]
        coffs[-1] = gfm[1]
        coffs[0] = gfm[7]
        coffs[1] = gfm[1]
        gpoly = Polynomial(coffs)
        #Checking that satisfies Goppa polynomial conditions e.g.Irreducible and elements of L not roots
        Zmodx = Zmod(2**m)
        Zmodx = [Polynomial(list(reversed(x))) for x in Zmodx]
        if is_reducable(gpoly, Zmodx):
            return 1    # try another polynomial check for each L not being a root
        #Generate parity matrix H
        [for elem in L]
        # calculate corrects gs!
        
        
    def encode(self):
    def decode(self):
        
        
