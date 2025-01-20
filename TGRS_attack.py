#!/usr/bin/env python
# coding: utf-8

# In[ ]:


reset()
import time
#SET UP
# a is used if the shortening is neccessary
q = 61
n = 60
k = 25
h= 19
t= 17
a=8

#galois field with q elements
Fq = GF(q)
Fqx = PolynomialRing(Fq, x)

# Create a set S of random evaluation points in the field
SSeq = []

#she took it as a list, but turned it into a set.
S = Set(SSeq)

#cardinality counts a,a as 1.
while S.cardinality() != n:
    rd = Fq.random_element()
    if rd != 0:
        SSeq.append(rd)
    S = Set(SSeq)

# Create a generator matrix of Reed Solomon code RS[n,k]
def genMatrix(S, k):
    n = len(S)
    #generated 0 matrix.
    M = matrix(Fq, k, n, 0)
#changed the entries to reed solomon(s,k)
    for i in range(k):
        for j in range(len(S)):
            M[i,j]=Fq(S[j]**i)
    return M


# Schur productprint(M.rank()) of two vectors
def schur(x, y):
    r = [x[i]*y[i] for i in range(n)]  # x ve y nin aynı boyutta vektörler olduğuna dikkat
    return r


############ created a zero row and added the rows, as zero does not contribute to the rank.
def schur_matrix(M,N):
    zero_vec = vector(Fq, n)
    MN= matrix(Fq, zero_vec) 
    for i in M:
        for j in N:
            MN=MN.stack( matrix(Fq, schur(i,j) ) )
    return MN
    
    

#Generator matrix of a twisted Reed Solomon code. Attention : eta = [1,1,..,1] for simplicity
def twistedGenMatrix(S, k, h, t):
    M = genMatrix(S,k)
    for j in range(len(S)):
        M[h,j]=Fq( (S[j]**h)  +(S[j]**(t+k-1))   )
    return M




######### created MONCODE to see if the results are correct at the end.
def Moncode(S, k, h, t):
    M = genMatrix(S,k)
    for j in range(len(S)):
        M[h,j]=Fq( 0)
    return M


M= twistedGenMatrix(S, k, h, t)
C= codes.LinearCode(M)

Moncode=Moncode(S, k, h, t)    
MONCODE=codes.LinearCode(Moncode)



# attack function:
# Input: a TGRS code (or shortening of a TGRS ) 
# Output: a basis of the monomial code
def attack(C):  
    M=C.basis()
    v1= C.random_element() 
    v2= C.random_element() 
    v3= C.random_element() 
    
    V=matrix([v1, v2, v3])
    S=schur_matrix(V,M)

    k=dim(C)
    count=0
    start=time.time()
    while S.rank()>(2*k+2) or V.rank()<3:
        
        v1= C.random_element() 
        v2= C.random_element() 
        v3= C.random_element() 
        V=matrix([v1, v2, v3])
        S=schur_matrix(V,M)

    Basis=matrix([v1, v2, v3])
    basis=3
 
    u=C.random_element() 
    V=matrix([v1, v2, u])
    S=schur_matrix(V,M)
    while basis<k:
        u=C.random_element()
        
        while S.rank()>(2*k+2) :
            u= C.random_element() 
            V=matrix([v1, v2, u])
            S=schur_matrix(V,M)
  
        Basis=Basis.stack(u)
        basis+=1   
        u=C.random_element()
        V=matrix([v1, v2, u])
        S=schur_matrix(V,M)
        
    return Basis


#Shortening of a matrix M at the positions in the array A
def shortening(A,M):
    z= vector(Fq, n)
    gen_matrix =matrix(Fq, z) 
    for a in range(n):
        if a not in A:
            z[a]=1
            gen_matrix=gen_matrix.stack(vector(Fq,z))
            z=vector(Fq, n)      
    D=codes.LinearCode(gen_matrix) 
    
    s=[]
    for d in D.basis():
        s.append(d)
    S=span(s)
    t=[]
    #t= BASIS of TGRS
    for d in M:
        t.append(d)
    T=span(t)
    #ST= shortened code
    ST=T.intersection(S)
    return(codes.LinearCode(ST))

begin=time.time()
if schur_matrix(M,M).rank()<n:
    attack(C)
else:
    shortening_array =  [0,..,a-1]
    D=attack(shortening(shortening_array,M))
    union_of_shortenings=D
    count=0
    while union_of_shortenings.rank()<k-2:
        shortening_array = [count*a ,.., (count+1)*a-1]
        D=attack(shortening(shortening_array,M))
        union_of_shortenings=union_of_shortenings.stack(D)

end=time.time()   
print('The implementtion of the attack take', end-begin, 'seconds.')   


#checking whether the monomial subcode is recovered
#recovered_code=codes.LinearCode(union_of_shortenings)
#recovered_code==MONCODE  



# we will compute ( C* (C*^2 )^perp  )^perp to recover the GRS code containing the monomial code

#R: Basis of monomial code
R=recovered_code.basis()

#R2: Basis of the dual of Schur Square
Schur_square=codes.LinearCode(schur_matrix(R,R))
Schur_square_perp=Schur_square.dual_code()
R2=Schur_square_perp.basis()

# GRS_code_dual= ( C* (C*^2 )^perp  ), GRS_code= ( C* (C*^2 )^perp  )^perp 
GRS_code_dual_generator=schur_matrix( R,R2 )
GRS_code_dual=codes.LinearCode( GRS_code_generator  )
GRS_code=GRS_code_dual.dual_code()




