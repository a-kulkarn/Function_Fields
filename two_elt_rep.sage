from mac_lane import * #/Users/JB/Dropbox/programming/Sage/Function_Fields/mac_lane
#from sage.rings.padics.padic_generic import pAdicGeneric


def ppio(a,b):
    '''
    Seperates the composite integer `a` into a product a = cn such that gcd(c,b) = 1.
    Returns (c,n)
    '''
    c = gcd(a,b)
    n = a//c
    m = c # not necessary
    g=  gcd(c,n)
    while g!=1:
            c = c*g
            n = n//g
            g = gcd(c,n)
    return c,n

def reduce(alpha,n):
    '''
    takes alpha and reduces it mod integer n
    '''
    den_alpha = denominator(alpha)
    red = (n*den_alpha).numerator()
    return alpha.parent([j.numerator() % red for j in numerator(alpha).list()])/den_alpha

class dedekind_integral_ideal(object):

    def __init__(self, min_gen,gen,is_prime,e,f,val,norm,ideal_basis):
        self.min_gen = min_gen
        self.gen = gen
        # (min_gen, gen) is the ideal in min_gen - normal representation
        #self.index = index
        self.is_prime = is_prime
        self.e = e
        self.f = f
        self.val = val
        self.norm = norm
        self.ideal_basis = ideal_basis

#	def __eq__(self,I):

    def inver(self):
        inv = self.gen^-1
        den = denominator(inv)
        gen = den/gcd(den,self.min_gen)*inv
        d = denominator(gen) # Should be denominator in maximal order
        min_gen = 1
        I = dedekind_integral_ideal(min_gen,reduce(d*gen,d^2),False,None,None,None,self.norm^-1,None)
        return dedekind_ideal(I,self.min_gen, False,None)     
       
    def __repr__(self):
        return"("+str(self.min_gen)+"," +str(self.gen)+")"
    
    def expand_support_gen(self,T):
        a = self.min_gen.numerator()
        #print"test",T.parent(),a.parent()
        t,r = ppio(T.numerator(),a)
        print"test1",t,r
        eins, u, v = xgcd(a*a,t)
        #print"test2"        
        beta = v*t*self.gen+u*a*a
        #print"test3"        
#        beta = reduce(beta,a^2)
        #print"test4"        
        return beta
    
    def __mul__(self,I):
        '''
        We assume that maximal order = equation order 
        '''
        # TODO catch multiplication by trivial ideal
        # TODO invert prod_gen wrt maximal order
        a = self.min_gen
        alpha = self.expand_support_gen(I.min_gen)
        b = I.min_gen
        beta = I.expand_support_gen(self.min_gen)
        prod_min = a*b
        print"what is zero",alpha,I.min_gen
        prod_gen = alpha*beta
        prod_min = gcd(denominator(prod_gen^-1).numerator(),prod_min.numerator())
        norm = self.norm*I.norm    
        return dedekind_integral_ideal(prod_min,prod_gen,False,None,None,None,norm,None)

#        self.factors = factors

def Valuation(K,p):
    '''
    Produces a valuation on p.parent() determined by prime elements p
    '''
    if not p.is_prime:
        raise Exception, "Error: p must be prime."
    if K in NumberFields():
        return pAdicValuation(QQ, p)
    elif K in FunctionFields():
        return FunctionFieldValuation(K.base_field(),p)
    elif isinstance(K,pAdicGeneric):
        return lambda x: K(x).valuation()
    else:
        raise NotImplementedError

def decompose_prime(K,p):
    '''
    Takes a prime number p.
    Returns the Dedekind ideal representation of all primes in K lying over p.
    '''
    v = Valuation(K,p)
    val_ext = v.extensions(K)
    normalized_uniformizer = all_normalized_uniformizer(val_ext)
    normalized_uniformizer = [reduce(elt, p^2) for elt in normalized_uniformizer]
 #   print"elt",normalized_uniformizer[1]
    list_of_prime_ideals = []
    for i in range(len(val_ext)):
        min_gen = p
        gen = normalized_uniformizer[i]
        is_prime = True
        e = val_ext[i]._base_valuation._initial_approximation.E()
        f = val_ext[i]._base_valuation._initial_approximation.F()
        val = e*val_ext[i]
        norm = p**f
        ideal_basis = None    
        P = dedekind_integral_ideal(min_gen,gen,is_prime,e,f,val,norm,ideal_basis)
        list_of_prime_ideals.append(P)
    return list_of_prime_ideals


def denominator(alpha):
    '''
    Returns the denominator of alpha, with respect to the maximal order.
    '''
    if alpha.parent() in NumberFields():
    	return alpha.denominator()
    # TODO: Sanitize code to ensure correct inputs.
    return lcm([co.denominator() for co in alpha.list()])


        
def numerator(alpha):
    '''
    Returns the numerator of alpha, with respect to the maximal order.
    '''
    return alpha*denominator(alpha)



        
def all_normalized_uniformizer(val_ext):
    """
    Given the list of valuations ``val_ext`` extending a valuation `v` on a field K,
    Return a list `L = [pi]` of uniformizers such that
        
        val_ext[i]( L[j] ) = delta_ij

    where delta_ij is the Kronecker delta.  
    """
    [w._base_valuation._improve_approximation() for w in val_ext]# should be improved
    s = len(val_ext)
    ram_ind = [w._base_valuation._initial_approximation.E() for w in val_ext]
    theta = val_ext[0].domain().0
    L = range(s)
    phi_pols = [w._base_valuation._approximation.phi()(val_ext[0].domain().0) for w in val_ext]
 #   print"phi_pols",phi_pols
    phi_power = [phi_pols[i]**ram_ind[i] for i in L]
    phi_val_mat = [[val_ext[i](phi)*ram_ind[i]  for i in L] for phi in phi_pols]
    normilized_uniformizer = []
#    correction_term = 0
    for ind in L:
        run = L[:ind]+L[ind+1:]
        local_sols = []

        for ind_to_zero in L:
            L_index = L[:ind_to_zero]+ L[ind_to_zero+1:]
            p_power = sum([phi_val_mat[j][ind_to_zero] for j in L_index])
            red_p = p**p_power
            reduce_pow = red_p*p^2
            if ram_ind[ind_to_zero] > 1:
                reduce_pow *= p


            elt_ind = prod([phi_power[l] for l in L_index])
            elt_red = elt_ind		#reduce(elt_ind,reduce_pow.numerator()) # has to be adjusted --> better reduction

            #elt_ind.parent([i.numerator() % reduce_pow.numerator() for i in elt_ind.list()])#.numerator()
            elt_red *= red_p**(-1)
            if ind == ind_to_zero:
                correction_term = elt_red
            else:
                local_sols.append(elt_red)

                    
            normalizer = correction_term*val_ext[ind].uniformizer()            
#            red = p
#            if ram_ind[ind] > 1:
#            	red *= p
#            tmp = reduce()	 
            normilized_uniformizer.append(sum(local_sols) + normalizer)
    return normilized_uniformizer
                


    
############## fractional ideals #######


class dedekind_ideal(object):

    def __init__(self, num_ideal, den, is_integral, factors):
        self.num_ideal = num_ideal
        self.den = den
        self.is_integral = is_integral
        self.factors = factors
        
    def __repr__(self):
        return str(self.num_ideal)+"/" +str(self.den)
    
	#def __invert__(self)
	
    def inver(self):
        num_ideal = self.num_ideal.inver()

        den = gcd(self.den,num_ideal.den)

        num_ideal.num_ideal.min_gen *= den

        num_ideal.num_ideal.gen *= den		
        num_ideal.num_ideal.norm *= den**num_ideal.num_ideal.gen.parent().degree()				
        num_ideal.den /= den 

        return num_ideal
		
	
#	def __eq__(self,I):
    
    def __mul__(self,I):
        ''' 
        We assume that maximal order = equation order
        Multiplies fractional ideals 
        '''
        # TODO catch multiplication by trivial ideal
        # TODO invert prod_gen wrt maximal order
        print"a",self.num_ideal,I.num_ideal
        num_ideal = self.num_ideal*I.num_ideal
        print"b"
        den = self.den*I.den
        #tmp = denominator(num_ideal.gen^-1)
        d = gcd([i.numerator() for i in numerator(num_ideal.gen).list()]) # TODO should be a better way
        print"d =",d,den,num_ideal.min_gen
        num_ideal.min_gen /= d
        num_ideal.gen /= d
        den /= d
        print"sad"
        num_ideal.norm /= d**num_ideal.gen.parent().degree()
        is_integral = False
        print"sad 1"
        if den.is_unit():
        	is_integral = True
        print"sad2"        
        return dedekind_ideal(num_ideal, den, is_integral, None)

#    def __div__(self,a):
#		"""
#		division of fractional ideals
#		"""
	
	
#        self.factors = factors

    
    
    
################## Testing ################    
'''

K.<a> = NumberField((x^3+3)*(x^3+6)+81)

Zx.<x>=PolynomialRing(ZZ)

Zx.<x>=PolynomialRing(ZZ)
p = 3

P = decompose_prime(K,p)[0]

p = 7

Q = decompose_prime(K,p)[0]

R = P*Q

p = 11

#L = decompose_prime(K,p)[0]

W = dedekind_ideal(P,1,false,None)

WW = dedekind_ideal(Q,5,false,None)

U = W*WW

PP = P.inver()

O = PP.inver()

Z = PP*W


k.<t> = FunctionField(GF(11))
kx.<x> = k[]
f = x^3+t^2*x+8
F.<theta> = k.extension(f)


p = t
P = decompose_prime(F,p)[0]

print"two_elt_test t", P.val(P.gen)

q = t+1

Q = decompose_prime(F,q)[1]

print"two_elt_test t+1", Q.val(Q.gen)


R = P*Q

print"after product",[P.val(R.gen),Q.val(R.gen)],R.min_gen
'''

'''
####### Test

K = F


p = t 
v = Valuation(K,p)
val_ext_t = v.extensions(K)
normalized_uniformizer_t = all_normalized_uniformizer(val_ext_t)
[w._base_valuation._improve_approximation() for w in val_ext_t]
wt,wwt = val_ext_t
print"test t",wt(normalized_uniformizer_t[0]),wwt(normalized_uniformizer_t[1])




p = t+1 
v = Valuation(K,p)
val_ext = v.extensions(K)
normalized_uniformizer = all_normalized_uniformizer(val_ext)
[w._base_valuation._improve_approximation() for w in val_ext]
w,ww = val_ext
print"test t+1",w(normalized_uniformizer[0]),ww(normalized_uniformizer[1])

'''


