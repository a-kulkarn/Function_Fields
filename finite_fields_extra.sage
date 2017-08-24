##### This file extands finite fields with some extra structure to built towers of finite
##### fields and to comunicate within relative extensions


########
#Helpfunctions
########

def split(a, n):
    k, m = divmod(len(a), n)
    return list((a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n)))


########

def Embed(self, alpha):
	"""
    Let F/K be a finite extension of GF
	returns embedding of alpha into F
	"""
	into_error(EE,[self,alpha])
	if alpha in ZZ:
		return self.ext_field(ZZ(alpha))
	if alpha.parent().is_prime_field():
		tmp = alpha.polynomial().list()
	else:
		tmp = vector(alpha)
	return sum([tmp[i]*self.gamma[i] for i in range(len(tmp))])



class GF_Extension(object):



    def __init__(self, F, K, pol):
    	self.ext_field = F
        self.base_field = K        
        self.res_field = K.extension(pol, 'xx') 
    	self.pol = pol
#    	print"huch 0"
        self.index = ZZ(F.degree()/K.degree())
    	self.prime_field = false
#    	print"huch0.1"
    	R.<x> = F[]	
    	mipo = 	R(K.0.minimal_polynomial())
    	gamma = mipo.roots()[0][0]
    	self.root = gamma
    	self.gamma = [gamma^exp for exp in range(K.degree())]

  #  	print"huch", len(self.gamma)
    	if F.is_prime_field():
   			M = matrix([(F.0^i*K.0^j).polynomial().list() for i in range(ZZ(F.degree()/K.degree())) for j in range(K.degree())])
    	else:
   			
	    	M = matrix([vector(F.0^i*Embed(self,K.0^j)) for i in range(ZZ(F.degree()/K.degree())) for j in range(K.degree())])
  #  	print"huch 2"
    	self.M = M^-1	
  #  	print"Mu"	
    	Ft.<t> = F[]
  #  	print"Mu1",pol
    	if F.is_prime_field():
    		f = Ft([i for i in pol.list()])
    	else:    	
	    	f = Ft([Embed(self,i) for i in pol.list()])
  #  	print"Mu2",f
    	z_g = f.roots()[0][0]
   # 	print"Mu3"
    	self.pol_root = z_g	
  #  	print"Mu4"
    	r = pol.degree()
    	
    	if F.is_prime_field():
    		M_1 = matrix([z_g^i for i in range(r)])
    		M_2 = matrix([F.0^i for i in range(r)])     	
    	else:
    		M_1 = matrix([elt_seq_ini(self, z_g^i) for i in range(r)])
    		M_2 = matrix([elt_seq_ini(self, F.0^i) for i in range(r)])     	
    #	print"Mu5"    		
    	self.trans_M = M_2*M_1^-1
    	self.trans_M_inv = self.trans_M^-1

    	## I want root of pol inside F
    #	print"Mu6"    	
    	pol_emped = R([Embed(self,co) for co in pol.list()])
    	
    #	print"Mu7"    	
    	self.root_pol = pol_emped.roots()[0][0]
    	
		


def elt_seq_ini(self, alpha):
#		"""
#		alpha has to be elt in F
#		returns coefficients of alpha over K
#		"""
	if not alpha in self.ext_field:
		raise Exception, "Error: alpha must be in self.ext_field."    	
	
	if self.ext_field.is_prime_field():
		return [alpha*self.M[0,0]]
	
	sol = vector(alpha)*self.M
	coeffs = split(sol,self.index)
	K = self.base_field
	K_basis = [K.0^j for j in range(K.degree())]
	return [sum([ coeff[j]*K_basis[j] for j in range(len(coeff))])   for coeff in coeffs]
	
def elt_seq(self, alpha):
#		"""
#		alpha has to be elt in F
#		returns coefficients of alpha over K
#		"""
	if not alpha in self.ext_field:
		raise Exception, "Error: alpha must be in self.ext_field."    	
	if self.ext_field.is_prime_field():# or len(self.ext_field) == len(self.base_field):
	#	print"bin im fall drin",alpha*self.M[0,0]
	#	into_error(EE,[self, alpha])
	#	return
		return [alpha*self.M[0,0]]
	if self.ext_field.order() == self.base_field.order():
		return elt_seq_ini(self, alpha)
		
	
	
	sol = vector(alpha)*self.M#*self.trans_M
	coeffs = split(sol,self.index)
	K = self.base_field
	K_basis = [K.0^j for j in range(K.degree())]
#	return [sum([ coeff[j]*K_basis[j] for j in range(len(coeff))])   for coeff in coeffs]	

	tmp = [sum([ coeff[j]*K_basis[j] for j in range(len(coeff))])   for coeff in coeffs]	
	return (vector(tmp)*self.trans_M).list()


def lift(self, alpha):
#		"""
#		alpha has to be elt in F
#		returns coefficients of alpha over K
#		"""
	if not alpha in self.ext_field:
		raise Exception, "Error: alpha must be in self.ext_field."    	
	
	coeffs = elt_seq(self,alpha)
	var = self.pol.parent().0
	return sum([ coeffs[j]*var^j for j in range(len(coeffs))])




	
def iso_to_gf(self,alpha):
	"""
	Input: an element alpha in R = self.base_field/self.pol (= K[t]/<g(t)>)
	Output: alpha embedded into self.ext_field = F
	"""
    
	coeffs = alpha.list()
	F_rel_basis = [self.pol_root^i for i in range(self.index)]
	return sum([coeffs[i]*F_rel_basis[i] for i in range(len(coeffs))])
	
	
	
	
def iso_to_res_field(self,alpha):
	"""
	Input: an element alpha in self.ext_field = F
	Output: alpha embedded into R = self.base_field/self.pol (= K[t]/<g(t)>)
	"""
    
    
	coeffs = elt_seq(self,alpha)
	return L.res_field((vector(coeffs)*L.trans_M).list())
	
		
	
class GF_Tower(object):


	def __init__(self,L_ext):
		self.prime_field = L_ext.base_field.base_ring()
		self.top_field = L_ext.ext_field
		self.fields = [L_ext]
		self.pols = []
		if not L_ext.base_field.is_prime_field():
			self.pols = [L_ext.base_field.0.minimal_polynomial()]
		self.pols.append(L_ext.pol)
		if L_ext.base_field.is_prime_field():
			self.basis = [L_ext.root^i for i in range(L_ext.index)]
			M = matrix([vector(self.basis) for bas in self.basis])
			self.transition_matrix = M^-1
		else:
			#root = Embed(L_ext,L_ext.base_field.0)
			#coefficient_bas = [root^exp for exp in range(L_ext.base_field.degree()) ]
			self.basis = [L_ext.ext_field.0^i*coeff for coeff in L_ext.gamma for i in range(L_ext.index) ]
			M = matrix([vector(bas) for bas in self.basis])
			self.transition_matrix = [M^-1]
			
			
def add_new_level(self, L_ext):
	self.fields.append(L_ext)  
	self.pols.append(L_ext.pol)
	bas_old =  [Embed(L_ext,i) for i in self.basis]	
	self.basis = [L_ext.ext_field.0^i*coeff for coeff in bas_old for i in range(L_ext.index) ]		
	M = matrix([vector(bas) for bas in self.basis])
	self.transition_matrix.append(M^-1) 
		
			
def change_rep(self,alpha):
	"""
	Input: an element alpha in finite field F/F_p
	Output: alpha in F with basis coming from GF_Tower
	"""
    
    
	coeffs = vector(alpha)*self.transition_matrix[-1]
	return sum([coeffs[i]*self.basis[i] for i in range(len(coeffs))])
''' 
def lift(self,alpha):
	"""
	Input: an element alpha in finite field F/F_p
	Output: alpha in F with basis coming from GF_Tower
	"""
    
    
	coeffs = vector(alpha)*self.transition_matrix
	return sum([coeffs[i]*self.basis[i] for i in range(len(coeffs))])

				
			
   	
	def add_new_level(self, F, K, pol):
		L = GF_Extension(F, K, pol)
		self.fields.append(L)    
'''
