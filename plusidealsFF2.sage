# -*- coding: utf-8 -*-

import pprint
from itertools import product
pp = pprint.PrettyPrinter(indent=4)




def is_prime_polynomial(g):
	"""
	Returns True if g is a prime polynomial, and False otherwise
	"""
	fac=factor(g)
	if len(list(fac)) == 1 and fac[0][1] == 1:
		return True
	else:
		return False	




def ab_slope(a, b, cloud):
    return QQ(cloud[b][1]-cloud[a][1])/QQ(cloud[b][0]-cloud[a][0])

def lower_convex_hull(cloud):
    """
    Computes the lower convex hull of a cloud of points.
    """
    # Special case where our cloud is a single point.
    if len(cloud) == 1:
        return [ Side(0, list2point(cloud[0]), list2point(cloud[0])) ]

    slopes = [ [0, 1, ab_slope(0, 1, cloud)] ]
    b = 2
    while b < len(cloud):
        for i in [0..len(slopes)-1]:
            a = slopes[i][0]
            slope = ab_slope(a, b, cloud)
            if slope <= slopes[i][2]:
                slopes = slopes[0:i]
                slopes.append([a, b, slope])
                break
        if slopes[-1][1] != b:
            a = b - 1
            slope = ab_slope(a, b, cloud)
            slopes.append([a, b, slope])
        b += 1

    slopes = [ Side(slope, list2point(cloud[a]), list2point(cloud[b])) for a, b, slope in slopes ]
    return slopes

def list2point(l):
    return Point(l[0], l[1])

class Point(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __unicode__(self):
        return '(%d, %d)' % (self.x, self.y,)

    def __repr__(self):
        return self.__unicode__()

    def list(self):
        return [self.x, self.y]

class Side(object):

    def __init__(self, slope, p1, p2):
        self.slope = slope
        self.p1 = p1
        self.p2 = p2

    def __unicode__(self):
        return '[%s, %s, %s]' % (unicode(self.slope), unicode(self.p1), unicode(self.p2),)

    def __repr__(self):
        return self.__unicode__()

    def width(self):
        return self.p2.x - self.p1.x 
    
    def height(self):
        return self.p1.y - self.p2.y

class MontesType(object):

    def __init__(self, F, p, varphi, omega, initial=True):
        self.parent = F
        self.pol = minimal_polynomial(F.0)
        self.prime_polynomial = p
        self.varphi = varphi        
        self.levels = [ ]
        self.sfl = [0, 0, 0, 0] # Single Factor Lifting
        self.phiadic = [F.base_field()(0) for i in range(0, 4)]	 
#        self.res_class_field = varphi.parent().base_ring()


        if initial is True:
            self.add_initial_level(varphi, omega, p)
    
    def copy(self):
        copy = MontesType(self.parent, self.prime_polynomial, self.varphi,
                          self.levels[0].omega, initial=False)

        for level in self.levels:
            copy.levels.append(level.copy())
        
        return copy

    def add_initial_level(self, varphi, omega, p):
        print"vorher"	        
        phi = self.pol.parent()(varphi)		#sollte richtig sein
	
        new_level = MontesTypeLevel(phi, omega, p)
 #       q = len(p.parent().base_ring())^(p.degree()*varphi.degree())
#        print"mal sehen inadd_ini "  ,q		
	
        new_level.Fq = list(varphi)[0].parent()

        print"nachher"
        into_error(E,[self, varphi, omega, p])                           
        print"waaaas"
        if varphi.degree() > 1:
        	new_level.z = new_level.Fq.0                                    
        else:
        	tmp = varphi.coeffs()[0] 	        	
        	new_level.z = -tmp.parent().base_ring()(tmp)
        new_level.Fqy = PolynomialRing(new_level.Fq, 'y0')
        new_level.inv_M = matrix([])

        self.levels.append(new_level)

    def add_new_level(self, phi, omega, u):
        s = len(self.levels)
        lvl_s = self.levels[-1]
        new_level = MontesTypeLevel(phi, omega, self.prime_polynomial)
        new_level.V = lvl_s.e*u
        new_level.prod_e = lvl_s.prod_e * lvl_s.e
        new_level.prod_f = lvl_s.prod_f * lvl_s.f        
        print"in Ad_new 1"        
#        q=(len(self.prime.parent().base_ring())^self.prime.degree())^new_level.prod_f
        print"in Ad_new 2"
        into_error(E,[self,lvl_s.res_pol])      
        #new_level.Fq = ResidueField(lvl_s.res_pol, name='w'+str(s),
        #                           conway=True, prefix='ww')
		#       print"mal sehen "  ,new_level.Fq 
        new_level.Fq = lvl_s.Fq.extension(lvl_s.res_pol, name='w'+str(s), prefix='ww') #conway=True,
        print"in Ad_new 3"
        new_level.Fqy = PolynomialRing(new_level.Fq, 'y'+str(s))
        #new_level.embedding = Hom(lvl_s.Fq, new_level.Fq).list()[0]
         
        if lvl_s.f > 1:
            
            lifted_res_pol = new_level.Fqy(
                    [new_level.Fq(c) for c in lvl_s.res_pol])

            # Take the first root of the residual polynomial in the new
            # level's Fq[y].
            new_level.z = lifted_res_pol.roots()[0][0]
            M = matrix([vector(new_level.z^j * lvl_s.Fq.gen()^k)
                        for j in range(lvl_s.f)
                        for k in range(lvl_s.prod_f)])
            new_level.inv_M = M^(-1)
        else:
            new_level.z = -list(lvl_s.res_pol)[0]
            new_level.inv_M = lvl_s.inv_M

        self.levels.append(new_level)
        return new_level

    def refinement(self, more_factors=False):
        r = len(self.levels) - 1
        if more_factors is True:
            self.lvl(r).refinements.append([self.lvl(r).phi,
                                            self.lvl(r).slope])
        
        #self.lvl(r).cutting_slope = ZZ(self.lvl(r+1).slope)
        self.lvl(r).phi = self.lvl(r+1).phi
        self.lvl(r).omega = self.lvl(r+1).omega

        self.remove_last_level()

    def add_last_level(self, phiadic, side, side_dev, last_psi=True):
        #tt.add_last_level(phiadic, sides[0], side_devs)

        r = len(self.levels)
        lvl_r = self.lvl(r)
        lvl_r.e = 1
        if r > 1:
            # FIXME: What does nur stand for?
            nur = sum([self.lvl(j).slope/self.lvl(j).prod_e for j in [1..r-1]])
            self.sfl[0] = floor((lvl_r.V/lvl_r.prod_e)-nur)

        if side.p1.x == 0:
            print"sloppy in add_last_level",side,side.slope        
            slope = -side.slope
            lvl_r.h = ZZ(slope)
            # FIXME: Why do we just take the first 2 elements?
          
            self.phiadic[0:2] = phiadic[0:2]
            if last_psi:
  #              print"inside2",r, side_dev
                res_pol = self.residual_polynomial(r, side_dev)
   #             print "last monic res.pol. = %s (%s)" % (res_pol, res_pol.monic())

                lvl_r.res_pol = res_pol.monic()

                lvl_r.log_gamma = lvl_r.log_phi - (lvl_r.h * lvl_r.log_pi)
  
        else:
            slope = +Infinity
        lvl_r.slope = slope
        
        return lvl_r

    def remove_last_level(self):
        self.levels.pop()

    def lvl(self, i):
        ''' 1 indexed access to array levels.'''
        assert i > 0, 'tt.lvl(i) is 1 indexed, so i must be > 0.'
        return self.levels[i-1]

    def __unicode__(self):
        return '(%s; %s)' % (self.varphi,
                             '; '.join([unicode(l) for l in self.levels]))

    def __repr__(self):
        return self.__unicode__()

    def rth_level(self):
        return self.levels[-1]


    #### Begin serious Montes methods ####

    ## Phi-Newton Polygons ##
    def phi_newton_polygon(self, i, phiadic):
        assert i <= len(self.levels)

        n = 0
        sides = []
        cloud = []
        side_devs = []
        all_devs = []
        for k in [0..len(phiadic)-1]:
            val = 0
            dev = [ ]

            val, dev = self.value(i, phiadic[k])
            #print "N_%d : v_%d(%s phi_%d) = %s + %d" % (i, i, phiadic[k], i, str(val), n)
            if abs(val) != Infinity:
                cloud.append([k, val + n])
                all_devs.append(dev)
            n += self.lvl(i).V

    #    print "%d. cloud: %s" % (i, cloud,)
        sides = lower_convex_hull(cloud)
        abscissas = [ p[0] for p in cloud ]    

        for side in sides:
            height = ZZ(side.p1.y) # [1][1]
            dev = [ ]

            # [1][0], [2][0]
            for j in range(side.p1.x, side.p2.x + side.slope.denominator(), side.slope.denominator()):
                try:
                    position = abscissas.index(j)
                except ValueError:
                   position = -1 
                
                if position > -1 and cloud[position][1] == height:
                    dev.append(all_devs[position])
                elif i == 1:
                    dev.append(0)
                else:
                    dev.append([])
                height += side.slope.numerator()
                
            dev.append(side.p1.list()) # list(side[1])
            side_devs.append(dev)    

        if len(sides) == 0:
            raise NotImplementedError, "Case for a zero side Newton Polygon has not been implemented"


        return sides, side_devs

    ## Value of polynomial a(x) in ZZ[x] ##
    
    def value(self, i, a):
        assert i <= len(self.levels)+1
        print"bin in value"
        val = +Infinity
        if a == 0:
            if i == 1:
                devs = 0
            else:
                devs = []
            return val, devs

        if i == 1:
            val = min([valuation(c, self.prime_polynomial) for c in list(a)])
            devs = a
        else:
            devs = []
            lvl_im1 = self.lvl(i-1)
            step = lvl_im1.V + lvl_im1.slope
            min_height = 0
            V_height = 0
            quot = a
            k = 0
            last = 0
            
            print"what is"
            while quot != 0 and min_height <= val: 
            	# muss eigentlich <= heißen
                quot, ak = quot.quo_rem(lvl_im1.phi)
                new_val, dev = self.value(i-1, ak)
                candidate = new_val + min_height
                if candidate <= val:
                    if candidate < val:
                        val = candidate
                        first_abscissa = k
                        first_ordinate = new_val + V_height
                        devs = [ dev ]
                    else:
                        for j in range(last+lvl_im1.e*2, k+1, lvl_im1.e):
                            if i-1 == 1:
                                devs.append(0)
                            else:
                                devs.append([])
                        devs.append(dev)
                    last = k
                min_height += step
                V_height += lvl_im1.V
                k += 1
            # FIXME: Why are we so certain that first_abscissa and
            #        first_ordinate will have been set?
            devs.append([first_abscissa, first_ordinate])
            val = ZZ(lvl_im1.e * val)

        return val, devs

    ## Residual Polynomial ##
    def residual_polynomial(self, i, side_devs):
        """
        Creates and returns the i-th residual polynomial of a polynomial (f).
        side_devs is a list of the multiadic expansions of the coefficients of
        f whose attached points in N_r(f) lie on teh side S of slope
        -levels[i].slope. The last element of side_devs is [s, u], where s,u
        are the coordinates of the left end-point on S.
        """
        into_error(E,[self,i, side_devs])    
        print"\n\n\n in residual_polynomial"        
        assert i <= len(self.levels)
        print"was issen? 1",i,side_devs
        lvl_i = self.lvl(i)
       	p=self.lvl(1).prime_polynomial.numerator()
        height = side_devs[-1][1]
        #print "height, side devs:", height, side_devs
        res_coeffs = [ ]
 #    	print"was is das",	i
        print"was issen? 2 lalei"
        for dev in side_devs[:-1]:

            if (i == 1 and dev == 0) or (i > 1 and len(dev) == 0): 
                print"was issen? 2.1"
                res_coeffs.append(lvl_i.Fq(0))
            elif i == 1:
                print"was issen? 2.2"
                # coefficients are polynomials too."
#		print"in new0000"
                coeff = dev // self.prime_polynomial^(height)
	
                k=self.lvl(1).Fq
                

                Fqy=self.lvl(1).Fqy
#                 
#                 if p.degree() > 1:
# 	           
# 	                coeff=Fqy([k(list(ii.numerator() % p.numerator())) for ii in list(coeff)]) 
#                             
#                 else:
# 	                coeff=Fqy([k(ii.numerator() % p) for ii in list(coeff)])                		                
                coeff = Fqy([ k.reduction_map()(ii.numerator()) for ii in list(coeff)]) 
                   		                
 #           		print"in res pol for", coeff
                print"vor eval"  
                into_error(E,[lvl_i,coeff])               
                res_coeffs.append(coeff(lvl_i.z))
                print"nach eval"
#                print"\n\n daten", dev, coeff,	coeff(lvl_i.z), lvl_i.z, lvl_i.z.parent()
                j = side_devs.index(dev)
                #print "%d. order Res.Pol. the easy way c_%d = (%d / %d) = %d (z_%d = %d)" % (i, j, dev, self.prime**height, dev // self.prime^(height), i-1, lvl_i.z)
            else:
            	print"was issen? 3"		
                lvl_im1 = self.lvl(i-1)
            	print"was issen? 3.1",	  height.parent() , height
                twist_exp = (dev[-1][0] - lvl_im1.inv_h*height) // lvl_im1.e
                #print "%d. order Res.Pol. twist exp: %s" % (i, twist_exp,)
            	print"\n\n in put respol",i-1, dev			
                # coefficients are polynomials too.
            
                coeff = self.residual_polynomial(i-1, dev)
                print"\n\n respol rekursive", coeff                                
                lifted_coeff = lvl_i.Fqy([lvl_i.Fq(c) for c in coeff])
        #        print "------"
         #       print "coeff:", coeff, "parent:", coeff.parent()
          #      print "z:", lvl_i.z, "parent:", lvl_i.z.parent()
           #     print "------"
                res_coeffs.append(lvl_i.z^(twist_exp) * lifted_coeff(lvl_i.z))
            	print"was issen? 5"	
                j = side_devs.index(dev)
                #print "%d. order Res.Pol. c_%d = %s (z_%d = %s)" % (
                #        i, j, lvl_i.z^(twist_exp)*coeff(lvl_i.z), i, lvl_i.z,)
            	print"was issen? 6"	

            height = height - lvl_i.h
#        print"\n \n\n in res pol with coeff", res_coeffs 
        res_pol = lvl_i.Fqy(res_coeffs)
 #       print "%d. order Res.Pol = %s" % (i, res_pol,)
        print"END residual_pol!!!!"
        return res_pol

    ## Representative ##
    def representative(self, omega):
        """
        Construct a representative phi of a type. A new level is added with
        phi and V.
        """

        s = len(self.levels)
        lvl_s = self.lvl(s)
        ef = lvl_s.e * lvl_s.f
        u = ef * lvl_s.V+lvl_s.f*lvl_s.h

        if s > 1:
            txp = -self.lvl(s-1).inv_h * (u // self.lvl(s-1).e)
            twist = lvl_s.z^(txp)
        else:
            twist = lvl_s.Fq(1)

        # Reductum from magma
        res_pol = twist * lvl_s.Fqy(list(lvl_s.res_pol)[:-1])
  #      print "%s = %s * Reductum(%s) = %s * %s" % (res_pol, twist, lvl_s.res_pol, twist, lvl_s.Fqy(list(lvl_s.res_pol)[:-1]))
     #   u += lvl_s.f * lvl_s.h
        phi0 = self.construct(s, res_pol, 0, u)
   #     print "%d. Repr. = %s + phi^ef" % (s, phi0,)

        phi0 = phi0 + lvl_s.phi^(ef)
        print"wo is der fehler etwa in add_new"        
        new_level = self.add_new_level(phi0, omega, u)
        print"doch nicht in add_new"        
        return new_level

    ## Construct ##
    def construct(self, i, res_pol, s, u):
        """
        This routine constructs a polynomial phi0 with integer coefficients
        such that:
          - deg phi0 < m_i + 1 and y^nu*R_i(phi0)(y) = res_pol(y), where
            nu = ord_y(res_pol).
          - The non-negative integers s,u are the coordinates of teh left
            endpoint of a segment of slope -lvl_i.slope supporting
            N_i(phi0)
        """
        assert i <= len(self.levels), "i (%d) must be <= #levels (%d)" % (i, len(self.levels),)
        lvl_i = self.lvl(i)
        assert res_pol.degree() < lvl_i.f, "res_pol is too large."
        assert u + s*lvl_i.slope >= lvl_i.f*(lvl_i.e*lvl_i.V + lvl_i.h), "the point (s, u) is too low"
        print"in construct"
        var = lvl_i.phi^(lvl_i.e)
        phi0 = 0
        height = u - res_pol.degree()*lvl_i.h
#        print "%d. Construct res.pol. coeffs: %s" % (i, unicode(list(reversed(list(res_pol)))),)
        print"in construct step 1"
        if i == 1:
            for a in reversed(list(res_pol)):
                print"in construct step 1.1"
                A= self.pol.coefficients()[0].numerator().parent()
                lift = A(list(a.polynomial()))

                phi0 = (phi0 * var) + (lift * self.prime_polynomial^height)
                height = height + lvl_i.h
        else:
            step = (lvl_i.e * lvl_i.V) + lvl_i.h
            new_V = u - (res_pol.degree() * step) - (s * lvl_i.V)
            lvl_im1 = self.lvl(i-1)
            print"in construct step 2"
            for a in reversed(list(res_pol)):
                # FIXME: What does pj stand for?
                pj = 0
                if a != 0:
                    txp, s_im1 = divmod(lvl_im1.inv_h*height, lvl_im1.e)
                    u_im1 = (new_V - (s_im1 * lvl_im1.h)) // lvl_im1.e
                    c = (a*lvl_i.z^txp)
                    # Doing the equivalent of (if Eltseq existed)
                    #         lvl_im1.Fqy(Eltseq(c, lvl_im1.Fq))
                    if c.parent().base_ring() == lvl_im1.Fq:
                        eltseq = list(c.polynomial())
                    else:
                        # FIXME: It may not be the most efficient way of doing
                        # this but it's pretty efficient for now. It is based
                        # on this answer:
                        # http://ask.sagemath.org/question/3398/representing-finite-field-elements-in-terms-of?answer=4540#4540
                        eltseq = None
                        w_im1 = lvl_im1.Fq.gen()
                        prod_f = lvl_im1.prod_f
                        m = vector(c) * lvl_i.inv_M
                        eltseq = [lvl_im1.Fq(m[j*prod_f:(j+1)*prod_f])
                                   for j in range(len(m)/prod_f)]
                        if c != sum([eltseq[j]*lvl_i.z^j for j in range(len(eltseq))]):
                            raise(Exception,
                                 "Eltseq calculated incorrectly for {0}"\
                                         .format(c))

       #             print "%d.   Eltseq(%s, %s) = %s" % (i, c, lvl_im1.Fq, eltseq)
                    new_res_pol = lvl_im1.Fqy(eltseq)
                    pj = self.construct(i-1, new_res_pol, s_im1, u_im1)
                phi0 = (phi0 * var) + pj
                new_V += step
                height += lvl_i.h
        print"in construct step rather End"
        phi0 = phi0 * lvl_i.phi^(s)
#        print "%d. Construct pol = %s" % (i, str(phi0))
        print"in construct step THE End"
        return phi0

    def single_factor_lifting(self, slope):
        """
        Perform single factor lifting on the representative of this type. The
        aim is to make self.rth_level().slope >= slope.
        """
        r = len(self.levels)
        lvl_r = self.rth_level()
        
        if lvl_r.slope >= slope:
            return
        if self.sfl[2] == 0:
            self.sfl_init()

        print "SFL: {0} (slope: {1})".format(self.sfl, self.rth_level().slope)
		
        pmode = 'terse'
        p = self.prime_polynomial
        exponent = self.sfl[0]
        nu = self.sfl[1]
        x0prec = self.sfl[2]
        x0num = self.phiadic[3]
        x0den = self.sfl[3]
        print"Test 1" 
        prod_e = lvl_r.prod_e
        h = lvl_r.h - lvl_r.cutting_slope
        print"cutting slope "  ,lvl_r.cutting_slope
        print"h= "  ,lvl_r.h
        last_h = slope - lvl_r.cutting_slope
        V = lvl_r.V + lvl_r.cutting_slope

        precision=nu + exponent + ceil((V+last_h)/prod_e)
        p_prec = p^precision
        print"Test 2"         
        pol_zp = coeff_mod(self.pol, p_prec)
        psinum_zp = coeff_mod(self.phiadic[2] , p_prec) #PsinumZp:=type[1]`Phiadic[3] mod p_prec;
        print"Test 3", last_h, h 
        print"x0prec", x0prec         
        path = path_of_precision(last_h, h)
        short_path = path_of_precision(h, x0prec)
        print"also hier 1"   
        newprecision= nu + exponent + ceil(h/prod_e) #Ceiling((V+path[2])/e)
        p_prec_new=p^newprecision;
        a1=coeff_mod(self.phiadic[1] , p_prec_new)		#type[1]`Phiadic[2] mod p_prec_new;
		#print "That precision:", nu + exponent + ceil(h/prod_e);
        print"also hier 2"        
        newprecision=nu + exponent + ceil((V+path[1])/prod_e)
        p_prec_new=p^newprecision
        phi=coeff_mod(lvl_r.phi , p_prec_new)
        psinum= coeff_mod(psinum_zp , p_prec_new )
        print"also hier 3", ((coeff_mod(self.phiadic[0] , p_prec_new ))*psinum) , phi
        
        print'ok' 


        a0num, a0den = self.cancel((coeff_mod(self.phiadic[0] , p_prec_new )*psinum) % phi, nu)
        print "a0 num, den:", [a0num, a0den]
        a1num, a1den = self.cancel((coeff_mod(a1 , p_prec_new)*psinum) % phi, nu)
        print "a1 num, den:", [a1num, a1den]

        print "--==--==--==--==--==--==--"

        for i in range(1, len(short_path)):
            low_precision = a1den + 2 * x0den  + ceil(short_path[i]/prod_e)
            print"waddehaddedudeda"
            x0num, x0den = self.inversion_loop([ a1num, a1den], x0num, x0den,
                                               phi, low_precision,p)
            print "x0 num, den:", x0num, x0den
            print '----'

        print "--==--==--==--==--==--==--"

        anum, aden = self.cancel((a0num*coeff_mod(x0num ,  p_prec_new)) % phi, x0den+a0den)
        phi = phi + anum

        print "between:", [anum, aden, phi]
        print "--==--==--==--==--==--==--"

        for i in range(1, len(path)-1):
            loop_prec = nu + exponent + ceil((V+path[i+1])/prod_e)
            
            ploop=p^loop_prec
            phi=coeff_mod(phi , ploop)
            Psinum = coeff_mod(psinum_zp , ploop)

            qq, c0 = pol_zp.quo_rem(phi)

            c1 = coeff_mod(qq , phi)

            c0num, c0den = self.cancel(coeff_mod(c0*psinum, ploop) % phi, nu)
            c1num, c1den = self.cancel(coeff_mod(c1*psinum, ploop) % phi, nu)
            #print "cs:", [c1, c0num, c0den, c1num, c1den]





            low_precision = c1den + 2 * x0den + ceil(path[i]/prod_e)
            x0num, x0den = self.inversion_loop([c1num, c1den], x0num, x0den,
                                               phi, low_precision,p)
            #print "x0 num, den", x0num, x0den

            xnum =x0num % ploop
            cnum, cden = self.cancel((c0num * coeff_mod(x0num , ploop)) % phi, x0den)
            phi = coeff_mod((phi + cnum) , ploop)

            print "phi:", phi
            print "--------------------------"


        print "--==--==--==--==--==--==--"

        self.sfl[2] = max(h, path[-2])
        lvl_r.phi = coeff_mod(phi,p_prec) #lvl_r.phi.parent()(phi)
        self.phiadic[3] = x0num
        self.sfl[3] = x0den
    
    
    def single_factor_lifting_Up(self, slope):
    	print"in SFL UPPI"    
    	if self.rth_level().slope >= slope:
    		return
    	print"in SFL UPPI 2"
    	self.single_factor_lifting(slope)
    	print"\n\n \n\n in SFL UPPI 2\n\n"
    	self.update_last_level()
    	
            
    def sfl_init(self):
        p = self.prime_polynomial.numerator()
        lvl_r = self.rth_level()
        prod_e = lvl_r.prod_e
        a1 = self.phiadic[1] # FIXME: Why [1]?

        A=p.parent()

        Ax.<x>=PolynomialRing(A)
        psinum = Ax(1)

        print"sdasd2"
        
        r = len(self.levels) - 1

        if r == 0:
            nu = min([valuation(a.numerator(), p) for a in a1])
            a1	=a1 // p^nu 
        
            F_0=FiniteField(len(A.base_ring())^p.degree(),conway=True,prefix='zini')

            F_0_y0=PolynomialRing(F_0, 'y0')

            if p.degree() > 1:
            	res_pol=F_0_y0([F_0(list(i.numerator() % p)) for i in list(a1)])
            else:
            	res_pol=F_0_y0([F_0(i.numerator() % p) for i in list(a1)])
 
            # Evaluate a1/p^nu in z_1 (this may be z_0)
            klass = res_pol(self.lvl(1).z)  
            print"klass",klass          
        else:
            val, dev = self.value(r+1, a1)
            res_pol = self.residual_polynomial(r, dev)

            logpsi = 0
            qq, s = (-val).quo_rem(prod_e)
            psinum, logpsi = self.prescribed_value(s)
            print'psinum',psinum.parent()
            nu = -logpsi[0] - qq
            vector = dev[-1][0] * self.lvl(r).log_phi + dev[-1][1]*self.lvl(r).log_pi
            klass = self.convert_logs(logpsi + vector)
            klass *= res_pol(self.lvl(r+1).z)

        self.phiadic[2] = psinum
        self.sfl[1] = nu
        self.sfl[2] = 1

        x0num = 0
        x0den = 0
        print"before classes",r,klass.parent()            
        x0num, x0den = self.local_lift(klass^(-1))
        print"before 2",r
        self.phiadic[3] = x0num
        self.sfl[3] = x0den

    def prescribed_value(self, value):
        """
        From +Ideals:
        If we are attached to the prime ideal P with Okutsu depth r, then
        logpsi=[a_0, ..., a_r] and psi=phi_1^a_1 ... phi_r^a_r, with
        v_P(p^a_0 psi(theta))=value.
        """
        print"\n\n parent hierrrr", self.lvl(1).phi.parent(),
        print"\n\n"       
        Ax=self.lvl(1).phi.parent()
        print"Ax",Ax
        psi = Ax(1)
        r = len(self.levels)
        logpsi = (ZZ^r)(0)
        qq, val = value.quo_rem(self.rth_level().prod_e)
        logpsi[0] = qq
        if val > 0:
            body = val
            for k in reverse(range(r-2)):
                jj = (self.levels[k].inv_h * body) % self.levels[k].e
                logpsi[k+1] = jj
                psi = psi * self.level[k].phi^jj
                res = (body - jj*self.levels[k].h) // self.levels[k].e
                body = res - jj*self.levels[k].V
            logpsi[0] += res

        return psi, logpsi

    def local_lift(self, clss):
        """
        From +Ideals:
        class should belong to the residue class field  type[r]`Fq. The output
        is a pair g,e such that g(theta)/p^e is a lift to a P-integral element
        in K and deg g(x)<n_P.
        """
        i = 1
        while clss not in self.lvl(i).Fq:
            i += 1

        if i == 1:
            Ax.<x>=PolynomialRing(self.prime_polynomial.numerator().parent())
            numlift = Ax(clss.polynomial())
            denlift = 0
        else:
            raise(NotImplemented, "local lift not implemented for i > 1.")

        return numlift, denlift

    def convert_logs(self, log):
        """"
        From +Ideals:
        log[1] is not used. The product of all Phi_i^log[i] for i>0 should have
        integer value M.

        The output is the class of this product divided by p^M.
        """

        vector = log
        z = 0
        klass = self.lvl(1).Fq.prime_subfield()(1)
        for i in reversed(range(len(vector)-1)):
            ti = vector[i+1] // self.levels[i].prod_e
            z = self.levels[i].z
            klass *= z^ti
            vector = vector - ti*self.levels[i].log_gamma

        return klass

    def cancel(self, poly, den):
        if poly == 0:
            return poly, 0
        print 'In cancel'
        cancel = min([den] + [ a.valuation(self.prime_polynomial) for a in poly])
        zq = poly[0].parent()
        print "poly, p^cancel: {0}, p^{1}".format(poly, cancel)
        num = poly.parent()([ poly.parent().base_ring()(c) / self.prime_polynomial^cancel for c in poly ])
        print "num: {0}".format(num)
        #num = poly / self.prime^cancel

        return num, den-cancel

    def inversion_loop(self, A, xnum, xden, phi, precision,p):
        anum = A[0]
        aden = A[1]
        print'in inverseloop'
#        zq = Zp(self.prime, precision, type='capped-abs', print_mode='terse')
#        zqt.<t> = PolynomialRing(zq)

        p_prec=p^precision
        phip = coeff_mod(phi,p_prec) 
        xnum = coeff_mod(xnum,p_prec)
        print"bin bis hier 1"	        
        x1num, x1den = self.cancel(2*self.prime_polynomial^(xden+aden) - (coeff_mod(anum,p_prec)*xnum) % phip, xden+aden)
        print"bin bis hier 2"
        xnum, xden = self.cancel((xnum*x1num) % phip, xden+x1den)
        xnum = coeff_mod(xnum,p_prec)
        return xnum, xden
        
    def update_last_level(self):
        
        pol = self.pol
        print"\n\n\n in UpdateLasteLevel"
        r = len(self.levels)
        lvl_r = self.rth_level()    
        qq,a_0 = pol.quo_rem(lvl_r.phi)
        if a_0 == 0:
            lvl_r.slope = +Infinity
        else:
            self.phiadic[0] = a_0
            self.phiadic[1] = qq % lvl_r.phi            
            sides, side_devs = self.phi_newton_polygon(r, self.phiadic)
            lvl_r.slope = -sides[0].slope
            lvl_r.h = -ZZ(sides[0].slope)	
            tmp = self.residual_polynomial(r,side_devs[0])	
            lvl_r.res_pol = self.residual_polynomial(r,side_devs[0]).monic()
            print"what 5" , tmp
            print"sdaa",lvl_r.res_pol             
            lvl_r.log_gamma = lvl_r.log_phi-lvl_r.h*lvl_r.log_pi
			 

class MontesTypeLevel(object):

    def __init__(self, phi, omega, p):
        self.phi = phi
        self.prime_polynomial = p

        self.e = None
        self.f = None
        self.h = None
        self.prod_e = 1
        self.prod_f = 1
        self.inv_h = None

        self.V = 0
        self.omega = omega
        self.cutting_slope = 0
        self.refinements = [ ]

        self.log_pi = vector([1, 0])
        self.log_phi = vector([0, 1])
        self.log_gamma = None

        self.Fq = None
        self.z = None
        self.Fqy = None
        self.inv_M = None
        #self.embedding = None

        self.slope = None
        self.res_pol = None

    def copy(self):
        copy = MontesTypeLevel(self.phi, self.omega, self.prime_polynomial)

        copy.phi = self.phi
        copy.prime_polynomial = self.prime_polynomial

        copy.e = self.e
        copy.f = self.f
        copy.h = self.h
        copy.prod_e = self.prod_e
        copy.prod_f = self.prod_f
        copy.inv_h = self.inv_h

        copy.V = self.V
        copy.omega = self.omega
        copy.cutting_slope = self.cutting_slope
        copy.refinements = self.refinements

        copy.log_pi = self.log_pi
        copy.log_phi = self.log_phi
        copy.log_gamma = self.log_gamma

        copy.Fq = self.Fq
        copy.z = self.z
        copy.Fqy = self.Fqy
        copy.inv_M = self.inv_M
        #copy.embedding = self.embedding

        copy.slope = self.slope
        copy.res_pol = self.res_pol

        return copy

    def __unicode__(self):
        phi = self.phi
        if self.slope is not None:
            slope = self.slope
        else:
            slope = '-'
        if self.res_pol is not None:
            res_pol = self.res_pol
        else:
            res_pol = '-'

        return '(%s, %s, %s)' % (str(phi), str(slope), str(res_pol))

    def __str__(self):
        return self.__unicode__()

    def __repr__(self):
        return self.__unicode__()

def montes(F, p, basis=False):

    if is_prime_polynomial(p) is False:
        raise Exception, "Error: p must be prime."

    f =minimal_polynomial(F.0)					#K.defining_polynomial()
    if f.is_monic() is False:
    #sollte monic_integral_model verwenden um normiertes Polynom zu erzeugen
        raise Exception, "Error: def.pol. of F must be monic."
    # TODO: Add requirements that coefficients of f are integers



    reps_OM = [ ]
    trees = [ ]
    total_index=0
    p=p.numerator()
    A=p.numerator().parent()
    F_0 = ResidueField(p.numerator(), name='z0',
                                   conway=True, prefix='ww')
    #F_0=FiniteField(len(A.base_ring())^p.degree(),conway=True,prefix='zini')

    F_0_y0=PolynomialRing(F_0, 'y0')
 
    res_pol = F_0_y0([F_0.reduction_map()(co.numerator()) for co in f.coeffs()])
   
    for fac in res_pol.factor():

#	print "Analysing irred. factor modulo p: %s" % (str(factor[0]),)
        tt = MontesType(F, p, fac[0], fac[1])
        tree, tree_index = montes_main_loop(F, p, tt)
        total_index += tree_index
	
        reps_OM += tree
        for pos in [0..len(reps_OM)-1]:
        	reps_OM[pos].position = pos
        	
        trees.append(tree)

    if len(reps_OM) == 1:
        reps_OM[0].rth_level().phi = f
        reps_OM[0].rth_level().slope = +Infinity

   # print "OM Representatives:"
   # print len(reps_OM)
    #for tt in reps_OM:
    #    print "   ", tt
    
    
    	
   
    return reps_OM, total_index
    
def montes_main_loop(K, p, tt):
	    
    total_index = 0
    leaves = [ ]
    type_stack = [ tt ]
    while len(type_stack) > 0:
     #   print "STARTING LOOP:\n  type stack:", type_stack
        tt = type_stack.pop()

        r = len(tt.levels)
        lvl_r = tt.rth_level()
  
#        print "\n\n\nAnalysing type of order: %d (%s)"  % (r, unicode(lvl_r),)
       # print "--------------------------------------------------------------------------------\n"
        print"phiadic",	tt.pol
        print"phiadic2",lvl_r.phi
        print"phiadic3"	,lvl_r.omega 	
        print"\n\n"               
        phiadic, quotients = phi_expansion(tt.pol,
                                           lvl_r.phi,
                                           lvl_r.omega)

        ## Phi-Newton Polygon

        sides, side_devs = tt.phi_newton_polygon(r, phiadic)
      #  print "Sides of Newton polygon: %s" % (sides,)         
        length_N = lvl_r.omega
        index_N = -lvl_r.cutting_slope * (length_N*(length_N-1) // 2)
        starting_prod_f = lvl_r.prod_f
       # print "index N: %d, length N: %d" % (index_N, length_N)
  
  
  
        if length_N == 1:
#            print"laste leven daten",phiadic,sides[0],side_devs[0]
            tt.add_last_level(phiadic, sides[0], side_devs[0])

#            print "  #### Found a factor of depth r = %d (%s):\n    %s\n" % (r-1, lvl_r.phi, tt)
            leaves.append(tt)
            sides = []





        # Create copies of tt for each side to use
        if len(sides) > 0:
            side_indexes_types = [(0, tt)]
            if len(sides) > 1:
                side_indexes_types += [(i, tt.copy()) for i in [1..len(sides)-1]]
                more_factors = True
            else:
                more_factors = False
            side_indexes_types.reverse()
            print side_indexes_types
        else:
            side_indexes_types = [ ]

        previous_h = 0

        for i, tt in side_indexes_types:
            #for i in reversed([0..len(sides)-1]):
            
            lvl_r = tt.rth_level()
            side = sides[i]
         #   print "Analysing side: %s" % (str(side),)
            print"side",side, side.slope.numerator()
            lvl_r.h = - side.slope.numerator()
            print"Das ist h",lvl_r.h
            lvl_r.e = side.slope.denominator()
            lvl_r.slope = - side.slope
            lvl_r.inv_h = inverse_mod(lvl_r.h, lvl_r.e)
            print"wo is pol"            
            lprime = (lvl_r.inv_h*lvl_r.h - 1) // lvl_r.e 	# should be (1- lvl_r.inv_h*lvl_r.h)
            new_pi = list( lvl_r.inv_h * lvl_r.log_phi - lprime*lvl_r.log_pi )
            new_pi.append(0)
            lvl_r.log_gamma = lvl_r.e*lvl_r.log_phi - lvl_r.h*lvl_r.log_pi

            #E = ZZ(side[2][0] - side[1][0]) # [2][0] - [1][0]
            #H = ZZ(side[1][1] - side[2][1]) # [1][1] - [2][1]
            E = ZZ(side.width())
            H = ZZ(side.height())
            index_N += (E * previous_h) + ((E*H-E-H+(E // lvl_r.e)) // 2)
       #     print "E: %d, H: %d, side: %s, index N: %d" % (E, H, str(side), index_N,)
            previous_h += H
            ## Residual polynomial
            res_pol = tt.residual_polynomial(r, side_devs[i])
            res_pol = res_pol.monic()
            factors = res_pol.factor()
#            print "Res.Pol. factors: %s" % (list(factors),)
            factors_types = [ (factors[0], tt) ]
            print"In Montes_loop fine 1"
            if len(factors) > 1:
                factors_types += [ (factors[i], tt.copy()) for i in [1..len(factors)-1] ]
                more_factors = True
            print"In Montes_loop fine 2"
            for factor, tt in factors_types:
          #      print "Analysing factor of the Res.Pol. %s" % (factor[0],)

           #     print 'old_tt:', tt
                lvl_r = tt.rth_level()

                omega = factor[1]
                lvl_r.res_pol = factor[0]
                lvl_r.f = lvl_r.res_pol.degree()

                ## Representative
                # The new level (lvl_rp1) is already part of the type.
                print"In Montes_loop before repre"                
                lvl_rp1 = tt.representative(omega)
                print"In Montes_loop after repre"                 
  #              print 'new_tt:', tt

                if lvl_r.phi.degree() == lvl_rp1.phi.degree():
                    # Non-optimal, refining level r
#                    print 'Refining, cutting slope: %s' % (str(lvl_rp1.slope),)
                    tt.refinement(more_factors=more_factors)
#                    print"hierrrrr",sides[0].p1,sides[0].slope                   
                    lvl_r.cutting_slope=-ZZ(sides[0].slope)
 #                   print "in refinments",lvl_r.cutting_slope
                else:
                    # Proceeding to higher order
 #                   print 'Proceeding to higher order'
                    
                    tt.log_pi = vector(new_pi)
                    tt.log_phi = -(lvl_rp1.V // lvl_r.e) * lvl_r.log_pi
                    tt.log_phi = vector(list(tt.log_phi) + [1])

                # Push the new or refined type onto the stack.
                type_stack.append(tt)
  #              print "after appending to stack: %s" % (str(type_stack),)

            ## End of `factors' for loop
        ## End of `sides' for loop
        
        total_index += starting_prod_f * index_N
#        print "Added %d * %d to the index (--> %d)" % (starting_prod_f, index_N, total_index,)

    return leaves, total_index

def phi_expansion(f, phi, omega):
    q = f
    coeffs = [ ]
    quos = [ ]
    Ax=f.parent()
    for j in [0..omega]:
        q, r = q.quo_rem(phi)
        coeffs.append(Ax(r))
        quos.append(Ax(q))

    return coeffs, quos

def path_of_precision(n, h):
    q = n
    path = [ n ]

    while q > h:
        q, a = q.quo_rem(2)
        q += a
        path.insert(0, q)

    return path

def change_precision(a, prec):
    if not isinstance(a.parent(), sage.rings.padics.local_generic.LocalGeneric):
        # We assume it's a polynomial.
        return a.parent()([ change_precision(c, prec) for c in a])

    if prec > a.precision_absolute():
        a = a.lift_to_precision(prec)
    else:
        a = a.add_bigoh(prec)

    return a



def index(F):

	f=minimal_polynomial(F.0)
	print"before dsc"
	dsc=f.discriminant().numerator()
	print"after dsc"
	dsc_prime=derivative(dsc)
	facs=list(factor(gcd(dsc,dsc_prime)))
	print"sqr fac"	
	local_index=[]
	for fac in facs:
		om,ind=montes(F,fac[0])
		local_index.append(ind*fac[0].degree())
	print"out"				
	return sum(local_index)	
	
	
def terms(f):

	x=f.parent().0
	L=list(f)
	if len(L) == 1:
		return L
		
	tmp=[]
	for i in [0..len(L)-1]:
		tmp.append(L[i]*x^i)
	
	return tmp
	
	

def infinity_representation(F):

	K=F.base_field()
	T=K.0
	A=T.numerator().parent()
	f=minimal_polynomial(F.0)	
	n=f.degree()
	coeff=[c.numerator() for c in f.coefficients()]
	cf=max([ceil(coeff[j].degree()/(-j+n)) for j in [0..len(coeff)-2]])
	coeff_list=list(T^(-n*cf)*f.parent()(	f(T^cf*f.parent().0))	)
	coeff_newf=[];
	for i in coeff_list:
		tmp=0
		
		for j in terms(i.numerator()):
			j=A(j)
			tmp=tmp+j.leading_coefficient()*T^(i.denominator().degree()-j.degree())
		
		coeff_newf.append(tmp)
		
	f_new=f.parent()(coeff_newf)
	return K.extension(f_new),cf;
	
	
def GENUS(F):

	n=F.degree()
	ind=index(F)
	FF,cf=infinity_representation(F)	
	om_inf,inf_ind=montes(FF,F.base_field().0)
	ind_new=-n-inf_ind+Integers()((cf*(n-1)*n/2))-ind;
	return ind_new+1


def coeff_mod(g,p):

	return g.parent()([i.numerator() % p.numerator() for i in list(g)])

def coeff_div(g,p):

	return g.parent()([i.numerator() // p.numerator() for i in list(g)])







class Error(object):

    def __init__(self):
        self.error = []
        
def into_error(self,elt):
	self.error.append(elt)        



class factional_ideal(object):

    def __init__(self, p,F):
        self.parent = F
        self.polynomial_generator = p




class PrimeIdeal(object):

    def __init__(self, tt):
        self.type = tt
        self.parent = tt.parent
        self.ramification_index = prod([i.e for i in tt.levels])
        self.residual_degree = tt.levels[len(tt.levels)-1].res_pol.degree()
        self.polynomial_generator = tt.prime_polynomial
        self.position = tt.position
        self.okutso_basis = []
        self.uniformizer = 1
        self.exponent = 1        
 
    def __unicode__(self):
        p = self.polynomial_generator
        position = self.position
        exponent = self.exponent
        if exponent == 0:
        	return ''
        elif exponent == 1:
        	return 'P(%s, %s)' % (str(p), str(position))        
        else:	
        	return 'P(%s, %s)^%s' % (str(p), str(position), str(exponent))

    def __str__(self):
        return self.__unicode__()

    def __repr__(self):
        return self.__unicode__()
        
    def __pow__(self,exponent):
        new_prime_ideal = copy(self)
        if exponent == 0:
    		new_prime_ideal.ramification_index = 0
    		new_prime_ideal.residual_degree = 0
    		new_prime_ideal.residual_degree = 0
    		new_prime_ideal.exponent = 0
    		new_prime_ideal.position = 0   
    		new_prime_ideal.polynomial_generator= 1 		    		    		    		    
        else:
	        new_prime_ideal.exponent *= exponent    	    
        return new_prime_ideal
 
    def copy(self):
        copy = PrimeIdeal(self.type)
 #      return copy
 
    def __mul__(P1,P2):
        return Ideal([P1,P2])
        
    def is_prime_ideal(self):
        return type(self) == PrimeIdeal and self.exponent ==1          

    def is_prime_ideal_power(self):
        return type(self) == PrimeIdeal
        
def primes_over_p(F,p):
    om_reps,_=montes(F,p)
    return [PrimeIdeal(i) for i in om_reps]	              
#     def __mul__(P1,P2):
#         
#         return self        

def sum2(iterable, start=0):
    return start + reduce(operator.add, iterable)

         
class Ideal(object):

    def __init__(self, prime_ideal_list):
    
        if len(Set([i.parent for i in prime_ideal_list])) > 1:
        	raise Exception, "Error: All prime ideals must belong to the same field."
        primes = [] 
        indices = []        
        for P in prime_ideal_list:
        	tmp = (P.polynomial_generator,P.position)
        	if tmp in indices:
        		l = indices.index(tmp)
        		primes[l].exponent += P.exponent 

        	else:
        		primes.append(P)	
        		indices.append(tmp)        		
        prime_ideal_list = primes		
        prime_ideal_list.sort()	
        self.prime_ideal_list = prime_ideal_list
        self.parent = prime_ideal_list[0].parent
        self.degree = sum([j.ramification_index*j.residual_degree*j.exponent 
        												for j in prime_ideal_list])        
        self.two_element_representation = []
        self.basis = []
        self.polynomial_support = list(Set([l.polynomial_generator 
														for l in prime_ideal_list]))
        self.index = 1 #muss noch angepasst werden

         															
    def __unicode__(self):

        list_string = []
        P_list = self.prime_ideal_list

        for j in range(len(P_list)-1):
			
            list_string.append(str(P_list[j]) + '*')
        list_string.append(str(P_list[len(P_list)-1]))        		
        return sum2(list_string, '')

    def __str__(self):
        return self.__unicode__()

    def __repr__(self):
        return self.__unicode__()
        
    def copy(self):
        copy = Ideal(self.prime_ideal_list)        
        
    def __pow__(self,exponent):
        new_Ideal = copy(self)
        return Ideal([i^exponent for i in new_Ideal.prime_ideal_list])
 
    def __eq__(I1,I2):
        return I1.prime_ideal_list == I2.prime_ideal_list

    def __mul__(I1,I2):
        return Ideal(I1.prime_ideal_list + I2.prime_ideal_list)
        
    def is_ideal(self):
        return type(self) == Ideal
 
    def is_prime_ideal(self):
        P = self.prime_ideal_list[0]		    
        if len(self.prime_ideal_list) == 1 and P.is_prime_ideal():
        	self = P
        	return True
        else:	
        	return False     
        	
    def is_prime_ideal_power(self):
        P = self.prime_ideal_list[0]		    
        if len(self.prime_ideal_list) == 1 and P.is_prime_ideal_power():
        	self = P
        	return True
        else:	
        	return False    
        	
        	
        	
        	
        	
################################### Valuation ##############################

def localize(alpha,p):
	R = alpha.element().parent()
	p = p.numerator()
	if alpha == 0:
		return 1,0,R(0)
	else:
		num = [i.numerator() for i in alpha.list()]
		val_num = min([j.valuation(p) for j in num])
		denom =  alpha.element().denominator()  				       	
		val_denom = denom.valuation(p)
		den = (denom/p^val_denom).numerator()
		
	return den, val_num-val_denom, R(num) // p^val_num	# Vielleicht statt //, coeff_div	        	         	  
	
	
	
	
def VValuation(pol,poly):
	ord = -1
	rem = 0
	p1 = pol.numerator()
	poly = poly.numerator()
	while rem == 0:
		p1, rem = p1.quo_rem(poly)
		ord+=1
	
	return ord	
	
	
	
def P_valuation(alpha,P):	
	if not P.is_prime_ideal:
        	raise Exception, "Error: P must be a prime ideal."

	F = P.parent
	if not alpha.parent() == F:	
        	raise Exception, "Error: Arguments should lie on the same function field."	

	if alpha == 0:
		return +Infinity, alpha
		
	if not type(P) == PrimeIdeal:
		P = P.prime_ideal_list[0]
		
	p = P.polynomial_generator
	tt = P.type
	lvl_0 = tt.levels[0]
	Fq = lvl_0.Fq
	Fqy = lvl_0.Fqy	
	e_P = P.ramification_index
	map = Fq.reduction_map()
	den,exp,numPol = localize(alpha,p)
	cua = exp*e_P
	pphi = [map(i.numerator()) for i in list(lvl_0.phi)]
	pnumPol = [map(i.numerator()) for i in list(numPol)]
	if VValuation(Fqy(pnumPol),Fqy(pphi)) == 0:
		return cua #,reduction
	res_pol = 0
#	z = 0
#	dev = [] 
#	val = 0
	value = 0
	i = 0	
#	r = len(tt.levels)	
	print"asdasdas 1"
 	while value == 0:
 		if i < len(tt.levels):
 			i+=1
 		else:
 			tt.single_factor_lifting_Up(2*tt.rth_level().h)
 		print"asdasdas 2", numPol,i,tt	
		val, dev = tt.value(i+1,numPol) # can be i+1
 		print"asdasdas 3"	 		
  		res_pol = tt.residual_polynomial(i,dev) # can be i
		if VValuation(res_pol,tt.lvl(i).res_pol) == 0:
			value = val*(e_P // (tt.lvl(i).e*tt.lvl(i).prod_e))   
	return value+cua
		
		
		
		
		
		
		
		 