load("~/Dropbox/programming/Sage/Function_Fields/finite_fields_extra.sage")

# -*- coding: utf-8 -*-

import pprint
from itertools import product
pp = pprint.PrettyPrinter(indent=4)




def is_prime_polynomial(g):
	"""
	Returns True if g is a prime polynomial, and False otherwise
	"""
	fac = g.factor()
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

    def __init__(self, F, p, varphi, omega, gf_ext_base, initial=True):
        self.parent = F
        self.pol = F.polynomial()	
        self.prime = p
        self.varphi = varphi        
        self.levels = [ ]
        self.sfl = [0, 0, 0, 0] # Single Factor Lifting
        self.phiadic = [F.base_field()(0) for i in range(0, 4)]
        self.local_index = 0	 
		
        self.gf_extension = gf_ext_base
        #FiniteField(len(p.parent().base_ring())^p.degree(),conway=True,prefix='zini') # has to be adjusted
        if initial is True:
            self.add_initial_level(varphi, omega, p)
    
    def copy(self):
        copy = MontesType(self.parent, self.prime, self.varphi,
                          self.levels[0].omega, self.gf_extension, initial=False)
        for level in self.levels:
            copy.levels.append(level.copy())
        return copy

    def add_initial_level(self, varphi, omega, p):
 #       print"in initial level" ,varphi       
        phi = self.pol.parent()(varphi)		#sollte richtig sein
        new_level = MontesTypeLevel(phi, omega, p)
#        print"ini level 1"
        F_0 = self.gf_extension.ext_field
#        print"in initial level 2"     
#        into_error(EE,[F_0,varphi])
	#the next line seems wrong if type has okutsu depth > 1 --> I am fixing it!!!
        F_1.<z0> = F_0.extension(varphi.degree())
#        print"in initial level 3"		

        new_level.gf_extension = GF_Extension(F_1, F_0, varphi)
        new_level.prod_f = varphi.degree()		
 #       print"in initial level 3.1"		
        new_level.Fq = F_1 #Fp.extension(varphi, name='z0', prefix='ww')
 #   	print"initial_level: 4"	        
        if varphi.degree() > 1:
        # Fehlerquelle 1: Die if abfrage is im orginal anders!
        # Fehlerquelle 2: Ich nehme in der nächsten Zeile schon entsprechende Nst von pol, muss aufpassen wenn ich später lifte und übersetzte!
        	new_level.z = new_level.gf_extension.root_pol  					#new_level.Fq.0                                    
        else:
        	new_level.z = -varphi.list()[0]#-varphi.coeffs()[0]
        new_level.Fqy = PolynomialRing(new_level.Fq, 'y1')
        new_level.inv_M = matrix([])
        self.levels.append(new_level)


        
    def add_new_level(self, phi, omega, u):
        s = len(self.levels)
        lvl_s = self.levels[-1]
        new_level = MontesTypeLevel(phi, omega, self.prime)
        new_level.V = lvl_s.e*u
        new_level.prod_e = lvl_s.prod_e * lvl_s.e
        new_level.prod_f = lvl_s.prod_f * lvl_s.f     


 #       q=(len(self.prime.parent().base_ring())^self.prime.degree())^new_level.prod_f
 #       print"in add_new_level",s,self.lvl(1).res_pol
#        print"\n\n\n\n\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"

   #     if s > 1:
#        print"this is s",s

        K = self.lvl(s).gf_extension.ext_field
        F = K.extension(lvl_s.res_pol.degree(), 'z'+str(s))#GF(len(K)^lvl_s.res_pol.degree(), 'z'+str(s))
        new_level.gf_extension = GF_Extension(F,K,lvl_s.res_pol)
         		
        new_level.Fq = new_level.gf_extension.ext_field
 #       print"mal sehen "  ,new_level.Fq 
        new_level.Fqy = PolynomialRing(new_level.Fq, 'y'+str(s))
        #new_level.embedding = Hom(lvl_s.Fq, new_level.Fq).list()[0]
        new_level.z = new_level.gf_extension.root_pol
        self.levels.append(new_level)
        return new_level

    def refinement(self, more_factors=False):



        r = len(self.levels) - 1
        if more_factors is True:

            self.lvl(r).refinements.append([self.lvl(r).phi,
                                            self.lvl(r).slope])
        
        #self.lvl(r).cutting_slope = ZZ(self.lvl(r+1).slope)  # Error?
        self.lvl(r).phi = self.lvl(r+1).phi
        self.lvl(r).omega = self.lvl(r+1).omega
#        print"factor",self.lvl(r+1).omega
#        self.lvl(r).gf_extension = self.lvl(r+1).gf_extension

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
#            print"sloppy in add_last_level",side,side.slope        
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
        #assert i <= len(self.levels)

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

        val = +Infinity
        if a == 0:
            if i == 1:
                devs = 0
            else:
                devs = []
            return val, devs
        if i == 1:

            val = min([valuation(c, self.prime) for c in a.numerator().list()])
    #        print"\n\n val",val
        
        
  #          val = min([valuation(c, self.prime) for c in list(a)])
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
#            print"daten in Value",step,lvl_im1.e*2
            while quot != 0 and min_height < val:
                quot, ak = quot.quo_rem(lvl_im1.phi)
#                print"ak",ak
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

        assert i <= len(self.levels)
        lvl_i = self.lvl(i)
       	p = self.lvl(1).prime.numerator()
        #height = side_devs[-1][1]
        gf_ext = self.lvl(1).gf_extension 
        k = gf_ext.ext_field
        res_coeffs = [k(0) for jj in [1..len(side_devs)] ]  #-1


        if i == 1:

        	height = side_devs[-1][1]
        	gf_base = self.gf_extension
        	for ind in range(len(side_devs)-1):	        
        		dev = side_devs[ind]
		        if not dev == 0:

		        	coeff = dev // self.prime^(height)


		        	Fqy = self.lvl(1).Fqy
#		        	print"new idea"
		        	tmp = [(i.numerator() % p).list() for i in list(coeff)]
		        	zz = gf_base.root_pol
#		        	print"new idea2"
		#        	print"new uuu",zz		        	
		        	tmp_2 = [ sum([Embed(gf_base,ll[jj])*zz^jj for jj in range(len(ll)) ]) for ll in tmp]
		        	coeff_embed = [Embed(gf_ext,co) for co in tmp_2 ]
		        	zz_2 = gf_ext.root_pol
		    #    	print"wwaewawe",Embed(gf_ext,tmp_2[0]).parent(),zz_2.parent()
		     #   	print"new idea3",		  sum([Embed(gf_ext,tmp_2[pp])*zz_2^pp for pp in range(len(tmp_2))]).parent()
	#	        	coeffs = []
		        	#coeffs = [sum([oo[l]*gf_ext.root_pol^l for l in range(len(oo))]) for oo in coeff_embed]    
	#	        	print"inne: 2",[lo.minimal_polynomial() for lo in coeffs]
	#	        	print"inne: 3",coeffs[0]
		#        	coeff = Fqy(coeffs)
#		        	print"inne: 1", [zu.minimal_polynomial() for zu in coeff_embed]
		
	#	        	print"inne: 4",coeff.parent()
		        	es = sum([Embed(gf_ext,tmp_2[pp])*zz_2^pp for pp in range(len(tmp_2))])
		        	res_coeffs[ind] = es#coeff(gf_ext.root_pol)
#		        	print"inne: 2", res_coeffs[ind].minimal_polynomial()
		        
		        height -= lvl_i.h	
	#	        print"asdas"
        else:
        	res_coeffs = [Embed(lvl_i.gf_extension,tr) for tr in res_coeffs]
        	lvl_im1 = self.lvl(i-1)
        	l_prime = (1 - lvl_im1.inv_h*lvl_im1.h) // lvl_im1.e      
        	
        	for ind2 in range(len(side_devs)-1):	        
        		dev = side_devs[ind2]
        		if not len(dev)==0:

		        	lvl_im1 = self.lvl(i-1)
		        	twist_exp = l_prime*dev[-1][0] - lvl_im1.inv_h*dev[-1][1] 
#		        	print"gehe in res_pol rein mit", dev		
		        	pj = self.residual_polynomial(i-1, dev)
#		        	print"gehe in res_pol raus",pj
#		        	print"unnnnnnnd", lvl_i.z.parent()
		        	tmp = pj.list()
		        	eval = [ Embed(lvl_i.gf_extension,tmp[m])*lvl_i.z^m for m in range(len(tmp))]
#		        	coeff_tmp = [coeff[ll]*lvl_i.z^ll for ll in range(len(pj))]
#		        	
#		        	coeff_embed = [Embed(self.lvl(i).gf_extension,c) for c in coeff]
#		        	lifted_coeff = lvl_i.Fqy(coeff_embed)## muss man dringend anpassen
		        	res_coeffs[ind2] = lvl_i.z^(twist_exp) * sum(eval)
#		        	print"embed works",res_coeffs[ind2].minimal_polynomial()

#        print"baue pol"#,	lvl_i.Fqy,[o.parent() for o in res_coeffs]
        res_pol = lvl_i.Fqy(res_coeffs)
#        print"out res_pol",res_pol

#        print"AAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\n\n\n"
        return res_pol

    ## Representative ##
    def representative(self, omega):
        """
        Construct a representative phi of a type. A new level is added with
        phi and V.
        """



 
 #       if len(self.levels) > 1:
  #      	return
        
        s = len(self.levels)
        lvl_s = self.lvl(s)
        ef = lvl_s.e * lvl_s.f
        u = ef * lvl_s.V+lvl_s.f*lvl_s.h

#        if s > 1:
#            txp = -self.lvl(s-1).inv_h * (u // self.lvl(s-1).e)
#            twist = lvl_s.z^(txp)
#        else:
#            twist = lvl_s.Fq(1)

        # Reductum from magma
#        res_pol = twist * lvl_s.Fqy(list(lvl_s.res_pol)[:-1])
        res_pol = lvl_s.Fqy(list(lvl_s.res_pol)[:-1]) 
#        print"\n\n skdjasdjalksjdlkasjdlasjdljasldjaslkdjas",lvl_s.z.minimal_polynomial()
 #       print"res_pol.mipo",[uu.minimal_polynomial() for uu in lvl_s.res_pol.list()]
  #      print "%s = %s * Reductum(%s) = %s * %s" % (res_pol, twist, lvl_s.res_pol, twist, lvl_s.Fqy(list(lvl_s.res_pol)[:-1]))
     #   u += lvl_s.f * lvl_s.h
     
        phi0 = self.construct(s, res_pol, 0, u)
   #     print "%d. Repr. = %s + phi^ef" % (s, phi0,)
        phi0 = phi0 + lvl_s.phi^(ef)

        new_level = self.add_new_level(phi0, omega, u)

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
#        print"in construct with",res_pol.parent()
        var = lvl_i.phi^(lvl_i.e)

        phi0 = 0
        height = u - res_pol.degree()*lvl_i.h

        if i == 1:
            for a in reversed(list(res_pol)):
#                print"Dataa", a, a.minimal_polynomial()

                a_tmp = elt_seq(self.lvl(1).gf_extension,a)
                A = self.prime.numerator().parent()         
                a_tmp2 = [A(elt_seq(self.gf_extension,re)) for re in a_tmp]
                Ax = self.pol.parent()
                lift = Ax(a_tmp2)   
     

                phi0 = (phi0 * var) + (lift * self.prime^height)

                height += lvl_i.h
        else:
 #           print"in construct data (i=2) ____________________________________"
            step = (lvl_i.e * lvl_i.V) + lvl_i.h
  #          print"\n step",step, lvl_i.e, lvl_i.V, lvl_i.h
            new_V = u - (res_pol.degree() * step) - (s * lvl_i.V)
   #         print"\n new_V",new_V;
            lvl_im1 = self.lvl(i-1)
            for a in reversed(list(res_pol)):
                pj = 0
                if a != 0:
                    txp, s_im1 = divmod(lvl_im1.inv_h*new_V, lvl_im1.e)
    #                print"txp, s_im1 ",txp, s_im1 
                    u_im1 = (new_V - (s_im1 * lvl_im1.h)) // lvl_im1.e
     #               print"u_im1 ",u_im1 
                    c = (a*lvl_i.z^txp)

                    gf_ext = lvl_i.gf_extension
      #              print"c ",[lo.minimal_polynomial() for lo in elt_seq(gf_ext,c)] 
                    
                    Fqy = lvl_im1.Fqy
                    new_res_pol = Fqy( elt_seq(gf_ext,c))
                    pj = self.construct(i-1, new_res_pol, s_im1, u_im1)
                phi0 = (phi0 * var) + pj
                new_V += step
            

        phi0 = phi0 * lvl_i.phi^(s)
       # print"in construct data (i=2) <><><><><><><><><><><><><><><>_ ende",   phi0
#        into_error(EE,phi0)
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
        p = self.prime
        exponent = self.sfl[0]
        nu = self.sfl[1]
        x0prec = self.sfl[2]
        x0num = self.phiadic[3]
        x0den = self.sfl[3]
        print"hier schon falsch?", x0num

        prod_e = lvl_r.prod_e
        h = lvl_r.h - lvl_r.cutting_slope
#        print"cutting slope "  ,lvl_r.cutting_slope
 #       print"h= "  ,lvl_r.h
        last_h = slope - lvl_r.cutting_slope
        V = lvl_r.V + lvl_r.cutting_slope

        precision=nu + exponent + ceil((V+last_h)/prod_e)
        p_prec = p^precision

        pol_zp = coeff_mod(self.pol, p_prec)
        psinum_zp = coeff_mod(self.phiadic[2] , p_prec) #PsinumZp:=type[1]`Phiadic[3] mod p_prec;
#        print"x0prec", x0prec         
        path = path_of_precision(last_h, h)
        short_path = path_of_precision(h, x0prec)
  #      print"also hier 1"   

        newprecision= nu + exponent + ceil(h/prod_e) #Ceiling((V+path[2])/e)
        p_prec_new=p^newprecision;
        a1=coeff_mod(self.phiadic[1] , p_prec_new)		#type[1]`Phiadic[2] mod p_prec_new;
		#print "That precision:", nu + exponent + ceil(h/prod_e);

        newprecision=nu + exponent + ceil((V+path[1])/prod_e)
        p_prec_new=p^newprecision
        phi=coeff_mod(lvl_r.phi , p_prec_new)
        psinum= coeff_mod(psinum_zp , p_prec_new )
#        print"also hier 3", ((coeff_mod(self.phiadic[0] , p_prec_new ))*psinum) , phi
        



        a0num, a0den = self.cancel((coeff_mod(self.phiadic[0] , p_prec_new )*psinum) % phi, nu)
#        print "a0 num, den:", [a0num, a0den]
        a1num, a1den = self.cancel((coeff_mod(a1 , p_prec_new)*psinum) % phi, nu)
#        print "a1 num, den:", [a1num, a1den], [a0num, a0den], short_path

        print "--==--==--==--==--==--==--"

        for i in range(1, len(short_path)):
            low_precision = a1den + 2 * x0den  + ceil(short_path[i]/prod_e)
            x0num, x0den = self.inversion_loop([ a1num, a1den], x0num, x0den,
                                               phi, low_precision,p)
 #           print "x0 num, den:", x0num, x0den
            print '----'

        print "--==--==--==--==--==--==--"
        print"ax0num, x0den",x0num, x0den
        print"################################"	
        anum, aden = self.cancel((a0num*coeff_mod(x0num ,  p_prec_new)) % phi, x0den+a0den)


        phi = phi + anum

#        print "between:", [anum, aden, phi]
        print "--==--==--==--==--==--==--"

        for i in range(1, len(path)-1):
            loop_prec = nu + exponent + ceil((V+path[i+1])/prod_e)
            ploop=p^loop_prec
            phi=coeff_mod(phi , ploop)
            Psinum = coeff_mod(psinum_zp , ploop)

            qq, c0 = pol_zp.quo_rem(phi)
            print"data",qq,phi;    
            c1 = coeff_mod(qq , phi)

            c0num, c0den = self.cancel(coeff_mod(c0*psinum, ploop) % phi, nu)
            c1num, c1den = self.cancel(coeff_mod(c1*psinum, ploop) % phi, nu)
#            print "cs:", [c1, c0num, c0den, c1num, c1den]





            low_precision = c1den + 2 * x0den + ceil(path[i]/prod_e)
            x0num, x0den = self.inversion_loop([c1num, c1den], x0num, x0den,
                                               phi, low_precision,p)
            print "x0 num, den", x0num, x0den

            xnum =x0num % ploop
            cnum, cden = self.cancel((c0num * coeff_mod(x0num , ploop)) % phi, x0den)
            phi = coeff_mod((phi + cnum) , ploop)

            print "phi:", phi
            print "--------------------------"


#        print "--==--==--==--==--==--==--"

        self.sfl[2] = max(h, path[-2])
        lvl_r.phi = coeff_mod(phi,p_prec) #lvl_r.phi.parent()(phi)
        self.phiadic[3] = x0num
        self.sfl[3] = x0den

        
    def sfl_init(self):
    	print"in SFL_init \n\n\n\n\n\n\n\n\n"
    	print"____________________________________________________________________________________"
        p = self.prime.numerator()
        lvl_r = self.rth_level()
        prod_e = lvl_r.prod_e
        a1 = self.phiadic[1] # FIXME: Why [1]?
        gf_ext = lvl_r.gf_extension
        A=p.parent()

        Ax.<x>=PolynomialRing(A)
        psinum = Ax(1)


        
        r = len(self.levels) - 1

        if r == 0:
            nu = min([valuation(a.numerator(), p) for a in a1])
            a1	=a1 // p^nu 
        
            F_0 = self.lvl(1).gf_extension.base_field

            F_0_y0 = self.lvl(1).Fqy


            coeff_embed = [[Embed(gf_ext,co) for co in list(i.numerator() % p)] for i in list(a1)]
            coeffs = [sum([co[l]*gf_ext.root_pol^l for l in range(len(co))]) for co in coeff_embed]
            klass = F_0_y0(coeffs)

            # Evaluate a1/p^nu in z_1 (this may be z_0)
#            print"klass",klass          
        else:
            val, dev = self.value(r+1, a1)
            print"val,dev",val,dev
            print"\n\n\n ___________0"
            res_pol = self.residual_polynomial(r, dev)
            print'res_pol',res_pol
            print"\n\n\n ___________1"            

            logpsi = 0
            qq, s = (-val).quo_rem(prod_e)
            psinum, logpsi = self.prescribed_value(s)
            print'psinum',psinum,logpsi
            print"\n\n\n ___________2"            
            nu = -logpsi[0] - qq
            
            vector = dev[-1][0] * self.lvl(r).log_phi + dev[-1][1]*self.lvl(r).log_pi
            print'vector',vector
            print"\n\n\n ___________3"            


            klass = self.convert_logs(logpsi + vector)
            print'klass',klass
            print"\n\n\n ___________4"                       
            tmp = res_pol.list()
            z_tmp = self.lvl(r+1).z
            embed_res_pol = sum([Embed(self.lvl(r+1).gf_extension,tmp[ex])*z_tmp^ex for ex in range(len(tmp))])
            klass = long_embed(self,z_tmp.parent(),klass)
            klass *= embed_res_pol
            print'klass',klass
            print"\n\n\n ___________5"              
            
        self.phiadic[2] = psinum
        self.sfl[1] = nu
        self.sfl[2] = 1
        print'psinum',psinum
        print"\n\n\n ___________6"  
        x0num = 0
        x0den = 0
        x0num, x0den = self.local_lift(klass^(-1))
        print'x0num, x0den',x0num, x0den, klass^(-1)
        print"\n\n\n ___________7"  

        self.phiadic[3] = x0num
        self.sfl[3] = x0den
        print"aus sfl_ini raus"

    def SFL(self, slope):
        
        self.single_factor_lifting(slope)
        self.update_last_level()

    def prescribed_value(self, value):
        """
        From +Ideals:
        If we are attached to the prime ideal P with Okutsu depth r, then
        logpsi=[a_0, ..., a_r] and psi=phi_1^a_1 ... phi_r^a_r, with
        v_P(p^a_0 psi(theta))=value.
        """
 #       print"\n\n parent hierrrr", self.lvl(1).phi.parent(),
 #       print"\n\n"       
        Ax = self.lvl(1).phi.parent()
 #       print"Ax",Ax
        psi = Ax(1)
        r = len(self.levels)
        logpsi = (ZZ^r)(0)
  #      print"check 1"
        qq, val = value.quo_rem(self.rth_level().prod_e)
   #     print"check 2",range(r-1)
        logpsi[0] = qq
        if val > 0:
            body = val
            for k in reversed(range(r-1)):
                jj = (self.levels[k].inv_h * body) % self.levels[k].e
                logpsi[k+1] = jj
                psi = psi * self.levels[k].phi^jj
                res = (body - jj*self.levels[k].h) // self.levels[k].e
                body = res - jj*self.levels[k].V
            logpsi[0] += res
    #    print"out prescribed_value"
        return psi, logpsi

    def local_lift(self, clss):
        """
        From +Ideals:
        class should belong to the residue class field  type[r]`Fq. The output
        is a pair g,e such that g(theta)/p^e is a lift to a P-integral element
        in K and deg g(x)<n_P.
        """
        i = 1
        print"in locallift _______________",clss.parent()
        
        print"\n\n\n\n"        
        while len(clss.parent()) != len(self.lvl(i).Fq):
            i += 1
        print"out i",i
        if not clss in self.lvl(i).Fq:
			clss = long_shift_down(self,self.lvl(i).Fq,clss)
#			into_error(EE,[clss,self])
			
        if i == 1:
            #Ax.<x>=PolynomialRing(self.prime.numerator().parent())
            A = self.prime.numerator().parent()
#            clss_coeff = elt_seq(self.lvl(1),clss)

            print"what", 		elt_seq(self.lvl(i).gf_extension,clss)[0].parent()	           
            clss_coeff = elt_seq(self.lvl(i).gf_extension,clss)
            tmp = [elt_seq(self.gf_extension,cl) for cl in clss_coeff ]
            numlift = A([i_1  for j_1 in range(len(tmp)) for i_1 in tmp[j_1]])
            denlift = 0
        else:
            #raise(NotImplemented, "local lift not implemented for i > 1.")
            expden = ceil(self.lvl(i).V/self.lvl(i).prod_e)			
            V = self.lvl(i).prod_e*expden
            log = V*self.lvl(i).log_pi
            log = vector(log.list()[:-1])
            newclss = self.convert_logs(log)    		
            H = V // self.lvl(i-1).e
            print"next step might go wrong 1"
            elt = self.lvl(i).z^(self.lvl(i-1).inv_h*H)*clss*newclss^(-1)
            print"next step might go wrong 2"
            c = self.lvl(i).Fq(elt)
            print"next step might go wrong 3"
            varphi = self.lvl(i-1).Fqy(c.polynomial().list())			
            print"next step might go wrong 4"            
            lift = self.construct(i-1, varphi, 0, H)
            print"next step might go wrong 5"            
            v1lift = min([rr.valuation(self.lvl(1).prime) for rr in lift.coefficients()])
            numlift = lift // self.lvl(1).prime^v1lift
            denlift = expden-v1lift
        print"\n\n\n\n"        
        print"out locallift _______________"


        return numlift, denlift

    def convert_logs(self, log):
        """"
        From +Ideals:
        log[1] is not used. The product of all Phi_i^log[i] for i>0 should have
        integer value M.

        The output is the class of this product divided by p^M.
        """

        vec = log
        z = 0
        klass = self.lvl(1).Fq.prime_subfield()(1)
        for i in reversed(range(len(vec))[1:]):

            ti = vec[i] // self.lvl(i).e

            z = Z_val(self,i)			#self.levels[i-1].z
            klass *= z^ti

            vec = vec - ti*self.levels[i-1].log_gamma

            vec = vector( vec.list()[:vec.degree()-1] );


        return klass

    def cancel(self, poly, den):
        if poly == 0:
            return poly, 0
#        print 'In cancel'
        cancel = min([den] + [ a.valuation(self.prime) for a in poly])
        zq = poly[0].parent()
#        print "poly, p^cancel: {0}, p^{1}".format(poly, cancel)
        num = poly.parent()([ poly.parent().base_ring()(c) / self.prime^cancel for c in poly ])
#        print "num: {0}".format(num)
        #num = poly / self.prime^cancel

        return num, den-cancel

    def inversion_loop(self, A, xnum, xden, phi, precision,p):
        anum = A[0]
        aden = A[1]
        print"ini ",xnum, xden
#        print'in inverseloop'
#        zq = Zp(self.prime, precision, type='capped-abs', print_mode='terse')
#        zqt.<t> = PolynomialRing(zq)

        p_prec=p^precision
        phip = coeff_mod(phi,p_prec) 
        xnum = coeff_mod(xnum,p_prec)
        print"beforehand", p_prec,phi
        print"\n\n",phip
        print"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        print"1",	anum

        anum_mod = coeff_mod(anum,p_prec)
        print"2",	anum_mod,xnum
        tmp = coeff_mod(anum_mod*xnum,p_prec)
        print"3",	tmp
        print"input cancel",(2*self.prime^(xden+aden)-tmp) % phip, xden+aden

        x1num, x1den = self.cancel((2*self.prime^(xden+aden)-tmp) % phip, xden+aden)
        print"in inversion_loop: x1num,x1den",x1num,x1den
        print"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
#        print"bin bis hier 2"
        xnum, xden = self.cancel((xnum*x1num) % phip, xden+x1den)
        xnum = coeff_mod(xnum,p_prec)
        return xnum, xden
        
    def update_last_level(self):
        
        pol = self.pol
#        print"\n\n\n in UpdateLasteLevel"
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
#            print"what 5" , tmp
#            print"sdaa",lvl_r.res_pol             
            lvl_r.log_gamma = lvl_r.log_phi-lvl_r.h*lvl_r.log_pi
			 

class MontesTypeLevel(object):

    def __init__(self, phi, omega, p):
        self.phi = phi
        self.prime = p

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
        self.gf_extension = None

        self.slope = None
        self.res_pol = None

    def copy(self):
        copy = MontesTypeLevel(self.phi, self.omega, self.prime)

        copy.phi = self.phi
        copy.prime = self.prime

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
        copy.gf_extension = self.gf_extension
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
    print"Montes for prime p=",p
    if is_prime_polynomial(p.numerator()) is False:
        raise Exception, "Error: p must be prime."


    f = F.polynomial()					#K.defining_polynomial()

    if f.is_monic() is False:
    #sollte monic_integral_model verwenden um normiertes Polynom zu erzeugen
        raise Exception, "Error: def.pol. of F must be monic."
    # TODO: Add requirements that coefficients of f are integral

    reps_OM = [ ]
    trees = [ ]
    total_index=0
    p=p.numerator()
    A=p.parent()

    F_base = p.parent().base_ring() 
    F_0=FiniteField(len(A.base_ring())^p.degree(),prefix='zini') #conway=True,
    #F_0.<zini> = F.constant_base_field().extension(p)
    F_0_y0=PolynomialRing(F_0, 'y0')
    gf_ext_base = GF_Extension(F_0,F_base,p.numerator())

    coeff_embed = [[Embed(gf_ext_base,co) for co in list(i.numerator() % p)] for i in list(f)]
    coeffs = [sum([co[l]*gf_ext_base.root_pol^l for l in range(len(co))]) for co in coeff_embed]
    res_pol=F_0_y0(coeffs)
#    else:
#		res_pol=F_0_y0([F_0(i.numerator() % p) for i in list(f)])
 
    
    pol_factor = res_pol.factor()

    for fac in pol_factor:

#	print "Analysing irred. factor modulo p: %s" % (str(factor[0]),)

        tt = MontesType(F, p, fac[0], fac[1], gf_ext_base)
        tree, tree_index = montes_main_loop(F, p, tt)
        total_index += tree_index	
        reps_OM += tree
        trees.append(tree)
    if len(reps_OM) == 1:
    	f = reps_OM[0].lvl(1).phi.parent()(f.list())
        reps_OM[0].rth_level().phi = f
        reps_OM[0].rth_level().slope = +Infinity

    for P in reps_OM:
		l_gen,log = P.prescribed_value(1)		#element_with_prescribed_value(P,1)
		P.local_generator = l_gen(K.0)*p^log[0]
		P.log_lg = log		
		P.local_index = total_index

    return reps_OM
    
def montes_main_loop(K, p, tt):

    total_index = 0
    leaves = [ ]
    type_stack = [ tt ]

    while len(type_stack) > 0:
     #   print "STARTING LOOP:\n  type stack:", type_stack
        tt = type_stack.pop()
        into_error(EE,tt.copy())

        r = len(tt.levels)
        lvl_r = tt.rth_level()


#        print"phi, omega",lvl_r.phi,lvl_r.omega
 #       print"\n\n\n",lvl_r.res_pol
  #      print"\n\n\n\n\n\n\n__________"
        phiadic, quotients = phi_expansion(tt.pol,
                                           lvl_r.phi,
                                           lvl_r.omega)
        ## Phi-Newton Polygon

        sides, side_devs = tt.phi_newton_polygon(r, phiadic)

#        print"\n\n ++++++++++",side_devs;
 #       print"asdasdasdasdasdasdasdasdasdasdad", phiadic;
#        print "Sides of Newton polygon: %s" % (sides,)         
        length_N = lvl_r.omega
        index_N = -lvl_r.cutting_slope * (length_N*(length_N-1) // 2)

       # starting_prod_f = lvl_r.prod_f
       # print "index N: %d, length N: %d" % (index_N, length_N)
#        print"In main loop: 1"  
  
        if length_N == 1:
#            print"laste leven daten",phiadic,sides[0],side_devs[0]
            tt.add_last_level(phiadic, sides[0], side_devs[0])

#            print "  #### Found a factor of depth r = %d (%s):\n    %s\n" % (r-1, lvl_r.phi, tt)
            leaves.append(tt)
            sides = []


#        print"In main loop: 2"

 #       print"\n\n\n\n\n\n\n\n\n --------------------------------------------------- A"
        # Create copies of tt for each side to use
        if len(sides) > 0:
            side_indexes_types = [(0, tt)]
            if len(sides) > 1:
                side_indexes_types += [(i, tt.copy()) for i in [1..len(sides)-1]]
                more_factors = True
            else:
                more_factors = False
        #    side_indexes_types.reverse()
         #   print side_indexes_types
        else:
            side_indexes_types = [ ]

        previous_h = 0
  #      print"In main loop: 3"
        for i, tt in side_indexes_types:
  #          print"\n\n\n\n\n\n\n\n\n --------------------------------------------------- 0"
            #for i in reversed([0..len(sides)-1]):
            
            lvl_r = tt.rth_level()
            side = sides[i]
            lvl_r.h = - side.slope.numerator()
            lvl_r.e = side.slope.denominator()
            lvl_r.slope = - side.slope
            lvl_r.inv_h = inverse_mod(lvl_r.h, lvl_r.e)
            lprime = (1-lvl_r.inv_h*lvl_r.h ) // lvl_r.e 	# should be (1- lvl_r.inv_h*lvl_r.h)
            new_pi = list( lvl_r.inv_h * lvl_r.log_phi - lprime*lvl_r.log_pi )
            new_pi.append(0)
            lvl_r.log_gamma = lvl_r.e*lvl_r.log_phi - lvl_r.h*lvl_r.log_pi

            E = ZZ(side.width())
            H = ZZ(side.height())
   #         print"index_N",E, previous_h, H, lvl_r.e
            index_N += (E * previous_h) + ((E*H-E-H+(E // lvl_r.e)) // 2)
            previous_h += H
            ## Residual polynomial
#            print"before res_pol", r, side_devs[i]

            res_pol = tt.residual_polynomial(r, side_devs[i])
  #          print"res_pol: ",res_pol            
            res_pol = res_pol.monic()
            factors = res_pol.factor()
 #           print"in main loop: for: 4"
            factors_types = [ (factors[0], tt) ]
    #        print"\n\n\n\n\n\n\n\n\n --------------------------------------------------- 1"
            if len(factors) > 1:

                factors_types += [ (factors[i], tt.copy()) for i in [1..len(factors)-1] ]
                more_factors = True
    

            for factor, tt in factors_types:
 #           	print"2. for loop", len(factors_types)

                lvl_r = tt.rth_level()

                omega = factor[1]
                lvl_r.res_pol = factor[0]
                lvl_r.f = lvl_r.res_pol.degree()

                lvl_rp1 = tt.representative(omega)
  #          	print"after rep"
     #           print"\n\n\n\n\n\n\n\n\n ---------------------------------------------------2"
                if lvl_r.phi.degree() == lvl_rp1.phi.degree():
	
	
	                    # Non-optimal, refining level r
   #                 print 'Refining, cutting slope: %s' % (str(lvl_rp1.slope),)
  #                  print"before refinement",  tt
  #                  print"\n\n\n\n\n\n\n___________________________________________________________"

      #              print"in refinement"
#                    into_error(EE,tt.copy())
                    tt.refinement(more_factors=more_factors)
#                    print"hierrrrr",sides[0].p1,sides[0].slope                   
                    lvl_r.cutting_slope=-ZZ(sides[0].slope)


#                    print"force quit"
#                    return
      #              print"\n\n\n\n\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

        #            print "after refinments",tt
         #           print"\n\n\n\n\n\n\nüüüüüüüüüüüüüüüüüüüüüüüüüüüüüüüü" 
                else:
                    # Proceeding to higher order
 #                   print 'Proceeding to higher order'
                    
                    tt.log_pi = vector(new_pi)
                    tt.log_phi = -(lvl_rp1.V // lvl_r.e) * lvl_r.log_pi
                    tt.log_phi = vector(list(tt.log_phi) + [1])
                # Push the new or refined type onto the stack.
       #         print"\n\n\n\n\n\n\n\n\n ---------------------------------------------------3"

                type_stack.append(tt)
        #        print"\n\n\n\n\n\n\n\n\n ---------------------------------------------------4"
  #              print "after appending to stack: %s" % (str(type_stack),)

            ## End of `factors' for loop
        ## End of `sides' for loop
        total_index += lvl_r.prod_f * index_N



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

	f = F.polynomial()
	K = f.parent().base_ring()
	A = f.parent().base_ring().0.numerator().parent()
	Ax.<x> = PolynomialRing(A)
	f = Ax([co.numerator() for co in f.list()]) 
	if not is_separable(f):
		raise Error, "Defining polynomial of function field is not separable"
	print"check if pol is separable"
	print"before dsc"
	dsc = f.discriminant().numerator()
	
	print"after dsc"
	dsc_prime = derivative(dsc)
	facs = list(factor(gcd(dsc,dsc_prime)))
	print"sqr fac"	
	local_index=[]
	index_divs = []
	for fac in facs:
		print"geht durch 0",fac[0],fac[0].parent()	

		om = montes(F,K(fac[0]))
		print"geht durch 1"	
		ind = om[0].local_index
		print"geht durch 2"
		local_index.append(ind*fac[0].degree())
		if ind > 0:
			index_divs.append(fac[0])
		
	print"out"				
	return sum(local_index),index_divs
	
	
def terms(f):

	x=f.parent().0
	L=list(f)
	if len(L) == 1:
		return L
		
	tmp=[]
	for i in [0..len(L)-1]:
		tmp.append(L[i]*x^i)
	
	return tmp
	
def Z_val(self,i):

#	print"in Z_val",i, len(self.levels)
	if i == len(self.levels):
		#respol = self.lvl(i).res_pol
		z = -list(self.lvl(i).res_pol)[0]
	else:
		z = self.lvl(i+1).z

	return z		

def infinity_representation(F):

	K = F.base_field()
	T = K.0
	A = T.numerator().parent()
	f = F.polynomial()	
	n = f.degree()
	coeff = [c.numerator() for c in f.list()]
	cf = max([ceil(coeff[j].degree()/(-j+n)) for j in [0..len(coeff)-2]])
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

	n= F.degree()
	ind = index(F)
	FF,cf = infinity_representation(F)	
	om_inf = montes(FF,F.base_field().0)
	inf_ind = om_inf[0].local_index
	ind_new=-n-inf_ind+Integers()((cf*(n-1)*n/2))-ind[0]
	return ind_new+1


def coeff_mod(g,p):
	if g.parent() == p.parent():
		return g % p
	else:	
		return g.parent()([i.numerator() % p.numerator() for i in list(g)])

def vvaluation(pol,poly):
	pol = pol.numerator()
	poly = poly.numerator()
	ord =- 1
	rem = 0
	pl = pol
	while rem == 0:
		pl, rem = pl.quo_rem(poly)
		ord += 1
	
	return ord


def den(alpha):
	return lcm([i.denominator() for i in alpha.list()])
	
def num(alpha):
	return alpha*den(alpha)



def localize(alpha,p):
	A = p.parent()
	if alpha == 0:		
		return 1,0,A(0)
	numm = [j.numerator()  for j in (num(alpha)).list()]
	valnum = min([x.valuation(p)  for x in numm])
	denom = A(den(alpha))
	valden = denom.valuation(p)
	denn = A(denom/p^valden)
	Ax.<x> = PolynomialRing(A)
	return	denn,valnum-valden,Ax(numm).quo_rem(p^valnum)[0]



def Z_val(self,i):
	if i == len(self.levels):
		#respol = self.lvl(i).res_pol
		z = -list(self.lvl(i).res_pol)[0]
	else:
		z = self.lvl(i).z

	return z	



def equalize_logs(log1,log2):
	d = log1.degree()-log2.degree()

	if not d == 0:
		tail = [0 for i in [1..abs(d)]]
		if d  > 0:
			log2 = vector(log2.list()+tail)
		else:
			log1 = vector(log1.list()+tail)
	return log1, log2



def P_valuation(alpha,P):
	r = len(P.levels)
	p = P.prime
	Fp = (P.lvl(1)).Fq	
	Fpy = (P.lvl(1)).Fqy	
	reduction = Fp(0)
	if alpha == 0:
		return Infinity, reduction
	denn, exp, numPol = localize(alpha,p)
	cua = exp*P.lvl(r).prod_e
	pphi =  [Fp(i.numerator() % p) for i in list(P.lvl(1).phi)]	
	pnumPol = [Fp(i.numerator() % p) for i in list(numPol)]
	if vvaluation(Fpy(pnumPol),Fpy(pphi)) == 0:
		print"insider", P.convert_logs(-cua*P.rth_level().log_gamma)
		reduction = P.convert_logs(-cua*P.rth_level().log_gamma)
		reduction *= (Fp(denn))^(-1)*Fpy(pnumPol)(P.lvl(1).z)
		return cua,	reduction
	
	res_pol = 0
	z = 0
	side_devs = []
	val = 0
	value = 0
	i = 0
	while value == 0:
		if i < r:
			i +=1	
		else:
			P.SFL(2*P.rth_level().h)	

		val, side_devs = P.value(i+1,numPol)	
		
		res_pol = P.residual_polynomial(i, side_devs)
		if vvaluation(res_pol,P.lvl(i).res_pol) == 0:
			value = val*(P.lvl(r).prod_e).quo_rem((P.lvl(i).e*P.lvl(i).prod_e))[0]
		
# if RED:

	log = side_devs[len(side_devs)-1][0]*P.lvl(i).log_phi+ side_devs[len(side_devs)-1][1]*P.lvl(i).log_pi
	P_log_gamma, log = equalize_logs(vector([P.rth_level().log_gamma[1]]),log)
 	reduction = P.convert_logs(log-(value+cua)*P_log_gamma)
	z = Z_val(P,i)
	reduction *= (Fp(denn))^(-1)*res_pol(z)

	return value+cua,reduction;
		
		
		
def simplify(self,beta,m):
	assert m >= 0
	assert beta in self.parent
	print"asdasds"
	p = self.prime
	den,exp,num = localize(beta,p)
	beta = self.parent(0)
	precision = ceil(m/self.rth_level().e)-exp
	if precision > 0:
		power = p^precision
		self.single_factor_lifting(precision*self.rth_level().e-self.rth_level().V)
		phi = coeff_mod(self.rth_level().phi, power)
		phi = num.parent()([i.numerator() for i in phi.list()])
		num = coeff_mod(coeff_mod(inverse_mod(den,power)*num , phi) ,power)
		beta = num(self.parent.0)*p^exp
	return beta


def long_embed(self,Fq,alpha):
	""""
    embeds alpha along a tower of finte fields to Fq
    """
	if alpha in Fq:
		return alpha
	for lvl in self.levels:
		gf_ext = lvl.gf_extension
		if alpha in gf_ext.base_field:
			alpha = Embed(gf_ext,alpha)	
			if alpha in Fq:
				return alpha

def long_shift_down(self,Fq,alpha):
	""""
    pushes alpha along a tower of finte fields down to Fq
    """
	if alpha in Fq:
		return alpha
	for lvl in reversed(self.levels):
		gf_ext = lvl.gf_extension
		if len(gf_ext.base_field) == len(gf_ext.ext_field):
			alpha = elt_seq(gf_ext,alpha)[0]	
			if alpha in Fq:
				return alpha



def is_separable(f):
	""""
    checks if f is separable 
    """
	f_prime = derivative(f)
	if gcd(f,f_prime) != 1:
		return False
	else:
		return True	





def reduction_mod_P(alpha,P):
	""""
    computes alpha mod P
    """
	
	val,red = P_valuation(alpha,P)
	assert val >= 0
	if val > 0:
		return P.rth_level().Fq(0)
	else:
		return red	



class Error(object):

    def __init__(self):
        self.error = []
        
def into_error(self,elt):
	self.error.append(elt)        



class factional_ideal(object):

    def __init__(self, p,F):
        self.parent = F
        self.polynomial_generator = p
