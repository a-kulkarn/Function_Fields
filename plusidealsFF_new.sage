# -*- coding: utf-8 -*-

import pprint
from itertools import product
pp = pprint.PrettyPrinter(indent=4)

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

    def __init__(self, K, p, varphi, omega, initial=True):
        self.parent = K
        self.pol = K.defining_polynomial()
        self.prime = p
#        self.local_generator = p
        self.varphi = varphi
        self.levels = [ ]
        self.local_index = 0
        self.log_lg = 0
        self.sfl = [0, 0, 0, 0] # Single Factor Lifting
        self.phiadic = [ZZ[x](0) for i in range(0, 4)]
#        print"\n\n, check ",varphi
        if initial is True:
            self.add_initial_level(varphi, omega, p)
    
    def copy(self):
        copy = MontesType(self.parent, self.prime, self.varphi,
                          self.levels[0].omega, initial=False)

        for level in self.levels:
            copy.levels.append(level.copy())
        
        return copy

    def add_initial_level(self, varphi, omega, p):
        phi = self.pol.parent()(varphi)

        new_level = MontesTypeLevel(phi, omega, p)
#        new_level.f = varphi.degree()
        new_level.prod_f = varphi.degree()
        Fp = varphi.parent().base_ring()
        if varphi.degree() >1:
	        new_level.Fq = Fp.extension(varphi, name='z0', prefix='ww')
	        new_level.z = new_level.Fq.0
        else:
	        new_level.Fq = Fp
	        new_level.z = -varphi.list()[0]

        #new_level.Fq = FiniteField(p^varphi.degree(), name='z0',
        #                           conway=True, prefix='ww')

        new_level.Fqy = PolynomialRing(new_level.Fq, 'y0')
        
       

        self.levels.append(new_level)

    def add_new_level(self, phi, omega, u):
        s = len(self.levels)
        lvl_s = self.levels[-1]
        new_level = MontesTypeLevel(phi, omega, self.prime)

        new_level.V = lvl_s.e*u
        new_level.prod_e = lvl_s.prod_e * lvl_s.e
        new_level.prod_f = lvl_s.prod_f * lvl_s.f
        if lvl_s.f > 1:
        	new_level.Fq = lvl_s.Fq.extension(lvl_s.res_pol, name='w'+str(s), prefix='ww')
        else:

        	new_level.Fq = lvl_s.Fq
        	
        #new_level.Fq = FiniteField(self.prime^new_level.prod_f, name='w'+str(s),
         #                          conway=True, prefix='ww')
        
        
        
        
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
#        print"in refi with",self.rth_level().log_phi
        if more_factors is True:
            self.lvl(r).refinements.append([self.lvl(r).phi,
                                            self.lvl(r).slope])
        self.lvl(r).cutting_slope = ZZ(self.lvl(r+1).slope)
        self.lvl(r).phi = self.lvl(r+1).phi
        self.lvl(r).omega = self.lvl(r+1).omega
#        print"in refinement: cuttingslope",self.lvl(r).cutting_slope
        self.remove_last_level()

    def add_last_level(self, phiadic, side, side_dev, last_psi=True):
        r = len(self.levels)
        lvl_r = self.lvl(r)
        lvl_r.e = 1

        if r > 1:
            # FIXME: What does nur stand for?
            nur = sum([self.lvl(j).slope/self.lvl(j).prod_e for j in [1..r-1]])
            self.sfl[0] = floor((lvl_r.V/lvl_r.prod_e)-nur)

        if side.p1.x == 0:
            slope = -side.slope

            lvl_r.h = ZZ(slope)
            # FIXME: Why do we just take the first 2 elements?
            self.phiadic[0:2] = phiadic[0:2]

            if last_psi:

                res_pol = self.residual_polynomial(r, side_dev)
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
#        print"Input phi_newton_polygon", i , phiadic	
#        print"\n\n\n ________________"
        n = 0
        sides = []
        cloud = []
        side_devs = []
        all_devs = []
        for k in [0..len(phiadic)-1]:
            val = 0
            dev = [ ]

            val, dev = self.value(i, phiadic[k])
#            print "N_%d : v_%d(%s phi_%d) = %s + %d" % (i, i, phiadic[k], i, str(val), n)
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
#        print"in value",i,a,a.parent()
 #       print"\n\n\n ____________________________"
        val = +Infinity
        if a == 0:
            if i == 1:
                devs = 0
            else:
                devs = []
            return val, devs

        if i == 1:
            val = min([valuation(c, self.prime) for c in list(a)])
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

            while quot != 0 and min_height <= val:
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

        assert i <= len(self.levels)

        lvl_i = self.lvl(i)

        height = side_devs[-1][1]
        res_coeffs = [ ]
        for dev in side_devs[:-1]:

            if (i == 1 and dev == 0) or (i > 1 and len(dev) == 0): 
                res_coeffs.append(lvl_i.Fq(0))
            elif i == 1:

                # coefficients are polynomials too.
                coeff = ZZ[x](dev) // self.prime^(height)


                res_coeffs.append(coeff(lvl_i.z))
                j = side_devs.index(dev)
            else:
                lvl_im1 = self.lvl(i-1)
                lprime = (1 - lvl_im1.inv_h*lvl_im1.h) // lvl_im1.e
 
                coeff = self.residual_polynomial(i-1, dev)

                lifted_coeff = lvl_i.Fqy([lvl_i.Fq(c) for c in coeff])
                twist_exp = lprime*dev[-1][0]-lvl_im1.inv_h*dev[-1][1]

                res_coeffs.append(lvl_i.z^(twist_exp) * lifted_coeff(lvl_i.z))

                j = side_devs.index(dev)

                                       
            height = height - lvl_i.h

        res_pol = lvl_i.Fqy(res_coeffs)

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
        u = ef * lvl_s.V #+lvl_s.f*lvl_s.h
        if s > 1:
            txp = -self.lvl(s-1).inv_h * (u // self.lvl(s-1).e)
            twist = lvl_s.z^(txp)
        else:
            twist = lvl_s.Fq(1)

        # Reductum from magma
        res_pol = lvl_s.Fqy(list(lvl_s.res_pol)[:-1])
        u += lvl_s.f * lvl_s.h
        phi0 = self.construct(s, res_pol, 0, u)

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
        var = lvl_i.phi^(lvl_i.e)
        phi0 = 0
        height = u - res_pol.degree()*lvl_i.h


        if i == 1:
            for a in reversed(list(res_pol)):
                lift = ZZ[x](list(a.polynomial()))
                phi0 = (phi0 * var) + (lift * self.prime^height)
                height = height + lvl_i.h

        else:
            step = (lvl_i.e * lvl_i.V) + lvl_i.h
            new_V = u - (res_pol.degree() * step) - (s * lvl_i.V)
            lvl_im1 = self.lvl(i-1)

            for a in reversed(list(res_pol)):
                # FIXME: What does pj stand for?
                pj = 0
                if a != 0:

                    txp, s_im1 = divmod(lvl_im1.inv_h*new_V, lvl_im1.e)
                    u_im1 = (new_V - (s_im1 * lvl_im1.h)) // lvl_im1.e

                    c = (a*lvl_i.z^txp)

                    if lvl_im1.f>1 :# c.parent().base_ring() == lvl_im1.Fq:
                        eltseq = list(c.polynomial())

                    else:
                        eltseq = c
                        w_im1 = lvl_im1.Fq.gen()
                        prod_f = lvl_im1.prod_f

            #            m = vector(c) * lvl_i.inv_M
 #                       print"asdasf 2",i,self.lvl(1).inv_M
            #            eltseq = [lvl_im1.Fq(m[j*prod_f:(j+1)*prod_f])
            #                       for j in range(len(m)/prod_f)]
            #            if c != sum([eltseq[j]*lvl_i.z^j for j in range(len(eltseq))]):
             #               raise(Exception,
             #                    "Eltseq calculated incorrectly for {0}"\
             #                            .format(c))


#                    into_error(EE,[lvl_im1,c])
                    new_res_pol = lvl_im1.Fqy(eltseq) #lvl_im1.Fqy([c])

                    pj = self.construct(i-1, new_res_pol, s_im1, u_im1)
                phi0 = (phi0 * var) + pj
                new_V += step
                height += lvl_i.h

        phi0 = phi0 * lvl_i.phi^(s)

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



        pmode = 'terse'
        p = self.prime
        exponent = self.sfl[0]
        nu = self.sfl[1]
        x0prec = self.sfl[2]
        x0num = self.phiadic[3]
        x0den = self.sfl[3]

        prod_e = lvl_r.prod_e
        h = lvl_r.h - lvl_r.cutting_slope
 #       print"in SFL: h",lvl_r.h , lvl_r.cutting_slope
        last_h = slope - lvl_r.cutting_slope
#        print"in SFL: last_h",last_h
        V = lvl_r.V + lvl_r.cutting_slope

        zp = Zp(p, nu + exponent + ceil((V+last_h)/prod_e), type='capped-abs',
                                                            print_mode=pmode)
        zpx = PolynomialRing(zp, 'X0')
        pol_zp = zpx(self.pol)
        psinum_zp = zpx(self.phiadic[2])

        path = path_of_precision(last_h, h)
        short_path = path_of_precision(h, x0prec)
 #       print"In SFL: path",path
#        print"with input",last_h,h
#        print"In SFL: shortpath",short_path
#        print"with input",h,x0prec
        zp = Zp(p, nu + exponent + ceil(h/prod_e), type='capped-abs',
                                              print_mode=pmode)

        zpx = PolynomialRing(zp, 'X')
        a1 = zpx(self.phiadic[1])

        zq = Zp(p, nu + exponent + ceil((V+path[1])/prod_e), type='capped-abs',
                                                        print_mode=pmode)
        zqt.<t> = PolynomialRing(zq)
        phi = zqt(lvl_r.phi)
        psinum = zqt(psinum_zp)
#        print"SFL 1"

        a0num, a0den = self.cancel((zqt(self.phiadic[0])*psinum) % phi, nu)
#print"SFL 1.1"
        a1num, a1den = self.cancel((zqt(a1)*psinum) % phi, nu)
        for i in range(1,len(short_path)):#range should start at 2 maybe?
            low_precision = a1den + 2 * x0den  + ceil(short_path[i]/prod_e)
            x0num, x0den = self.inversion_loop([ a1num, a1den], x0num, x0den,
                                               phi, low_precision)
        anum, aden = self.cancel((a0num*zqt(x0num)) % phi, x0den+a0den)
        phi = phi + anum
  #      print"in SFL_ after shortpath: phi",ZZ[x](phi)
        for i in range(1, len(path)-1):
            loop_prec = nu + exponent + ceil((V+path[i+1])/prod_e)
          
            zq = Zp(p, loop_prec, type='capped-abs', print_mode=pmode)
            zqt.<t> = PolynomialRing(zq)
            phi = change_precision(zqt(phi), loop_prec)
            psinum = zqt(psinum_zp)
            qq, c0 = zqt(pol_zp).quo_rem(phi)

            c1 = qq % phi

         
            c0num, c0den = self.cancel((c0*psinum) % phi, nu)
 #           print"c0num, c0den",c0num, c0den
            c1num, c1den = self.cancel((c1*psinum) % phi, nu)

            low_precision = c1den + 2 * x0den + ceil(path[i]/prod_e)
            x0num, x0den = self.inversion_loop([c1num, c1den], x0num, x0den,
                                               phi, low_precision)
            xnum = change_precision(zqt(x0num), low_precision)
#            print"in SFL: input for cnum",c0num , zqt(x0num), phi, x0den
            cnum, cden = self.cancel((c0num * zqt(ZZ[x](x0num))) % phi, x0den)

       #     print"in SFL: cnum",ZZ[x](cnum)		
            phi = phi + cnum
        #    print"in SFL: phi in path loop",ZZ[x](phi)		
#            print"_____________________________"
 
        self.sfl[2] = max(h, path[-2])
        lvl_r.phi = lvl_r.phi.parent()([Integers()(i) for i in phi.list()])
                      
        self.phiadic[3] = ZZ[x](x0num)
        self.sfl[3] = x0den
#        print"phi",phi
        
    def sfl_init(self):
        p = self.prime
        lvl_r = self.rth_level()
        print"so weit"
        prod_e = lvl_r.prod_e
        a1 = self.phiadic[1] # FIXME: Why [1]?
        psinum = ZZ[x](1)
#        print"aha 1"
        r = len(self.levels) - 1
 #       print"aha 2",r
        if r == 0:
#            print"r=0",min([valuation(a, p) for a in a1])
        	
            nu = min([valuation(a, p) for a in a1])

            # Evaluate a1/p^nu in z_1 (this may be z_0)
            klass = ZZ[x](a1 // p^nu)(self.lvl(1).z)
            
        else:
            val, dev = self.value(r+1, a1)
            res_pol = self.residual_polynomial(r, dev)

            logpsi = 0
            qq, s = (-val).quo_rem(prod_e)
            psinum, logpsi = self.prescribed_value(s)

            nu = -logpsi[0] - qq

            vector = dev[-1][0] * self.lvl(r).log_phi + dev[-1][1]*self.lvl(r).log_pi
            #print"logpsi + vector",logpsi + vector
#            into_error(EE,[self,logpsi + vector])
            klass = self.convert_logs(logpsi + vector)

            klass *= res_pol(self.lvl(r+1).z)

        self.phiadic[2] = psinum
        self.sfl[1] = nu
        self.sfl[2] = 1

        x0num = 0
        x0den = 0
#        print"aha 5.0,local_lift",klass^-1 
#        into_error(EE,[self,klass^-1])
        x0num, x0den = self.local_lift(klass^(-1))
#        print"aha 5" 
        self.phiadic[3] = x0num
        self.sfl[3] = x0den

    def prescribed_value(self, value): #log_lg or log not necessary as an input since it gets initialized
        """
        From +Ideals:
        If we are attached to the prime ideal P with Okutsu depth r, then
        logpsi=[a_0, ..., a_r] and psi=phi_1^a_1 ... phi_r^a_r, with
        v_P(p^a_0 psi(theta))=value.
        """
        psi = ZZ[x](1)
        r = len(self.levels)
        logpsi = (ZZ^r)(0)
#        print"logpsi",logpsi
        qq, val = value.quo_rem(self.rth_level().prod_e)
        logpsi[0] = qq
        if val > 0:
            body = val

#            into_error(EE,self)
            for k in reversed(range(r-1)):
                jj = (self.levels[k].inv_h * body) % self.levels[k].e
                logpsi[k+1] = jj
                psi = psi * self.levels[k].phi^jj
                res = (body - jj*self.levels[k].h) // self.levels[k].e
                body = res - jj*self.levels[k].V
            logpsi[0] += res

        return psi, logpsi
        
    def SFL(self, slope):
        

        self.single_factor_lifting(slope)
#        print"in SFL: before UpdateLastLevel",self.rth_level().phi
        self.update_last_level()        
        
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
            numlift = ZZ[x](clss.polynomial())
            denlift = 0
        else:
            #raise(NotImplemented, "local lift not implemented for i > 1.")
            expden = ceil(self.lvl(i).V/self.lvl(i).prod_e)			
            V = self.lvl(i).prod_e*expden
            log = V*self.lvl(i).log_pi
            log = vector(log.list()[:-1])
            newclss = self.convert_logs(log)    		
            H = V // self.lvl(i-1).e
            elt = self.lvl(i).z^(self.lvl(i-1).inv_h*H)*clss*newclss^(-1)
            c = self.lvl(i).Fq(elt)
#            print"it works?"
            varphi = self.lvl(i-1).Fqy(c.polynomial().list())			
            lift = self.construct(i-1, varphi, 0, H)
            v1lift = min([rr.valuation(self.lvl(1).prime) for rr in lift.coefficients()])
            numlift = lift // self.lvl(1).prime^v1lift
            denlift = expden-v1lift


				
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

        cancel = min([den] + [ a.valuation(self.prime) for a in poly])
        zq = poly[0].parent()
#        print "poly, p^cancel: {0}, p^{1}".format(poly, cancel)
        num = poly.parent()([ ZZ(c) / self.prime^cancel for c in poly ])
#        print "num: {0}".format(num)
        #num = poly / self.prime^cancel

        return num, den-cancel

    def inversion_loop(self, A, xnum, xden, phi, precision):
        anum = A[0]
        aden = A[1]

        zq = Zp(self.prime, precision, type='capped-abs', print_mode='terse')
        zqt.<t> = PolynomialRing(zq)
        phip = zqt(phi)
        xnum = zqt(ZZ[x](xnum))
#        print"in Inversion_loop",zqt(xnum),zqt,zqt(ZZ[x](xnum))
 #       print"data1"
#        into_error(EE,[zqt(anum),xnum])
        x1num, x1den = self.cancel(2*self.prime^(xden+aden) - (zqt(anum)*zqt(ZZ[x](xnum))) % phip, xden+aden)
#        print"data2"
 #       print"After Cancel x1num,x1den",x1num,x1den
        xnum, xden = self.cancel((xnum*x1num) % phip, xden+x1den)
 #       print" \n\n\n\n output invloop \n\n\n\n",xnum

        return xnum, xden
        
    def update_last_level(self):
        
        pol = self.pol

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

def montes(K, p, basis=False):

    if p.is_prime is False:
        raise Exception, "Error: p must be prime."
    f = K.defining_polynomial()
    if f.is_monic() is False:
        raise Exception, "Error: def.pol. of K must be monic."
    # TODO: Add requirements that coefficients of f are integers


    reps_OM = [ ]
    trees = [ ]
    total_index=0
    res_pol = PolynomialRing(FiniteField(p), 'y0')(f)

    for factor in res_pol.factor():

        tt = MontesType(K, p, factor[0], factor[1])
        tree, tree_index = montes_main_loop(K, p, tt)
        total_index += tree_index
        reps_OM += tree
        trees.append(tree)

    if len(reps_OM) == 1:
        reps_OM[0].rth_level().phi = f
        reps_OM[0].rth_level().slope = +Infinity

	## Produce local generators
    for P in reps_OM:
		l_gen,log = P.prescribed_value(1)		#element_with_prescribed_value(P,1)
		P.local_generator = l_gen(K.0)*p^log[0]
		P.log_lg = log				
		P.lvl(1).local_index = total_index

    return reps_OM


def element_with_prescribed_value(P,val):
	"""
	Produces element a with V_P(a)=val
	"""

	tmp = P.prescribed_value(val)
	print"Psi,logLG",tmp
	if ramification_index(P) ==1:
		return P.parent(P.prime),tmp[1]

		

	r = len(P.levels)
	psi = prod([P.lvl(i).phi^tmp[1][i-1] for i in range(1,r+1) ])
	return psi(P.parent.0),tmp[1]
	 



def montes_main_loop(K, p, tt):
    print"tt",tt, len(tt.levels)
    total_index = 0
    leaves = [ ]
    run = 0
    type_stack = [ tt ]
    while len(type_stack) > 0:
        tt = type_stack.pop()
        r = len(tt.levels)
        lvl_r = tt.rth_level()

        print"Input of phi_expansion",lvl_r.phi,lvl_r.omega
        print"_________________________ \n\n\n"
        phiadic, quotients = phi_expansion(tt.pol,
                                           lvl_r.phi,
                                           lvl_r.omega)
        print"\n\n ****************************** \n\n"                                           
        ## Phi-Newton Polygon
        sides, side_devs = tt.phi_newton_polygon(r, phiadic)     
 #       print"outout phi_newton_polygon",        sides, side_devs
#        print"_________________________ \n\n\n"

        length_N = lvl_r.omega
        index_N = -lvl_r.cutting_slope * (length_N*(length_N-1) // 2)
        starting_prod_f = lvl_r.prod_f
        if length_N == 1:
 
            tt.add_last_level(phiadic, sides[0], side_devs[0])
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

        else:
            side_indexes_types = [ ]

        previous_h = 0
        for i, tt in side_indexes_types:
            #for i in reversed([0..len(sides)-1]):
            lvl_r = tt.rth_level()
            side = sides[i]


            lvl_r.h = - side.slope.numerator()
            lvl_r.e = side.slope.denominator()
            lvl_r.slope = - side.slope
            lvl_r.inv_h = inverse_mod(lvl_r.h, lvl_r.e)
            
            lprime = (lvl_r.inv_h*lvl_r.h - 1) // lvl_r.e

            new_pi = list( lvl_r.inv_h * lvl_r.log_phi - lprime*lvl_r.log_pi )
            new_pi.append(0)

            lvl_r.log_gamma = lvl_r.e*lvl_r.log_phi - lvl_r.h*lvl_r.log_pi

            E = ZZ(side.width())
            H = ZZ(side.height())
            index_N += (E * previous_h) + ((E*H-E-H+(E // lvl_r.e)) // 2)

            previous_h += H
            res_pol = tt.residual_polynomial(r, side_devs[i])
            print"this is res_pol", res_pol

            #fix for multiqudratic fields: In general one can not do that because the residual degree might be bigger than 2
      
            res_pol = res_pol.monic()
            factors = res_pol.factor()


            
            factors_types = [ (factors[0], tt) ]

            if len(factors) > 1:
                factors_types += [ (factors[i], tt.copy()) for i in [1..len(factors)-1] ]
                more_factors = True

            for factor, tt in factors_types:
                lvl_r = tt.rth_level()

                omega = factor[1]
                lvl_r.res_pol = factor[0]
                lvl_r.f = lvl_r.res_pol.degree()

                ## Representative
                # The new level (lvl_rp1) is already part of the type.
                if r == 4: 
                	into_error(EE, [tt,omega])
#                	return
                lvl_rp1 = tt.representative(omega)

                if lvl_r.phi.degree() == lvl_rp1.phi.degree():
                    # Non-optimal, refining level r
                    tt.refinement(more_factors=more_factors)
                    tt.rth_level().cutting_slope = -ZZ(side.slope)
                else:
                    # Proceeding to higher order
                    
                    #tt.log_pi = vector(new_pi)
                    tt.rth_level().log_pi = vector(new_pi)

                    #print"log_pi",    tt.log_pi, new_pi ,tt.rth_level().log_pi               
                    #tt.log_phi = -(lvl_rp1.V // lvl_r.e) * lvl_r.log_pi	
                    
                    
                    
                    #tt.rth_level().log_phi = -(lvl_rp1.V // lvl_r.e) * lvl_r.log_pi
                    #tt.rth_level().log_phi = vector(new_pi)
                    
                    vect =-tt.rth_level().V*tt.rth_level().log_pi
                    vect[r+1]=1
                    tt.rth_level().log_phi=vect 
                    
                    #print"rth_level().log_phi data", lvl_rp1.V , lvl_r.e, lvl_r.log_pi
                    #ta[r+1]`logPi:=Vector(newPi);
                    #tt.log_phi = vector(list(tt.log_phi) + [1])
                   # tt.rth_level().log_phi = vector(list(tt.rth_level().log_phi) + [1])
#                    print"log_phi,log_pi",tt.log_phi,tt.log_pi					
#                    print"\n ________________\n",tt					                    	





                # Push the new or refined type onto the stack.
#                print"starting for loop",tt.rth_level().log_phi                
                type_stack.append(tt)

            ## End of `factors' for loop
        ## End of `sides' for loop
        
        total_index += starting_prod_f * index_N
#        print "Added %d * %d to the index (--> %d)" % (starting_prod_f, index_N, total_index,)

    return leaves, total_index

def phi_expansion(f, phi, omega):
    q = f
    coeffs = [ ]
    quos = [ ]
    for j in [0..omega]:
        q, r = q.quo_rem(phi)
        coeffs.append(ZZ[x](r))
        quos.append(ZZ[x](q))

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




def den(alpha):
	return lcm([i.denominator() for i in alpha.list()])
	
def num(alpha):
	return alpha*den(alpha)



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
	
	
def localize(alpha,p):
	A = p.parent()

	if alpha == 0:		
		return 1,0,A(0)

	numm = [ZZ(j)  for j in (num(alpha)).list()]		#[j.numerator()  for j in (num(alpha)).list()]

	valnum = min([x.valuation(p)  for x in numm])
	denom = A(den(alpha))

	valden = denom.valuation(p)
	denn = A(denom/p^valden)
	Ax.<x> = PolynomialRing(A)
	return	denn,valnum-valden,Ax(numm).quo_rem(p^valnum)[0]



def Z_val(self,i):

#	print"in Z_val",i, len(self.levels)
	if i == len(self.levels):
		#respol = self.lvl(i).res_pol
		z = -list(self.lvl(i).res_pol)[0]
	else:
		z = self.lvl(i+1).z

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
	print"Fpy",Fpy
	reduction = Fp(0)
	if alpha == 0:
		return Infinity, reduction
#	print"in P_val: in localize",alpha,p	
	denn, exp, numPol = localize(alpha,p)
#	print"in PVal:den,exp,numPol",denn,exp,numPol
	cua = exp*ramification_index(P)
#	print"in PVal:cua",cua

	pphi =  [Fp(i.numerator() % p) for i in list(P.lvl(1).phi)]	
	pnumPol = [Fp(i.numerator() % p) for i in list(numPol)]

	if vvaluation(Fpy(pnumPol),Fpy(pphi)) == 0:

#		print"reduction",reduction
		reduction = P.convert_logs(-cua*P.log_lg)
#		print"reduction",reduction
#		print"input reduction",Fp(denn)^(-1),pnumPol(P.lvl(1).z)
		reduction *= (Fp(denn))^(-1)*Fpy(pnumPol)(P.lvl(1).z)
#		print"reduction",reduction
		return cua,	reduction
	
	res_pol = 0
	z = 0
	side_devs = []
	val = 0
	value = 0
	i = 0

	while value == 0:
#		print"in while loop with i=",i
		if i < r:
			i +=1	

		else:

			P.SFL(2*P.rth_level().h)	
			print"in PValuation, SFL with 2*P.rth_level().h", P.rth_level().phi,P.rth_level().slope


		val, side_devs = P.value(i+1,numPol)
	

		res_pol = P.residual_polynomial(i, side_devs)
		if vvaluation(res_pol,P.lvl(i).res_pol) == 0:
			value = val*(P.lvl(r).prod_e).quo_rem((P.lvl(i).e*P.lvl(i).prod_e))[0]
		
# if RED:

	log = side_devs[len(side_devs)-1][0]*P.lvl(i).log_phi+ side_devs[len(side_devs)-1][1]*P.lvl(i).log_pi
	P.log_lg, log = equalize_logs(P.log_lg,log)	#vector([P.rth_level().log_gamma[1]])


# 	print"in PVal: before ConvertLogs",log-(value+cua)*P.log_lg
 	reduction = P.convert_logs(log-(value+cua)*P.log_lg)
 #	print"in PVal: after ConvertLogs",reduction
	z = Z_val(P,i)
#	print"In P_val : output Z_val", z
#	print"In Pvsl: reduction data",reduction
	reduction *= (Fp(denn))^(-1)*res_pol(z)

	return value+cua,reduction;
	




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


def inertia_degree(self):
	return self.rth_level().prod_f

def ramification_index(self):
	return self.rth_level().prod_e


def find_degree_less2_prime(K,primes,n):
	""""
    Finds prime ideal with inertia degree <= 2
    """
	p = n.next_prime()
	inertia_deg = 3
	used_primes = []
	prime = 0
	while inertia_deg > 2:
		om_rep = montes(K,p)[0]

		for P in om_rep:
			inertia_deg = inertia_degree(P)

			if inertia_deg <= 2:
				om_prime = om_rep
				break
		
		used_primes.append(p)
		p = (p+1).next_prime()
		while p in used_primes:
			p = (p+1).next_prime()

	return om_prime


def find_prime_ideals(K,primes,tries):
	""""
    Collects > tries many prime ideals of inertia degree <= 2
    """    
	prime_ideals = []

	start_prime =  random_prime(6) #should choose prime number which is not in primes 
									#and is small
	p_Primes = find_degree_less2_prime(K,primes,p-1)



	primes.append(p_Primes[0].prime)
	for P in p_Primes:
		if inertia_degree(P) <= 2:
			prime_ideals.append(P)
			
	run = len(prime_ideals)		
			
	while run < tries:
		p_Primes = find_degree_less2_prime(K,primes,prime_ideals[-1].prime)
		primes.append(p_Primes[0].prime)

		for P in p_Primes:
			if inertia_degree(P) <= 2:
				prime_ideals.append(P)
		run = len(prime_ideals)		
		
	return prime_ideals
	
	
def is_square(alpha,prime_ideals):
	""""
    Checks if alpha mod P is a square for all P in prime_ideals
    """    
	K = alpha.parent()
	for P in prime_ideals:

		alpha_P = reduction_mod_P(alpha,P)
		if not alpha_P == 0 and not alpha_P.is_square():

			return False
			
	return true	
		
		
def is_square_default(alpha,primes):
	""""
    Checks if alpha mod P is a square for sufficient many prime ideals P
    """    
	K = alpha.parent()

	tries = 2*K.degree()	

	prime_ideals = find_prime_ideals(K,primes,tries)

	return 	is_square(alpha,prime_ideals)


def reduction(alpha,P,m):
	"""
	Computes reduction map O_K-->O_K/P^m
	"""
	if m <= 0:
		raise Exception, "Error: Third argument should be positive."
	beta = alpha
	if beta == 1:
		return [P.rth_level().Fq(1)] +[P.rth_level().Fq(0) for j in [1..m-1]]
	K = alpha.parent()
#	print"initial beta",beta
	beta,P = shrink(beta,P,m)
#	print"after shrink",beta	
	value, red = P_valuation(beta,P)
#	print"after PValuation",value,red
	if value < 0:
		raise Exception, "Error: First argument should be P-integral."
	cla = [P.rth_level().Fq(0) for j in [1..m]]
	run = 0
	while value < m and run < 400:
		run +=1

		
		cla[value] = red
		if value == m-1:
#			print"drinne"
			value = m

		else:
#			print"*****************";
		
			num,den = P.local_lift(red)
			lo_li = K(num)/P.prime^den
#			print"LocalLift?",lo_li
#			print"!!!!!!!!!!!\n local_lift with",beta,P.local_generator,value
#			print"\n\n"			
			beta-=lo_li*P.local_generator^value
#			print"after LocalLift",beta
#			print"___________________"; 
#			print"in reduction: calling shrink with",beta
			beta,P = shrink(beta,P,m)
#			print"----^________^-------\n \n that is beta",beta
			value,red = P_valuation(beta,P)
#			print"value,red",value,red
			
	return cla
    
    
        
def shrink(beta,P,m):
	"""
	Replace beta by an element num/p^power congruent to beta mod P^m such that 
	num is given by a polynomial of degree < e_Pf_P, 
	The element beta is assumed (without checking) to be P-integral. 
	"""
	if beta == 0:
		return 0,P
	p = P.prime
#	print"In shrink: beta",beta
	den, exp,num = localize(beta,p)
#	print"In shrink: den, exp,num",den, exp,num
	beta = P.parent(0)
	e_P = ramification_index(P)
	precision = ceil(m/e_P)-exp
	R = P.rth_level().phi.parent()

	if precision > 0:		
		power = p^precision
#		print"In shrink: SFL with ",precision*e_P-P.rth_level().V;

		P.SFL(precision*e_P-P.rth_level().V)
		phi = R([i % power for i in list(P.rth_level().phi)])
#		print"In shrink: Input num",inverse_mod(den, power),num , phi, power
		num = R([j% power for j in list((inverse_mod(den, power)*num) % phi)])
#		print"In shrink: num",num
		beta = num(K.0)*p^exp
	
	return beta,P
	
	

class Error(object):

    def __init__(self):
        self.error = []
        
def into_error(self,elt):
	self.error.append(elt)        

def random_integer(bound):
	sig = ZZ.random_element(2)
	return (-1)^sig*ZZ.random_element(bound)

def random_element(K,bound):

	return K([random_integer(bound) for i in range(K.degree())])