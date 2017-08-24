load("~/Dropbox/programming/Sage/Function_Fields/plusidealsFF.sage")

'''

##### Example 1

EE = Error()

Fq.<base> = GF(13)

K.<t> = FunctionField(Fq)

Ax.<x> =  K[]

p = t+2

n = 5

k = 7

r = 2

f = (x+sum([p^i for i in range(r+1)]))^n+p^k

F.<theta> = K.extension(f)

om_rep = montes(F,p)

g = GENUS(F)

print"g is",g





##################

EE = Error()

Fq.<base> = GF(13^2)

K.<t> = FunctionField(Fq)

Ax.<x> =  K[]

p = t^2+base*t+1

k = 7

f = (x^2-2*x+4)^3+p^k

F.<theta> = K.extension(f)

g = GENUS(F)

print"g is",g


########## Gutes Bsp



EE = Error()

Fq.<base> = GF(13^2)

K.<t> = FunctionField(Fq)

Ax.<x> =  K[]

p =  t^2+base*t+1

E_1 = x^2+p

E_2 = E_1^2+(p-1)*p^3*x

E_3 = E_2^3+p^11

E_4 = E_3^3+p^29*x*E_2

#E_5 = E_4^2+(p-1)*p^42*x*E_1*E_3^2

#E_6 = E_5^2+p^88*x*E_3*E_4

print"berechnet"

F.<theta> = K.extension(E_3)

FF, cf = infinity_representation(F)

om_inf = montes(FF,t)

print"local index at inf",om_inf[0].local_index

g = GENUS(F)

print"g is",g


###############
 
 
 
EE = Error()

Fq.<base> = GF(11^2)

K.<t> = FunctionField(Fq)

Ax.<x> =  K[]

p =  t+1#t^2+base*t+1

n = 2

l = 12

k = 11

f = sum([x^i for i in range(l)])+p^k

F.<theta> = K.extension(f)

#om_rep = montes(F,p)


FF,cf = infinity_representation(F)	

om_rep = montes(FF,t)

g = GENUS(F)

print"g is",g



######################


EE = Error()

Fq.<base> = GF(11)

K.<t> = FunctionField(Fq)

Ax.<x> =  K[]

p = t+1#t^2+base*t+1

k = 2

f = ((x^6+4*p*x^3+3*p^2*x^2+4*p^2   )^2+p^6)^3+p^k

F.<theta> = K.extension(f)

q  = t^8 + 8*t^7 + 6*t^6 + t^5 + t^4 + 10*t^2 + 7*t + 8

g = GENUS(F)

print"g is",g




######################


EE = Error()

Fq.<base> = GF(17^4)

K.<t> = FunctionField(Fq)

Ax.<x> =  K[]


f = x^9+x^5*(t^23+12*t^8)+t^123*x^4+(t^12+1)^12+2

F.<theta> = K.extension(f)


#p = t^6 + 12*base^3 + 11*base^2 + 15*base + 1

#om_rep = montes(F, p)


time g = GENUS(F)

print"g is",g

'''
###################### (Sage can not compute Disc(f))


EE = Error()

Fq.<base> = GF(5)

K.<t> = FunctionField(Fq)

Ax.<x> =  K[]

p = t

g = x^11+t^13+3*x^7+(2*(t+1)^8+2*t)^7*x^6+(2*t^8+t^6+2*t)^45*x^3+(t+1)^5*t^4*(2*x+x)+t^12+2
f = g^2+x^4*(t+2)^4-t^12

F.<theta> = K.extension(f)

FF,cf = infinity_representation(F) 

om_inf =montes(FF,t)

#time g = GENUS(F)

#print"g is",g
