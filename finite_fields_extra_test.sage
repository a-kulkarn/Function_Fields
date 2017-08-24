load("~/Dropbox/programming/Sage/Function_Fields/finite_fields_extra.sage")



K.<z0> = GF(11)
tmp.<t> = K[]
pol = t^6+6*t^4+9*t^2+1 
F.<z1> = GF(len(K)^pol.degree())
L_ext = GF_Extension(F, K, pol)

zz = L_ext.root_pol

elt_seq(L_ext,zz)
'''
magma:

K := GF(11);

tmp<t> := PolynomialRing(K);

pol := t^6+6*t^4+9*t^2+1 ;

F := ext<K|pol>;


'''





'''

K = F
tmp2.<t_1> = K[]

pol2 = t_1^2 + (z1^3 + 2*z1)*t_1 + 2*z1^2 + z1 + 1

F.<z2> = GF(len(K)^pol2.degree())
L_ext2 = GF_Extension(F, K, pol2)








Tow = GF_Tower(L_ext)

add_new_level(Tow,L_ext2)

alpha = F.0^12

test = change_rep(Tow, alpha)

print"1test",test


L_top = Tow.fields[-1]

elt = L_top.ext_field.random_element()

ww = lift(L_top,elt)

print"check",ww



alpha = F.0^12

test = change_rep(Tow, alpha)

print"1test",test

Tow.add_new_level(F,K,pol)



####### lvl 2

K = F
tmp.<t> = K[]
pol = t^3 + 2*t + 1
F = GF(len(K)^pol.degree())
Tow.add_new_level(F,K,pol)



#### isoliertes Bsp

K = GF(3^4)
tmp.<t> = K[]
pol = t^3 + 2*t + 1
F = GF(len(K)^pol.degree())

L = GF_Extension(F, K, pol)






check = true
elt = 0

for alpha in F:
	test = iso_to_res_field(L,alpha)
	print"alpga",alpha	
	if not iso_to_gf(L,test) == alpha:
		elt = alpha
		check = False
		break
		
		
print"__________________\n"
print" checki", check		


for i in range(1000):

	aa = I.random_element()
	print"i",i
	test = iso_to_gf(L,aa)

	if not iso_to_res_field(L,test) == aa:
		elt = alpha
		check = False
		break
		
		
print" checki", check		






aa = I.0
test = iso_to_gf(L,aa)
iso_to_res_field(L,test)
aa

'''


