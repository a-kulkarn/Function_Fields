Attach("/Users/JB/Mathematik/Programming/Magma/Diss/Helpfunctions/Helpfunctions.m");
Attach("~/Dropbox/programming/Function Fields/+IdealsFF/+IdealsFF.m");

Z := Integers();
PutInZ([1]);
Z`Error := [**];


///////////// First Example!!!
/*


Fq := GF(13);

A<t> := PolynomialRing(Fq);

Ax<x> :=  PolynomialRing(A);

p :=  t + 2;

n := 5;

k := 7;

r := 2;

f := (x+&+[p^i : i in [0..r]])^n+p^k;

F<theta> := FunctionField(f);



g := GENUS(F);

print"Genus of g =",g;


///////////// Second Example!!!





Fq := GF(13^2);

A<t> := PolynomialRing(Fq);

Ax<x> :=  PolynomialRing(A);

p := t^2+Fq.1*t+1;

n := 5;

k := 7;

r := 2;

f := (x^2-2*x+4)^3+p^k;



F<theta> := FunctionField(f);

Montes(F,p);

g := GENUS(F);

print"Genus of g =",g;

//########## Gutes Bsp




Fq<base> := GF(13^2);

A<t> := PolynomialRing(Fq);

Ax<x> :=  PolynomialRing(A);

p := t^2+base*t+1;

E_1 := x^2+p;

E_2 := E_1^2+(p-1)*p^3*x;

E_3 := E_2^3+p^11;

E_4 := E_3^3+p^29*x*E_2;

E_5 := E_4^2+(p-1)*p^42*x*E_1*E_3^2;

E_6 := E_5^2+p^88*x*E_3*E_4;

F<theta> := FunctionField(E_3);


g := GENUS(F);

print"Genus of g =",g;


///////////// 


Fq<base> := GF(11^2);

A<t> := PolynomialRing(Fq);

Ax<x> :=  PolynomialRing(A);

p := t+1;//t^2+base*t+1;

n := 2;

l := 12;

k := 11;

f := &+[x^i : i in [0..l-1]]+p^k;

F<theta> := FunctionField(f);


Montes(F,p);


g := GENUS(F);

print"Genus of g =",g;



///////////// Second Example!!!




Z := Integers();

Z`Error := [**];

Fq := GF(11);

A<t> := PolynomialRing(Fq);

Ax<x> :=  PolynomialRing(A);

p := t+1;

k := 2;

f := ((x^6+4*p*x^3+3*p^2*x^2+4*p^2   )^2+p^6)^3+p^k;


F<theta> := FunctionField(f);

q := t^8 + 8*t^7 + 6*t^6 + t^5 + t^4 + 10*t^2 + 7*t + 8;



g := GENUS(F);

print"Genus of g =",g;





///////////// New Examples:!!!


Fp:=FiniteField(17,4);

Fpt<t>:=PolynomialRing(Fp);
Kx<T> := RationalFunctionField(Fp);

KxT<x> := PolynomialAlgebra(Kx);


f:=x^9+x^5*(t^23+12*t^8)+t^123*x^4+(t^12+1)^12+2;

F<y>:= FunctionField( f);

g := GENUS(F);

print"Genus of g =",g;


*/
/////////////

Fp:=FiniteField(5);

Fpt<t>:=PolynomialRing(Fp);
Kx<T> := RationalFunctionField(Fp);

KxT<x> := PolynomialAlgebra(Kx);


g:=x^11+t^13+3*x^7+(2*(t+1)^8+2*t)^7*x^6+(2*t^8+t^6+2*t)^45*x^3+(t+1)^5*t^4*(2*x+x)+t^12+2;
f:=g^2+x^4*(t+2)^4-t^12;
F<y>:= FunctionField( f);

g := GENUS(F);

print"Genus of g =",g;
