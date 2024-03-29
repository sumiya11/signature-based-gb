% cyclic-5 system in lex

load_package groebner;

system := {
15*b**4*c*d**2+6*b**4*c**3+21*b**4*c**2*d-144*b**2*c-8*b**2*c**2*e
-28*b**2*c*d*e-648*b**2*d+36*b**2*d**2*e+9*b**4*d**3-120,

30*c**3*b**4*d-32*d*e**2*c-720*d*b**2*c-24*c**3*b**2*e-432*c**2*b**2
+576*e*c-576*d*e+16*c*b**2*d**2*e+16*d**2*e**2+16*e**2*c**2+9*c**4*b**4+5184
+39*d**2*b**4*c**2+18*d**3*b**4*c-432*d**2*b**2+24*d**3*b**2*e-16*c**2*b**2*d*e
-240*c,

216*d*b**2*c-162*d**2*b**2-81*c**2*b**2+5184+1008*e*c-1008*d*e
+15*c**2*b**2*d*e-15*c**3*b**2*e-80*d*e**2*c+40*d**2*e**2+40*e**2*c**2,

261+4*d*b**2*c-3*d**2*b**2-4*c**2*b**2+22*e*c-22*d*e
}$

vars := {b, c, d, e}$
torder(vars, revgradlex)$

share system;

lisp;
st := time();

algebraic;
gb := groebner(system)$

lisp;
prin2t({"TIME", time() - st});

end;
