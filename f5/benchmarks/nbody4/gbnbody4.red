% cyclic-5 system in lex

load_package groebner;

system := {
(b-d)*(bb-dd)-2*ff+2,
(b-d)*(bb+dd-2*ff)+2*(bb-dd),
(b-d)^2-2*(b+d)+f+1,
bb^2*b^3-1,dd^2*d^3-1,ff^2*f^3-1
}$

vars := {bb,dd,ff,b,d,f}$
torder(vars, revgradlex)$

share system;

lisp;
st := time();

algebraic;
gb := groebner(system)$

lisp;
prin2t({"TIME", time() - st});

end;
