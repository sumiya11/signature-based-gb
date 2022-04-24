% cyclic-5 system in lex

load_package f5;

system := {
3*(b-d)*(bb-dd)+(bb+dd)-6*ff+4,
3*(b-d)*(bb+dd-2*ff)+5*(bb-dd),
3*(b-d)^2-6*(b+d)+4*f+3,
bb^2*b^3-1,dd^2*d^3-1,ff^2*f^3-1
}$

vars := {bb,dd,ff,b,d,f}$

share system, vars;

lisp;
st := time();

algebraic;
gb := f5(system, vars, 'revgradlex)$

lisp;
prin2t({"TIME", time() - st});

end;
