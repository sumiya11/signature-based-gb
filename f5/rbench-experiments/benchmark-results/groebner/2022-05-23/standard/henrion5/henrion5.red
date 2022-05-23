% cyclic-5 system in lex

load_package groebner;

system := {
2*f1*f2*f3*f4*f5-9823275,
21/5*f1*f2*f4*f5+16/5*f1*f3*f4*f5+9/5*f2*f3*f4*f5+24/5*f1*f2*f3*f5+5*f4*f3*f1*f2-4465125,
14/5*f4*f5*f1+14/5*f4*f5*f2+8/5*f3*f4*f5+18/5*f1*f2*f5+24/5*f1*f3*f5+18/5*f2*f3*f5+4*f3*f1*f2+6*f1*f2*f4+6*f3*f4*f1+4*f2*f3*f4-441486,
7/5*f4*f5+12/5*f5*f1+12/5*f5*f2+12/5*f5*f3+3*f1*f2+4*f3*f1+4*f4*f1+3*f2*f3+4*f4*f2+3*f3*f4-15498,
6/5*f5+2*f4+2*f3+2*f2+2*f1-215,
f1+2*f2+3*f3+4*f4+5*f5+6*t
};

vars := {f1,f2,f3,f4,f5,t};
torder(vars, revgradlex)$

gb := groebner(system)$

end;
