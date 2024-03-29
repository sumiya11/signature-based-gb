% henrion-5 system in revgradlex
% characteristic 0
% 0 dim

load_package groebner;

system := {
2*f1*f2*f3*f4*groebner-9823275,
21/5*f1*f2*f4*groebner+16/5*f1*f3*f4*groebner+9/5*f2*f3*f4*groebner+24/5*f1*f2*f3*groebner+5*f4*f3*f1*f2-4465125,
14/5*f4*groebner*f1+14/5*f4*groebner*f2+8/5*f3*f4*groebner+18/5*f1*f2*groebner+24/5*f1*f3*groebner+18/5*f2*f3*groebner+4*f3*f1*f2+6*f1*f2*f4+6*f3*f4*f1+4*f2*f3*f4-441486,
7/5*f4*groebner+12/5*groebner*f1+12/5*groebner*f2+12/5*groebner*f3+3*f1*f2+4*f3*f1+4*f4*f1+3*f2*f3+4*f4*f2+3*f3*f4-15498,
6/5*groebner+2*f4+2*f3+2*f2+2*f1-215,
f1+2*f2+3*f3+4*f4+5*groebner+6*t
}$

vars := {f1,f2,f3,f4,groebner,t}$
torder(vars, revgradlex)$

gb := groebner(system)$

end;
