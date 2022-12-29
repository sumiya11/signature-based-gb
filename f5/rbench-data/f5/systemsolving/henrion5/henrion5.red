% henrion-5 system in revgradlex
% characteristic 0
% 0 dim

load_package f5;

system := {
2*f1*f2*f3*f4*f88-9823275,
21/5*f1*f2*f4*f88+16/5*f1*f3*f4*f88+9/5*f2*f3*f4*f88+24/5*f1*f2*f3*f88+5*f4*f3*f1*f2-4465125,
14/5*f4*f88*f1+14/5*f4*f88*f2+8/5*f3*f4*f88+18/5*f1*f2*f88+24/5*f1*f3*f88+18/5*f2*f3*f88+4*f3*f1*f2+6*f1*f2*f4+6*f3*f4*f1+4*f2*f3*f4-441486,
7/5*f4*f88+12/5*f88*f1+12/5*f88*f2+12/5*f88*f3+3*f1*f2+4*f3*f1+4*f4*f1+3*f2*f3+4*f4*f2+3*f3*f4-15498,
6/5*f88+2*f4+2*f3+2*f2+2*f1-215,
f1+2*f2+3*f3+4*f4+5*f88+6*t
}$

vars := {f1,f2,f3,f4,f88,t}$
torder(vars, revgradlex)$

gb := f5(system)$

end;
