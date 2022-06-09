% eco-7 system in revgradlex
% characteristic 0
% 0 dim

load_package f5;

system := {
        x1*x2*x7 + x1*x7 + x2*x3*x7 + x3*x4*x7 + x4*x5*x7 + x5*x6*x7 - 1,
        x1*x3*x7 + x2*x4*x7 + x2*x7 + x3*x5*x7 + x4*x6*x7 - 2,
        x1*x4*x7 + x2*x5*x7 + x3*x6*x7 + x3*x7 - 3,
        x1*x5*x7 + x2*x6*x7 + x4*x7 - 4,
        x1*x6*x7 + x5*x7 - 5,
        x6*x7 - 6,
        x1 + x2 + x3 + x4 + x5 + x6 + 1
}$

vars := {x1, x2, x3, x4, x5, x6, x7}$
torder(vars, revgradlex)$

gb := f5(system)$

end;
