% eco-10 system in revgradlex
% characteristic 0
% 0 dim

load_package groebner;

system := {
        x0*x1*x9+x1*x2*x9+x2*x3*x9+x3*x4*x9+x4*x5*x9+x5*x6*x9+x6*x7*x9+x7*x8*x9+x0*x9-1,
        x0*x2*x9+x1*x3*x9+x2*x4*x9+x3*x5*x9+x4*x6*x9+x5*x7*x9+x6*x8*x9+x1*x9-2,
        x0*x3*x9+x1*x4*x9+x2*x5*x9+x3*x6*x9+x4*x7*x9+x5*x8*x9+x2*x9-3,
        x0*x4*x9+x1*x5*x9+x2*x6*x9+x3*x7*x9+x4*x8*x9+x3*x9-4,
        x0*x5*x9+x1*x6*x9+x2*x7*x9+x3*x8*x9+x4*x9-5,
        x0*x6*x9+x1*x7*x9+x2*x8*x9+x5*x9-6,
        x0*x7*x9+x1*x8*x9+x6*x9-7,
        x0*x8*x9+x7*x9-8,
        x8*x9-9,
        x0+x1+x2+x3+x4+x5+x6+x7+x8+1
}$

vars := {x0, x1, x2, x3, x4, x5, x6, x7, x8, x9}$
torder(vars, revgradlex)$

gb := groebner(system)$

end;
