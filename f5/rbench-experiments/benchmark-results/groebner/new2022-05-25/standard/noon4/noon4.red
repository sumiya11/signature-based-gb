% noon-4 system in revgradlex
% characteristic 0
% 0 dim
%
% reference TODO

load_package groebner;

system := {
        10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 - 11*x1 + 10,
        10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 - 11*x2 + 10,
        10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 - 11*x3 + 10,
        10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 - 11*x4 + 10
}$

vars := {x1, x2, x3, x4}$
torder(vars, revgradlex)$

gb := groebner(system)$

end;
