% cyclic-6 system in revgradlex
% characteristic 0
% 0 dim
%

load_package groebner;

system := {
        z1 + z2 + z3 + z4 + z5 + z6,
        z1*z2 + z1*z6 + z2*z3 + z3*z4 + z4*z5 + z5*z6,
        z1*z2*z3 + z1*z2*z6 + z1*z5*z6 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6,
        z1*z2*z3*z4 + z1*z2*z3*z6 + z1*z2*z5*z6 + z1*z4*z5*z6 + z2*z3*z4*z5 + z3*z4*z5*z6,
        z1*z2*z3*z4*z5 + z1*z2*z3*z4*z6 + z1*z2*z3*z5*z6 + z1*z2*z4*z5*z6 + z1*z3*z4*z5*z6 + z2*z3*z4*z5*z6,
        z1*z2*z3*z4*z5*z6 - 1
}$

vars := {z1,z2,z3,z4,z5,z6}$
torder(vars, revgradlex)$

gb := groebner(system)$

end;
