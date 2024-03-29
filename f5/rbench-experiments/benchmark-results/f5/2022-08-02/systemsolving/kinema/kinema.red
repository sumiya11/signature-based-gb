% kinema system in revgradlex
% characteristic 0

load_package f5;
1;

system := {
    z1^2 - 12*z1 + z2^2 + z3^2 - 68,
    z4^2 + z5^2 - 12*z5 + z6^2 - 68,
    z7^2 + z8^2 - 24*z8 + z9^2 - 12*z9 + 100,
    z1*z4 - 6*z1 + z2*z5 + z3*z6 - 6*z5 - 52,
    z1*z7 - 6*z1 + z2*z8 + z3*z9 - 12*z8 - 6*z9 + 64,
    z4*z7 + z5*z8 - 6*z5 + z6*z9 - 12*z8 - 6*z9 + 32,
    2*z2 + 2*z3 - z4 - z5 - 2*z6 - z7 - z9 + 18,
    z1 + z2 + 2*z3 + 2*z4 + 2*z6 - 2*z7 + z8 - z9 - 38,
    z1 + z3 - 2*z4 + z5 - z6 + 2*z7 - 2*z8 + 8
}$

vars := {z1, z2, z3, z4, z5, z6, z7, z8, z9}$
torder(vars, revgradlex)$

gb := f5(system)$

end;
