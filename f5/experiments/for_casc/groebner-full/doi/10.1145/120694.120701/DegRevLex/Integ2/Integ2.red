% Giovini, A., Mora, T., Niesi, G., Robbiano, L. and Traverso, C., "One sugar cube, please" or
% Selection strategies in the Buchberger algorithm, J. of the ACM 1991, pp. 49-54,
% https://doi.org/10.1145/120694.120701
%
% 3.7 Another example from integer programming

load_package groebner;

system := {-y**82*a + x**32*z**23,
   x**45 - y**13*z**21*b,
   y**33*z**12 - x**41*c,
   -y**33*z**12*d + x**22,
   x**5*y**17*z**22*e - 1,
   x*y*z*t - 1}$

vars := {e, d, c, b, a, t, z, y, x}$
torder(reverse vars, revgradlex)$

gb := groebner(system)$

end;
