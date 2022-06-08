% Giovini, A., Mora, T., Niesi, G., Robbiano, L. and Traverso, C., "One sugar cube, please" or
% Selection strategies in the Buchberger algorithm, J. of the ACM 1991, pp. 49-54,
% https://doi.org/10.1145/120694.120701
%
% 3.6 An example from integer programming

load_package f5;

system := {x**2*y*z**4 - t,
   x**5*y**7 - z**2*u,
   -x**3*z*v + y**2,
   -z**5*w + x*y**3}$

vars := {w, v, u, t, z, y, x}$
torder(reverse vars, revgradlex)$

gb := f5(system)$

end;
