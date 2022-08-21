% Giovini, A., Mora, T., Niesi, G., Robbiano, L. and Traverso, C., "One sugar cube, please" or
% Selection strategies in the Buchberger algorithm, J. of the ACM 1991, pp. 49-54,
% https://doi.org/10.1145/120694.120701
%
% 3.2 Cyclic roots

load_package f5;

system := {x + y + z + t + u,
   x*y + y*z + z*t + t*u + u*x,
   x*y*z + y*z*t + z*t*u + t*u*x + u*x*y,
   x*y*z*t + y*z*t*u + z*t*u*x + t*u*x*y + u*x*y*z,
   x*y*z*t*u - 1}$

vars := {u, t, z, y, x}$
torder(reverse vars, revgradlex)$

gb := f5(system)$

end;
