% Giovini, A., Mora, T., Niesi, G., Robbiano, L. and Traverso, C., "One sugar cube, please" or
% Selection strategies in the Buchberger algorithm, J. of the ACM 1991, pp. 49-54,
% https://doi.org/10.1145/120694.120701
%
% 3.1 Gerdt examples, example 2

load_package f5;

system := {35*y**4 - 30*x*y**2 - 210*y**2*z + 3*x**2 + 30*x*z - 105*z**2 + 140*y*t - 21*u,
   5*x*y**3 - 140*y**3*z - 3*x**2*y + 45*x*y*z - 420*y*z**2 + 210*y**2*t - 25*x*t + 70*z*t + 126*y*u}$

vars := {u, t, z, y, x}$
torder(reverse vars, lex)$

gb := f5(system)$

end;
