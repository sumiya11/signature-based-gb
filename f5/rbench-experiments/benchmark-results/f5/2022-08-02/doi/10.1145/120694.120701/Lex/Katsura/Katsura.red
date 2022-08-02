% Giovini, A., Mora, T., Niesi, G., Robbiano, L. and Traverso, C., "One sugar cube, please" or
% Selection strategies in the Buchberger algorithm, J. of the ACM 1991, pp. 49-54,
% https://doi.org/10.1145/120694.120701
%
% 3.5 The Katsura-4 example

load_package f5;
1;

system := {2*x**2 + 2*y**2 + 2*z**2 + 2*t**2 + u**2 - u, x*y + 2*y*z + 2*z*t + 2*t*u - t,
   2*x*z + 2*y*t + t**2 + 2*z*u - z, 2*x*t + 2*z*t + 2*y*u - y,
   2*x + 2*y + 2*z + 2*t + u - 1}$

vars := {u, t, z, y, x}$
torder(reverse vars, lex)$

gb := f5(system)$

end;
