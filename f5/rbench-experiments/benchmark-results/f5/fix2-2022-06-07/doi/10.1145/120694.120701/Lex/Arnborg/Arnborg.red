% Giovini, A., Mora, T., Niesi, G., Robbiano, L. and Traverso, C., "One sugar cube, please" or
% Selection strategies in the Buchberger algorithm, J. of the ACM 1991, pp. 49-54,
% https://doi.org/10.1145/120694.120701
%
% 3.3 Arnborg-Lazard system

load_package f5;
on f5sugar;

system := {x**2*y*z + x*y**2*z + x*y*z**2 + x*y*z + x*y + x*z + y*z,
   x**2*y**2*z + x**2*y*z + x*y**2*z**2 + x*y*z + x + y*z + z,
   x**2*y**2*z**2 + x**2*y**2*z + x*y**2*z + x*y*z + x*z + z + 1}$

vars := {z, y, x}$
torder(reverse vars, lex)$

gb := f5(system)$

end;
