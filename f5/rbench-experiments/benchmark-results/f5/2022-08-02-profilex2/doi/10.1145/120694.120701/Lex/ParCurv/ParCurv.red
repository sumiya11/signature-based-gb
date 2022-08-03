% Giovini, A., Mora, T., Niesi, G., Robbiano, L. and Traverso, C., "One sugar cube, please" or
% Selection strategies in the Buchberger algorithm, J. of the ACM 1991, pp. 49-54,
% https://doi.org/10.1145/120694.120701
%
% 3.4 A parametric curve

load_package f5;

system := {x**31 - x**6 - x - y,
   x**8 - z, x**10 - t}$

vars := {t, z, y, x}$
torder(reverse vars, lex)$

gb := f5(system)$
in "~/signature-based-gb/f5/sed.red";

end;
