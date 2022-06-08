% Boege, W., Gebauer, R., Kredel, H., Some Examples for Solving Systems of Algebraic Equations by
% Calculating Groebner Bases, J. Symbolic Computation (1986) 1, 83-96,
% https://doi.org/10.1016/S0747-7171(86)80014-1
%
% Hairer, Runge-Kutta 1, 05.11.83

load_package groebner;

system := {+ C2 - A21,
   + C3 - A31 - A32,
   + B1 + B2 + B3 - 1,
   + B2*C2 + B3*C3 - 1/2,
   + B2*C2**2 + B3*C3**2 - 1/3,
   + B3*A32*C2 - 1/6}$

vars := {C2, C3, B3, B2, B1, A21, A32, A31}$
torder(reverse vars, lex)$

gb := groebner(system)$

end;
