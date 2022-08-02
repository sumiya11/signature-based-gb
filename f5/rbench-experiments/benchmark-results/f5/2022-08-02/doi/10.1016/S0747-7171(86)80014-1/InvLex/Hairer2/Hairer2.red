% Boege, W., Gebauer, R., Kredel, H., Some Examples for Solving Systems of Algebraic Equations by
% Calculating Groebner Bases, J. Symbolic Computation (1986) 1, 83-96,
% https://doi.org/10.1016/S0747-7171(86)80014-1
%
% Hairer, Runge-Kutta 2, 05.11.83

load_package f5;
1;

system := {+ B1 + B2 + B3 + B4 - 1,
   + B2*C2 + B3*C3 + B4*C4 - 1/2,
   + B2*C2**2 + B3*C3**2 + B4*C4**2 - 1/3,
   + B3*A32*C2 + B4*A42*C2 + B4*A43*C3 - 1/6,
   + B2*C2**3 + B3*C3**3 + B4*C4**3 - 1/4,
   + B3*C3*A32*C2 + B4*C4*A42*C2 + B4*C4*A43*C3 - 1/8,
   + B3*A32*C2**2 + B4*A42*C2**2 + B4*A43*C3**2 - 1/12,
   + B4*A43*A32*C2 - 1/24,
   + C2 - A21,
   + C3 - A31 - A32,
   + C4 - A41 - A42 - A43}$

vars := {C2, C3, C4, B4, B3, B2, B1, A21, A31, A32, A41, A42, A43}$
torder(reverse vars, lex)$

gb := f5(system)$

end;
