% Boege, W., Gebauer, R., Kredel, H., Some Examples for Solving Systems of Algebraic Equations by
% Calculating Groebner Bases, J. Symbolic Computation (1986) 1, 83-96,
% https://doi.org/10.1016/S0747-7171(86)80014-1
%
% Hairer, Runge-Kutta 3, 10.11.83

load_package f5;
on f5sugar;

system := {+ B2*C2 + B3*C3 + B4*C4 + B5*C5 - 1/2,
   + B2*C2**2 + B3*C3**2 + B4*C4**2 + B5*C5**2 - 1/3,
   + B3*A32*C2 + B4*A42*C2 + B4*A43*C3 + B5*A52*C2 + B5*A53*C3 + B5*A54*C4 - 1/6,
   + B2*C2**3 + B3*C3**3 + B4*C4**3 + B5*C5**3 - 1/4,
   + B3*C3*A32*C2 + B4*C4*A42*C2 + B4*C4*A43*C3 + B5*C5*A52*C2 + B5*C5*A53*C3 + B5*C5*A54*C4 - 1/8,
   + B3*A32*C2**2 + B4*A42*C2**2 + B4*A43*C3**2 + B5*A52*C2**2 + B5*A53*C3**2 + B5*A54*C4**2 - 1/12,
   + B4*A43*A32*C2 + B5*A53*A32*C2 + B5*A54*A42*C2 + B5*A54*A43*C3 - 1/24,
   + B2*C2**4 + B3*C3**4 + B4*C4**4 + B5*C5**4 - 1/5,
   + B3*C3**2*A32*C2 + B4*C4**2*A42*C2 + B4*C4**2*A43*C3 + B5*C5**2*A52*C2 + B5*C5**2*A53*C3 + B5*C5**2*A54*C4 - 1/10,
   + B3*C2**2*A32*C3 + B4*C2**2*A42*C4 + B4*C3**2*A43*C4 + B5*C2**2*A52*C2 + B5*C3**2*A53*C5 + B5*C4**2*A54*C5 - 1/15,
   + B4*C4*A43*A32*C2 + B5*C5*A53*A32*C2 + B5*C5*A54*A42*C2 + B5*C5*A54*A43*C3 - 1/30,
   + B3*A32**2*C2**2 + B4*A42**2*C2**2 + 2*B4*A42*C2*A43*C3 + B4*A43**2*C3**2 + B5*A52**2*C2**2 + B5*A53**2*C3**2 + B5*A54**2*C4**2 + 2*B5*A52*C2*A53*C3 + 2*B5*A52*C2*A54*C4 + 2*B5*A53*C3*A54*C4 - 1/20,
   + B3*A32*C2**3 + B4*A42*C2**3 + B4*A43*C3**3 + B5*A52*C2**3 + B5*A53*C3**3 + B5*A54*C4**3 - 1/20,
   + B4*A43*C3*A32*C2 + B5*A53*C3*A32*C2 + B5*A54*C4*A42*C2 + B5*A54*C4*A43*C3 - 1/40,
   + B4*A43*A32*C2**2 + B5*A53*A32*C2**2 + B5*A54*A42*C2**2 + B5*A54*A43*C3**2 - 1/60,
   + B5*A54*A43*A32*C2 - 1/120}$

vars := {C2, C3, C4, C5, B2, B3, B4, B5, A32, A42, A43, A52, A53, A54}$
torder(reverse vars, gradlex)$

gb := f5(system)$

end;
