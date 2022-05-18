% Correctness tests of f5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sanity checks

torder({x1}, lex);

f5({x1});

f5({x1}, {x1}, lex);

torder({x1, x2}, revgradlex);

f5({x1, x2^2});

f5({x1 + 1, x1});

f5({x2, x1, x1, x1, x2, x2, x1, x2, x1, x2}, {x1, x2}, lex);

f5({x1 + x2, (x1 + x2)^2, (x1 + x2)^3}, {x1, x2}, lex);

f5({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

f5({x1*x2 + 1, x2*x3 + 1}, {x1, x2, x3}, lex);

f5({x1 + x2 + x3, x1*x2 + x2*x3 + x1*x3, x1*x2*x3 - 1}, {x1, x2, x3}, lex);

f5({10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10}, {x1, x2}, lex);

noon3 := {10*x1*x2^2 + 10*x1*x3^2 - 11*x1 + 10,
          10*x1^2*x2 + 10*x2*x3^2 - 11*x2 + 10,
          10*x1^2*x3 + 10*x2^2*x3 - 11*x3 + 10}$
f5(noon3, {x1, x2, x3}, revgradlex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for different combinations of switches

off f5integers;
f5({5x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

on f5statistics;
on f5integers;
f5({5x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

f5(noon3, {x1, x2, x3}, revgradlex);

f5({x1 + x2, x1*x2 + 100}, {x1, x2}, lex);
f5({4x1 + x2, 1234x1*x2 + 1e5}, {x1, x2}, lex);

off f5statistics;

% the number 4194319 is special because it is the default prime
% number used in modular reduction
f5({x1 + x2, x1*x2 + 4194319}, {x1, x2}, lex);
f5({x1 + 4194329*x2, x1*x2 + 4194319}, {x1, x2}, lex);

off f5integers;
f5({x1 + x2, x1*x2 + 1e10}, {x1, x2}, lex);
f5({x1 + x2, 2347624x1*x2 + 1e100}, {x1, x2}, lex);

on f5integers;
f5({x1 + x2, x1*x2 + 1e10}, {x1, x2}, lex);
f5({x1 + x2, 2347624x1*x2 + 1e100}, {x1, x2}, lex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test from groebner Reduce package

vars := {q1,q2,q3,q4,q5,q6}$
system := {q1,
          q2**2 + q3**2 + q4**2,
          q4*q3*q2,
          q3**2*q2**2 + q4**2*q2**2 + q4**2*q3**2,
          q6**2 + 1/3*q5**2,
          q6**3 - q5**2*q6,
          2*q2**2*q6 - q3**2*q6 - q4**2*q6 + q3**2*q5 - q4**2*q5,
          2*q2**2*q6**2 - q3**2*q6**2 - q4**2*q6**2 - 2*q3**2*q5*q6
          + 2*q4**2*q5*q6 - 2/3*q2**2*q5**2 + 1/3*q3**2*q5**2
          + 1/3*q4**2*q5**2,
          - q3**2*q2**2*q6 - q4**2*q2**2*q6 + 2*q4**2*q3**2*q6 -
          q3**2*q2**2*q5 + q4**2*q2**2*q5,
          - q3**2*q2**2*q6**2 - q4**2*q2**2*q6**2 + 2*q4**2*q3**2*q6**2
          + 2*q3**2*q2**2*q5*q6 - 2*q4**2*q2**2*q5*q6 + 1/3*q3**2*q2**2
          *q5**2 + 1/3*q4**2*q2**2*q5**2 - 2/3*q4**2*q3**2*q5**2,
          - 3*q3**2*q2**4*q5*q6**2 + 3*q4**2*q2**4*q5*q6**2
          + 3*q3**4*q2**2*q5*q6**2 - 3*q4**4*q2**2*q5*q6**2
          - 3*q4**2*q3**4*q5*q6**2 + 3*q4**4*q3**2*q5*q6**2
          + 1/3*q3**2*q2**4*q5**3 - 1/3*q4**2*q2**4*q5**3
          - 1/3*q3**4*q2**2*q5**3 + 1/3*q4**4*q2**2*q5**3 + 1/3*q4**2
            *q3**4*q5**3 - 1/3*q4**4*q3**2*q5**3}$

f5(system, vars, lex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests from Groebner.jl Julia package
%   https://github.com/sumiya11/Groebner.jl

% Cyclic-5
system := {
  x1 + x2 + x3 + x4 + x5,
  x1*x2 + x1*x3 + x1*x4 + x1*x5 + x2*x3 + x2*x4 + x2*x5 + x3*x4 + x3*x5 + x4*x5,
  x1*x2*x3 + x1*x2*x4 + x1*x2*x5 + x1*x3*x4 + x1*x3*x5 + x1*x4*x5 + x2*x3*x4 + x2*x3*x5 + x2*x4*x5 + x3*x4*x5,
  x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5,
  x1*x2*x3*x4*x5 - 1
}$
vars := {x1, x2, x3, x4, x5}$
f5(system, vars, lex);

% Sparse-5
system := {
  x1^2*x2^2*x3^2*x4^2*x5^2 + 3*x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + 3*x2^2 + x3^2 + x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + 3*x3^2 + x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + 3*x4^2 + x5^2 + 5,
  x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + x4^2 + 3*x5^2 + 5
}
vars := {x1, x2, x3, x4, x5}$
f5(system, vars, revgradlex);

on f5integers;

% Katsura-5
system := {
  x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2,
  2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5,
  2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5,
  2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 - x3,
  2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 - x4,
  x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 - 1
}$
vars := {x0,x1,x2,x3,x4,x5}$
f5(system, vars, revgradlex);

end;  % of file
