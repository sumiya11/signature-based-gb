
% Tests that were previously discovered to fail

load_package f5;

on f5fullreduce;

procedure regressionMessage(system);
  {"Regression", system};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s9_1, regular? NA

% car on nil error

vars := {a, b, c, d, e, f, g, h}$
system := {-e*g - 2*d*h,
        9*e + 4*b,
        -4*c*h - 2*e*f - 3*d*g,
        -7*c + 9*a - 8*f,
        -4*d*f - 5*c*g - 6*h - 3*e,
        -5*d - 6*c*f - 7*g + 9*b,
        9*d + 6*a - 5*b,
        9*c - 7*a + 8}$

% f5(system, vars, lex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bug in pair selection

vars := {x1,x2,x3}$
system := {x1*x2 + x2, x1^2 - x2^2, x1^2 + x2 + x3}$

gb := f5(system, vars, lex);
truegb := {x1^2 + x2 + x3, x1*x2 + x2, x1*x3 + x3, x2^2 + x2 + x3, x2*x3 - x3, x3^2 + 2*x3};

if not (gb = truegb) then
  regressionMessage("Bug in pair selection");

end; % of file
