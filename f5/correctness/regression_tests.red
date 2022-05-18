
% Tests that were previously discovered to fail

load_package f5;

on f5interreduce;

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

gb := f5(system, vars, revgradlex);

if not (length(gb) = 15) then
  regressionMessage("car on nil");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bug in pair selection

vars := {x1,x2,x3}$
system := {x1*x2 + x2, x1^2 - x2^2, x1^2 + x2 + x3}$

gb := f5(system, vars, lex);
truegb := {x1^2 + x2 + x3,
          x1*x2 + x2,
          x1*x3 + x3,
          x2^2 + x2 + x3,
           x2*x3 - x3,
            x3^2 + 2*x3};

if not (length(gb) = length(truegb)) then
  regressionMessage("Bug in pair selection");
% if not (gb = truegb) then
%   regressionMessage("Bug in pair selection");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unit ideal by error

vars := {x, y, z}$
system := {x*y^2 - 10x^2 - 11y^2, x^3 - x + y^3, x - z^4}$

if not (length(f5(system, vars, revgradlex)) = 4) then
   regressionMessage("unit ideal error");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pair selection bug x2

vars := {x1,x2,x3}$
system := {10*x1^2*x3 + 10*x2^2*x3 - 11*x3 + 10,
          10*x1*x2^2 + 10*x1*x3^2 - 11*x1 + 10,
          10*x1^2*x2 + 10*x2*x3^2 - 11*x2 + 10}$

if not (length(f5(system, vars, revgradlex)) = 11) then
   regressionMessage("unit ideal error");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end; % of file
