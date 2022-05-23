
load_package f5;

on f5interreduce;

load_package dipoly;

procedure errorMessage(system);
  {"Wrong Answer", system};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simple sanity check tests

% Without torder:

gb := f5({x + z + y});
if not (gb = {x + y + z}) then
   errorMessage("Term order problem");

gb := f5({x + z + y}, {z, y, x}, lex);
if not (gb = {z + y + x}) then
   errorMessage("Term order problem");

gb := f5({x + z + y^2}, {z, y, x}, revgradlex);
if not (gb = {y^2 + z + x}) then
   errorMessage("Term order problem");

% With torder:

torder({z, x, y}, lex);
gb := f5({x + y^2 + z});
if not (gb = {z + x + y^2}) then
   errorMessage("Term order problem");

gb := f5({x + y^2 + z});
if not (gb = {z + x + y^2}) then
   errorMessage("Term order problem");

% Weighted order
torder({x,y,z},weighted,{3, 5, 4});
gb := f5({x + y + z});
if not (gb = {y + z + x}) then
   errorMessage("Term order problem");

% Block order
torder({x, y, z}, lexgradlex, 1);
gb := f5({x + y + z^2});
if not (gb = {x + z^2 + y}) then
   errorMessage("Term order problem");

% Graded order
torder({x, y, z}, graded, {1, 1, 2}, lex);
gb := f5({x + y^4 + z^3});
if not (gb = {z^3 + y^4 + x}) then
   errorMessage("Term order problem");

ans := f5({x1}, {x1}, lex);
if not (ans = {x1}) then
  errorMessage({x1});

f5({x1, x2}, {x1, x2}, lex);

f5({x2, x1, x1, x1, x2, x2, x1, x2, x1, x2}, {x1, x2}, lex);

f5({x1 + x2, (x1 + x2)^2, (x1 + x2)^3}, {x1, x2}, lex);

f5({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

f5({x1*x2 + 1, x2*x3 + 1}, {x1, x2, x3}, lex);

f5({x1 + x2 + x3, x1*x2 + x2*x3 + x1*x3, x1*x2*x3 - 1}, {x1, x2, x3}, lex);

ans := f5({10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10}, {x1, x2}, lex)$
if not (length(ans) = 2) then
  errorMessage("noon2 lex");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Several corner cases tests

system := {x*y + x^2 - y^2, x^3*y^2 + y^5, x*y^6 - x^7};
vars := {x, y};

gb := f5(system, vars, lex);
if not (length(gb) = 3) then
  errorMessage(system);
if not (gb = {x^2 + x*y - y^2, x*y^4, y^6}) then
  errorMessage(system);


system := {(x + y), (x + y)^5, (x + y)^100};
vars := {x, y};

gb := f5(system, vars, lex);
if not (gb = {x + y}) then
  errorMessage(system);


system := {5x + 10y, 4x^4 - 8y^4 - 100x^2*y^2};
vars := {x, y};

gb := f5(system, vars, lex);
if not (length(gb) = 2) then
  errorMessage(system);
if not (gb = {x + 2y, y^4}) then
  errorMessage(system);


system := {5x + 10y, 4x^4 - 8y^4 - 100x^2*y^2};
vars := {x, y};

gb := f5(system, vars, revgradlex);
if not (length(gb) = 2) then
  errorMessage(system);
if not (gb = {y^4, x + 2y}) then
  errorMessage(system);


system := {x1^5 - x2^2 - 3, x1*x2^4 + x2^8};
vars := {x1, x2};

gb := f5(system, vars, lex);
if not (length(gb) = 3) then
  errorMessage(system);
if not (gb = {x1^5 - x2^2 - 3, x1*x2^4 + x2^8, x2^24 + x2^6 + 3*x2^4}) then
  errorMessage(system);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Several known systems tests

% Test from `groebner` Reduce package

vars := {q1,q2,q3,q4,q5,q6}$
system := {q1,
          q2^2 + q3^2 + q4^2,
          q4*q3*q2,
          q3^2*q2^2 + q4^2*q2^2 + q4^2*q3^2,
          q6^2 + 1/3*q5^2,
          q6^3 - q5^2*q6,
          2*q2^2*q6 - q3^2*q6 - q4^2*q6 + q3^2*q5 - q4^2*q5,
          2*q2^2*q6^2 - q3^2*q6^2 - q4^2*q6^2 - 2*q3^2*q5*q6
          + 2*q4^2*q5*q6 - 2/3*q2^2*q5^2 + 1/3*q3^2*q5^2
          + 1/3*q4^2*q5^2,
          - q3^2*q2^2*q6 - q4^2*q2^2*q6 + 2*q4^2*q3^2*q6 -
          q3^2*q2^2*q5 + q4^2*q2^2*q5,
          - q3^2*q2^2*q6^2 - q4^2*q2^2*q6^2 + 2*q4^2*q3^2*q6^2
          + 2*q3^2*q2^2*q5*q6 - 2*q4^2*q2^2*q5*q6 + 1/3*q3^2*q2^2
          *q5^2 + 1/3*q4^2*q2^2*q5^2 - 2/3*q4^2*q3^2*q5^2,
          - 3*q3^2*q2^4*q5*q6^2 + 3*q4^2*q2^4*q5*q6^2
          + 3*q3^4*q2^2*q5*q6^2 - 3*q4^4*q2^2*q5*q6^2
          - 3*q4^2*q3^4*q5*q6^2 + 3*q4^4*q3^2*q5*q6^2
          + 1/3*q3^2*q2^4*q5^3 - 1/3*q4^2*q2^4*q5^3
          - 1/3*q3^4*q2^2*q5^3 + 1/3*q4^4*q2^2*q5^3 + 1/3*q4^2
            *q3^4*q5^3 - 1/3*q4^4*q3^2*q5^3}$

gb := f5(system, vars, lex)$

if not (length(gb) = 20) then
  errorMessage("Bad length for model system");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s9_1, regular? NA

vars := {a, b, c, d, e, f, g, h}$
system := {-e*g - 2*d*h,
        9*e + 4*b,
        -4*c*h - 2*e*f - 3*d*g,
        -7*c + 9*a - 8*f,
        -4*d*f - 5*c*g - 6*h - 3*e,
        -5*d - 6*c*f - 7*g + 9*b,
        9*d + 6*a - 5*b,
        9*c - 7*a + 8}$

gb := f5(system, vars, revgradlex)$

if not (length(gb) = 15) then
  errorMessage("Bad length for model system");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% weispfenning94, regular? no

vars := {x, y, z, h}$
system := {y^4 + x*y^2*z + x^2*h^2 - 2*x*y*h^2 + y^2*h^2 + z^2*h^2,
          x*y^4 + y*z^4 - 2*x^2*y*h^2 - 3*h^5,
          -x^3*y^2 + x*y*z^3 + y^4*h + x*y^2*z*h - 2*x*y*h^3}$

% f5(system, vars, lex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% liu, regular? yes

off f5integers;

vars := {x,y,z,t,a,h}$
system := {y*z - y*t - x*h + a*h,
            z*t - z*x - y*h + a*h,
            t*x - y*t - z*h + a*h,
            x*y - z*x - t*h + a*h}$

f5(system, vars, revgradlex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sym33, regular? yes

vars := {h, x, y, z}$
system := {y*z^3 + h^3*x - 2*h^4,
          x^3*z + h^3*y - 2*h^4,
          x*y^3 + h^3*z - 2*h^4}$

f5(system, vars, revgradlex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eco5, regular? NA

system := {
        x1*x2*x5 + x1*x5 + x2*x3*x5 + x3*x4*x5 - 1,
         x1*x3*x5 + x2*x4*x5 + x2*x5 - 2,
          x1*x4*x5 + x3*x5 - 3,
           x4*x5 - 4,
        x1 + x2 + x3 + x4 + 1
};

vars := {x1, x2, x3, x4, x5};

gb := f5(system, vars, revgradlex);
if not (length(gb) = 11) then
  errorMessage("eco 5");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% root-7, regular? no

system := {x1 + x2 + x3 + x4 + x5 + x6 + x7, x1*x2 + x1*x3 + x1*x4 + x1*x5 + x1*x6
+ x1*x7 + x2*x3 + x2*x4 + x2*x5 + x2*x6 + x2*x7 + x3*x4 + x3*x5 + x3*x6 + x3*x7 + x4*x5 + x4*x6 + x4*x7 + x5*x6 + x5*x7
+ x6*x7, x1*x2*x3 + x1*x2*x4 + x1*x2*x5 + x1*x2*x6 + x1*x2*x7 + x1*x3*x4 + x1*x3*x5 + x1*x3*x6 + x1*x3*x7 + x1*x4*x5 + x1*x4*x6 + x1*x4*x7 + x1*x5*x6 + x1*x5*x7 + x1*x6*x7 + x2*x3*x4 + x2*x3*x5 + x2*x3*x6 + x2*x3*x7 + x2*x4*x5 + x2*x4*x6 +
x2*x4*x7 + x2*x5*x6 + x2*x5*x7 + x2*x6*x7 + x3*x4*x5 + x3*x4*x6 + x3*x4*x7 + x3*x5*x6 + x3*x5*x7 + x3*x6*x7 + x4*x5*x6 + x4*x5*x7 + x4*x6*x7 + x5*x6*x7, x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x3*x6 + x1*x2*x3*x7 + x1*x2*x4*x5 + x1*x2*x4*x6 + x1*x2*x4*x7 + x1*x2*x5*x6 + x1*x2*x5*x7 + x1*x2*x6*x7 + x1*x3*x4*x5 + x1*x3*x4*x6 + x1*x3*x4*x7 + x1*x3*x5*x6 + x1*x3*x5*x7 + x1*x3*x6*x7 + x1*x4*x5*x6 + x1*x4*x5*x7 + x1*x4*x6*x7 + x1*x5*x6*x7 + x2*x3*x4*x5 + x2*x3*x4*x6 + x2*x3*x4*x7 + x2*x3*x5*x6 + x2*x3*x5*x7 + x2*x3*x6*x7 + x2*x4*x5*x6 + x2*x4*x5*x7 + x2*x4*x6*x7 + x2*x5*x6*x7 + x3*x4*x5*x6 + x3*x4*x5*x7 + x3*x4*x6*x7 + x3*x5*x6*x7 + x4*x5*x6*x7, x1*x2*x3*x4*x5 + x1*x2*x3*x4*x6 + x1*x2*x3*x4*x7 + x1*x2*x3*x5*x6 + x1*x2*x3*x5*x7 + x1*x2*x3*x6*x7 + x1*x2*x4*x5*x6 + x1*x2*x4*x5*x7 + x1*x2*x4*x6*x7 + x1*x2*x5*x6*x7 + x1*x3*x4*x5*x6 + x1*x3*x4*x5*x7 + x1*x3*x4*x6*x7 + x1*x3*x5*x6*x7 + x1*x4*x5*x6*x7 + x2*x3*x4*x5*x6 + x2*x3*x4*x5*x7 + x2*x3*x4*x6*x7 + x2*x3*x5*x6*x7 + x2*x4*x5*x6*x7 + x3*x4*x5*x6*x7, x1*x2*x3*x4*x5*x6 + x1*x2*x3*x4*x5*x7 + x1*x2*x3*x4*x6*x7 + x1*x2*x3*x5*x6*x7 + x1*x2*x4*x5*x6*x7 + x1*x3*x4*x5*x6*x7 + x2*x3*x4*x5*x6*x7, x1*x2*x3*x4*x5*x6*x7 - 1}$

vars := {x1, x2, x3, x4, x5, x6, x7}$

gb := f5(system, vars, lex);
truegb := {x1 + x2 + x3 + x4 + x5 + x6 + x7, x2^2 + x2*x3 + x2*x4 + x2*x5 + x2*x6 + x2*x7 + x3^2 + x3*x4 + x3*x5 + x3*x6 + x3*x7 + x4^2 + x4*x5 + x4*x6 + x4*x7 + x5^2 + x5*x6 + x5*x7 + x6^2 + x6*x7 + x7^2, x3^3 + x3^2*x4 + x3^2*x5 + x3^2*x6 + x3^2*x7 + x3*x4^2 +
x3*x4*x5 + x3*x4*x6 + x3*x4*x7 + x3*x5^2 + x3*x5*x6 + x3*x5*x7 + x3*x6^2 + x3*x6*x7 + x3*x7^2 + x4^3 + x4^2*x5 + x4^2*x6 + x4^2*x7 + x4*x5^2 + x4*x5*x6 + x4*x5*x7 + x4*x6^2 + x4*x6*x7 + x4*x7^2 + x5^3 + x5^2*x6 + x5^2*x7 + x5*x6^2 + x5*x6*x7 + x5*x7^2 + x6^3 + x6^2*x7 + x6*x7^2 + x7^3, x4^4 + x4^3*x5 + x4^3*x6 + x4^3*x7 + x4^2*x5^2 + x4^2*x5*x6 + x4^2*x5*x7
+ x4^2*x6^2 + x4^2*x6*x7 + x4^2*x7^2 + x4*x5^3 + x4*x5^2*x6
+ x4*x5^2*x7 + x4*x5*x6^2 + x4*x5*x6*x7 + x4*x5*x7^2 + x4*x6^3 + x4*x6^2*x7 + x4*x6*x7^2 + x4*x7^3 + x5^4 + x5^3*x6 + x5^3*x7 + x5^2*x6^2 + x5^2*x6*x7 + x5^2*x7^2 + x5*x6^3 + x5*x6^2*x7 + x5*x6*x7^2 + x5*x7^3 + x6^4 + x6^3*x7 + x6^2*x7^2 +
x6*x7^3 + x7^4, x5^5 + x5^4*x6 + x5^4*x7 + x5^3*x6^2 + x5^3*x6*x7 + x5^3*x7^2 + x5^2*x6^3 + x5^2*x6^2*x7 + x5^2*x6*x7^2
+ x5^2*x7^3 + x5*x6^4 + x5*x6^3*x7 + x5*x6^2*x7^2 + x5*x6*x7^3 + x5*x7^4 + x6^5 + x6^4*x7 + x6^3*x7^2 + x6^2*x7^3 + x6*x7^4 + x7^5, x6^6 + x6^5*x7 + x6^4*x7^2 + x6^3*x7^3 + x6^2*x7^4 + x6*x7^5 + x7^6, x7^7 - 1}$
% if not (gb = truegb) then
  % errorMessage("root7 lex");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% noon-4

system := { 10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 - 11*x1 + 10,
            10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 - 11*x2 + 10,
            10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 - 11*x3 + 10,
            10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 - 11*x4 + 10}$
vars := {x1, x2, x3, x4}$

gb := f5(system, vars, revgradlex);

if not (length(gb) = 28) then
  errorMessage("noon 4 revgradlex");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

off f5interreduce;

% katsura-3

system := {
    x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2,
    2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3,
    2*x0*x2 + x1^2 + 2*x1*x3 - x2,
    x0 + 2*x1 + 2*x2 + 2*x3 - 1
}$
vars := {x0,x1, x2, x3}$

gb := f5(system, vars, revgradlex);

if not (length(gb) = 7) then
  errorMessage("katsura 3 revgradlex");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% katsura-4

system := {
  x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2,
  2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4,
  2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2,
  2*x0*x3 + 2*x1*x2 + 2*x1*x4 - x3,
  x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1}$
vars := {x0, x1, x2, x3, x4}$

gb := f5(system, vars, revgradlex);

if not (length(gb) = 13) then
  errorMessage("katsura 4 revgradlex");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% katsura-5

system := { x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2,
            2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5,
            2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5,
            2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 - x3,
            2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 - x4,
            x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 - 1}$
vars := {x0, x1, x2, x3, x4, x5}$

gb := f5(system, vars, revgradlex);

if not (length(gb) = 22) then
  errorMessage("% katsura 5 revgradlex");

end; % eof
