
load_package f5;

procedure errorMessage(system);
  {"Wrong Answer", system};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Algebraic mode interface

gb := f5({0});
if not (gb = {0}) then
   errorMessage("special cases problem");

gb := f5({1});
if not (gb = {1}) then
   errorMessage("special cases problem");

gb := f5({x, 3});
if not (gb = {1}) then
   errorMessage("special cases problem");

gb := f5({0, x, 0});
if (gb = {1}) then
   errorMessage("special cases problem");

off f5parametric;
on f5fractionfree;

gb := f5({x});
if not (gb = {x}) then
   errorMessage("fractionfree");

gb := f5({x^2 + x, y});
if not (gb = {x^2 + x, y}) then
   errorMessage("fractionfree");

gb := f5({x^2 + x, y});
if not (gb = {x^2 + x, y}) then
  errorMessage("fractionfree");

gb := f5({(a + 1)*x + 1, a*y}, {x, y}, lex);
if not (gb = {(a + 1)*x + 1, y}) then
  errorMessage("fractionfree,params");

gb := f5({a*x^2 + a^2*x, a*y}, {x, y}, lex);
if not (gb = {x^2 + a*x, y}) then
  errorMessage("fractionfree,params");

gb := f5({a*x^2 + x/(b^3 - 8), y/a}, {x, y}, lex);
if not (gb = {a*b^3*x^2 - 8*a*x^2 + x, y}) then
  errorMessage("fractionfree,params");

gb := f5({x/a});
if not (gb = {x}) then
  errorMessage("fractionfree,params");

gb := f5({a*x - 1, -a*y - 1, -a*z + 1, a*w + 1}, {x,y,z,w}, lex);
if not (gb = {a*x - 1, a*y + 1, a*z - 1, a*w + 1}) then
  errorMessage("fractionfree,params");

off f5fractionfree;

gb := f5({12*x - 11, 13*y + 14, -5*z - 5});
if not (gb = {x - 11/12, y + 14/13, z + 1}) then
  errorMessage("");

gb := f5({x/a});
if not (gb = {x}) then
  errorMessage("params");

gb := f5({a*x - 1, -a*y - 1, -a*z + 1, a*w + 1}, {x,y,z,w}, lex);
if not (gb = {x - 1/a, y + 1/a, z - 1/a, w + 1/a}) then
  errorMessage("params");

gb := f5({a*x + 1/(b + 1)*y, x*y - c}, {x, y}, lex);
if not (gb = {x + 1/(a*b + a)*y, y^2 + a*b*c + a*c}) then
  errorMessage("params");

on f5parametric;
on f5fractionfree;

gb := f5({2*x + 3, -2*y + 3, 2*z - 3, -2*w - 3}, {x,y,z,w}, lex);
if not (gb = {2*x + 3, 2*y - 3, 2*z - 3, 2*w + 3}) then
  errorMessage("fractionfree");

gb := f5({-5*x + 10});
if not (gb = {x - 2}) then
  errorMessage("fractionfree");

gb := f5({(-5*x + 10)/3 + y, y - 8/11});
if not (gb = {5x - 3*y - 10, 11*y - 8}) then
  errorMessage("fractionfree");

gb := f5({a*x});
if not (gb = {a*x}) then
  errorMessage("fractionfree,params");

gb := f5({a*x^4 + a*x^3 - a^2*b*x^2 - a*b^10*x + a^13}, {x}, lex);
if not (gb = {a*x^4 + a*x^3 - a^2*b*x^2 - a*b^10*x + a^13}) then
  errorMessage("fractionfree,params");

gb := f5({(a*x - 4*a)/a}, {x}, lex);
if not (gb = {x - 4}) then
  errorMessage("fractionfree,params");

gb := f5({2*a^2*x^2 - 4*a*x, -a*y^2 - a*y + a^11}, {x,y}, lex);
if not (gb = {a^2*x^2 - 2*a*x, a*y^2 + a*y - a^11}) then
  errorMessage("fractionfree,params");

gb := f5({a*x*y - a, a*x + a^2}, {x,y}, lex);
if not (gb = {a*x + a^2, a^3*y + a^2}) then
  errorMessage("fractionfree,params");

on f5parametric;
off f5fractionfree;

gb := f5({2*a*x*y - a, 3*a*x + 6*a^2}, {x,y}, lex);
if not (gb = {a*x + 2*a^2, a^3*y + a^2/4}) then
   errorMessage("fractionfree,params");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% f5 over Zp

% off f5parametric;
% off f5fractionfree;

% setmod 7;
% on modular;

% gb := f5({x + y}, {x, y}, lex);
% if not (gb = {x + y}) then
%    errorMessage("on modular");

% off modular;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Symbolic mode interface

f1 := x*y - 1; f2 := y^2 - x;
g1 := x^2 - y; g2 := x*y - 1; g3 := y^2 - x;

share x,y,z,w,f1,f2,g1,g2,g3;

lisp;

procedure errorMessageSym(system);
  {"Wrong Answer", system};

torder({{'list, w,z,y,x}, 'lex});

gb := f5_groebnerq({simp x, simp 0, simp y, simp z, simp 0, simp w});
if not (gb = {simp w, simp z, simp y, simp x}) then
  errorMessageSym("symmode,sq");

gb := f5_groebnerf({numr simp x, numr simp y, numr simp z, numr simp w});
if not (gb = {numr simp w, numr simp z, numr simp y, numr simp x}) then
  errorMessageSym("symmode,sf");

torder({{'list, x,y,z,w}, 'gradlex});

gb := f5_groebnerq({simp f1, simp f2});
if not (gb = {simp g1, simp g2, simp g3}) then
  errorMessageSym("symmode,sq");

gb := f5_groebnerp({'w, 'x, 'z, 'y});
if not (gb = {'x,'y,'z,'w}) then
  errorMessageSym("symmode,lp");

torder({{'list, x,y,z,w}, 'lex});

gb := f5_groebnerp({'(plus (expt x 2) (times 2 z) (expt w 3) (expt y 4))});
if not (gb = {'(plus (expt x 2) (expt y 4) (times 2 z) (expt w 3))}) then
  errorMessageSym("symmode,lp");

algebraic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic interface

% share x,y,z,w,f1,f2,g1,g2,g3;

% lisp;

% torder({{'list, w,z,y,x}, 'lex});

% gb of {x, 0, y, z, 0, w} in Standard Quotient
% f5_groebnerq({simp x, simp 0, simp y, simp z, simp 0, simp w});

% gb of {x, y, z, w} in Standard Form
% f5_groebnerf({numr simp x, numr simp y, numr simp z, numr simp w});

% gb of {w, x, z, y} as in Lisp Prefix
% f5_groebnerp({'w, 'x, 'z, 'y});
% algebraic;


end; % eof
