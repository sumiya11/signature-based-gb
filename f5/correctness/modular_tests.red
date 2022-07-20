
% Modular computations Tests

load_package f5;

procedure modularMessage(system);
  {"Modular Error", system};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rational reconstruction tests

lisp;

procedure modularMessageSym(system);
  {"Modular Error", system};

p := 1031;
truerecs := {0 . 1, 1 . 1, 2 . 1, 3 . 1, 4 . 1, 5 . 1, 6 . 1, 7 . 1, 8 . 1, 9 . 1,
            10 . 1, 11 . 1, 12 . 1, 13 . 1, 14 . 1, 15 . 1, 16 . 1, 17 . 1, 18 . 1,
            19 . 1, 20 . 1, 21 . 1, (-19) . 46, (-19) . 44, 1 . 43, (-6) . 41, (-17) . 39, (-5) . 38, 5 . 37, (-16) . 35,
            (-11) . 34, (-8) . 33, (-7) . 32, (-8) . 31, (-11) . 30, (-16) . 29, 13 . 29, 5 . 28, (-5) . 27,
            (-17) . 26, 9 . 26, (-6) . 25, 19 . 25, 1 . 24, (-19) . 23, 4 . 23, (-19) . 22, 3 . 22, 2 . 43, (-2) . 21};

for i := 0:48 do <<
  truerec := pop(truerecs);
  rec := mod_reconstruction(i, p);
  if not (rec = truerec) then
    prin2t(modularMessage("fail of rational reconstruction"))
>>;

algebraic;

on f5modular;

gb := f5({x, 0, y, 5, 0});
if not (gb = {1}) then
   modularMessage("special case error");

gb := f5({x, 0, y, 0, 0});
if not (gb = {x, y}) then
   modularMessage("special case error");

gb := f5({x});
if not (gb = {x}) then
   modularMessage("basic error");

on f5fractionfree;

gb := f5({x + 3, y - 5});
if not (gb = {x + 3, y - 5}) then
   modularMessage("fractionfree,basic error");

gb := f5({4*x^2 - 8, x*y^4 + 3*y^3});
if not (gb = {x^2 - 2, 3*x*y^3 + 2*y^4, 2*y^5 - 9*y^3}) then
   modularMessage("fractionfree,basic error");

gb := f5({x^2*y - 12345*x - 1234567890, x*y^2 + 67890*y + 1234567890});
if not (gb = {376335277782482*x - 24458161643*y^2  + 376443191958587*y + 1962400599426105,
              13717421*y^3  - 60523935*y^2  - 2031892985625*y - 16935087500211690}) then
   modularMessage("fractionfree,big arithmetic error");

gb := f5({234324355345345435*x - 3242343242342342342341});
if not (gb = {234324355345345435*x - 3242343242342342342341}) then
   modularMessage("fractionfree,big arithmetic error");

off f5fractionfree;

gb := f5({234324355345345435*x - 3242343242342342342341});
if not (gb = {x - 3242343242342342342341/234324355345345435}) then
   modularMessage("big arithmetic error");

end; % of file
