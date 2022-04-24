
% Modular computations Tests

load_package f5;

on f5fullreduce;

procedure modularMessage(system);
  {"Modular Error", system};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rational reconstruction tests

lisp;

in "C:\data\projects\mpi\signature-based-gb\f5\f5.red"$

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
    prin2t(modularMessage("faile of rational reconstruction"))
>>;


end; % of file
