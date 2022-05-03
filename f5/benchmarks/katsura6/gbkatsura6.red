% cyclic-5 system in lex

load_package groebner;

setmod 65537;
on modular;

system := {
  x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2,
   2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6,
    2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5 + 2*x4*x6,
     2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 - x3,
     2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 - x4,
      2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 - x5,
      x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 - 1
};

vars := {x0, x1, x2, x3, x4, x5, x6}$

torder(vars, revgradlex)$

share system;

lisp;
st := time();

algebraic;
gb := groebner(system)$

lisp;
prin2t({"TIME", time() - st});

end;
