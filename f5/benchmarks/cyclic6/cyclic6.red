% cyclic-5 system in lex

in "C:\data\projects\mpi\signature-based-gb\f5\f5.red"$

system := {
x1 + x2 + x3 + x4 + x5 + x6, x1*x2 + x1*x3 + x1*x4 + x1*x5 + x1*x6 + x2*x3 + x2*x4 + x2*x5 + x2*x6 + x3*x4 + x3*x5 + x3*x6 + x4*x5 + x4*x6 + x5*x6, x1*x2*x3 + x1*x2*x4 + x1*x2*x5 + x1*x2*x6 + x1*x3*x4 + x1*x3*x5 + x1*x3*x6 + x1*x4*x5 + x1*x4*x6 + x1*x5*x6
+ x2*x3*x4 + x2*x3*x5 + x2*x3*x6 + x2*x4*x5 + x2*x4*x6 + x2*x5*x6 + x3*x4*x5 + x3*x4*x6 + x3*x5*x6 + x4*x5*x6, x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x3*x6 + x1*x2*x4*x5 + x1*x2*x4*x6 + x1*x2*x5*x6 + x1*x3*x4*x5 + x1*x3*x4*x6 + x1*x3*x5*x6 + x1*x4*x5*x6 + x2*x3*x4*x5 + x2*x3*x4*x6 + x2*x3*x5*x6 + x2*x4*x5*x6 + x3*x4*x5*x6, x1*x2*x3*x4*x5 + x1*x2*x3*x4*x6 + x1*x2*x3*x5*x6 + x1*x2*x4*x5*x6 + x1*x3*x4*x5*x6 + x2*x3*x4*x5*x6, x1*x2*x3*x4*x5*x6 + 1
}$

vars := {x1, x2, x3, x4, x5, x6}$

share system, vars;

lisp;
st := time();

algebraic;
gb := f5(system, vars, 'lex)$

lisp;
prin2t({"TIME", time() - st});

end;
