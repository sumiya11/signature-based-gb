
in "C:\data\projects\mpi\signature-based-gb\f5\f5.red"$

% eco-5 system in degrevlex

system := {
        x1*x2*x5 + x1*x5 + x2*x3*x5 + x3*x4*x5 - 1, x1*x3*x5 + x2*x4*x5 + x2*x5 - 2, x1*x4*x5 + x3*x5 - 3, x4*x5 - 4,
        x1 + x2 + x3 + x4 + 1
};

vars := {x1, x2, x3, x4, x5};

share system, vars;

lisp;
st := time();

algebraic;
gb := f5(system, vars, 'revgradlex)$

lisp;
prin2t({"TIME", time() - st});

end;