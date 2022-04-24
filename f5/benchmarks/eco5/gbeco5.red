% eco-5 system in degrevlex

load_package groebner;

system := {
        x1*x2*x5 + x1*x5 + x2*x3*x5 + x3*x4*x5 - 1, x1*x3*x5 + x2*x4*x5 + x2*x5 - 2, x1*x4*x5 + x3*x5 - 3, x4*x5 - 4,
        x1 + x2 + x3 + x4 + 1
};

vars := {x1, x2, x3, x4, x5};
torder(vars, revgradlex)$

share system;

lisp;
st := time();

algebraic;
gb := groebner(system)$

lisp;
prin2t({"TIME", time() - st});

end;
