% katsura-6 system in degrevlex

load_package f5;

off f5fractionfree;

setmod 65537;
on modular;

system := {x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 - 1,
    2*x1*x6 + 2*x2*x5 + 2*x2*x7 + 2*x3*x4 - x6,
    2*x1*x5 + 2*x2*x4 + 2*x2*x6 + x3^2 + 2*x3*x7 - x5,
    2*x1*x4 + 2*x2*x3 + 2*x2*x5 + 2*x3*x6 + 2*x4*x7 - x4,
    2*x1*x3 + x2^2 + 2*x2*x4 + 2*x3*x5 - x3 + 2*x4*x6 + 2*x5*x7,
    2*x1*x2 + 2*x2*x3 - x2 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7,
    x1^2 - x1 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2};

vars := {x1, x2, x3, x4, x5, x6, x7}$

share system, vars;

lisp;
st := time();

algebraic;
gb := f5(system, vars, 'revgradlex)$

lisp;
prin2t({"TIME", time() - st});

end;
