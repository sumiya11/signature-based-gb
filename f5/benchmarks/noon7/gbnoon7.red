% noon-6 system in degrevlex

load_package groebner;

system := {10*x1*x2^2 + 10*x1*x3^2 +
10*x1*x4^2 + 10*x1*x5^2 + 10*x1*x6^2 + 10*x1*x7^2 - 11*x1 + 10, 10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 + 10*x2*x5^2 + 10*x2*x6^2 + 10*x2*x7^2 - 11*x2
+ 10, 10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 + 10*x3*x5^2 + 10*x3*x6^2 + 10*x3*x7^2 - 11*x3 + 10, 10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 + 10*x4*x5^2 +
10*x4*x6^2 + 10*x4*x7^2 - 11*x4 + 10, 10*x1^2*x5 + 10*x2^2*x5 + 10*x3^2*x5 + 10*x4^2*x5 + 10*x5*x6^2 + 10*x5*x7^2 - 11*x5 + 10, 10*x1^2*x6 + 10*x2^2*x6 + 10*x3^2*x6 + 10*x4^2*x6 + 10*x5^2*x6 + 10*x6*x7^2 - 11*x6 + 10, 10*x1^2*x7 + 10*x2^2*x7 + 10*x3^2*x7 + 10*x4^2*x7 + 10*x5^2*x7 + 10*x6^2*x7 -
11*x7 + 10};

vars := {x1, x2, x3, x4, x5, x6, x7};
torder(vars, revgradlex)$

share system;

lisp;
st := time();

algebraic;
gb := groebner(system)$

lisp;
prin2t({"TIME", time() - st});

end;
