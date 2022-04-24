
load_package f5;

% ojika-4 system in degrevlex

system := {
        x1+x3*x1^3+x1*x3*x2^2-x1*x3,
        10*x2-2*x2*x3*x1^2-x3*x2^3-x2*x3,
        -6*x3^2*x1^4-3*x1^2*x2^2*x3^2-x3^2*x1^2+28*x3*x1^2 - 3*x3^2*x2^4+2*x3^2*x2^2+7*x3*x2^2+x3^2-11*x3+10
};

vars := {x1, x2, x3};

share system, vars;

lisp;
st := time();

algebraic;
gb := f5(system, vars, 'revgradlex)$

lisp;
prin2t({"TIME", time() - st});

end;
