% sparse-5 system in degrevlex

load_package f5;

system := {
x1^2*x2^2*x3^2*x4^2*x5^2 + 3*x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + x4^2 + x5^2 + 5,
x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + 3*x2^2 + x3^2 + x4^2 + x5^2 + 5,
x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + 3*x3^2 + x4^2 + x5^2 + 5,
x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + 3*x4^2 + x5^2 + 5,
x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x1*x2*x3*x4*x5 + x2^2 + x3^2 + x4^2 + 3*x5^2 + 5
};

vars := {x1, x2, x3, x4, x5};

basis := f5(system, vars, revgradlex);

end;
