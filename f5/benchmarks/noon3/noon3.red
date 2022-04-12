% noon-3 system in degrevlex

load_package f5;

system := {
    10*x1*x2^2 + 10*x1*x3^2 - 11*x1 + 10,
    10*x1^2*x2 + 10*x2*x3^2 - 11*x2 + 10,
    10*x1^2*x3 + 10*x2^2*x3 - 11*x3 + 10
};

vars := {x1, x2, x3};

ord := gradlex;

basis := f5(system, vars, revgradlex);

end; % end of file
