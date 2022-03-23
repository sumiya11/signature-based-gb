
% Simple correctness tests of f5

in "f5.red";

% how to load this file hm

f5({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

f5({x1*x2 + 1, x2*x3 + 1}, {x1, x2, x3}, lex);

f5({x1 + x2 + x3, x1*x2 + x2*x3 + x1*x3, x1*x2*x3 - 1}, {x1, x2, x3}, lex);

f5({10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10}, {x1, x2}, revgradlex);

end;  % of file


