

in "bb.red";


gb := buchberger({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

gb := buchberger(
        {x1 + x2 + x3, x1*x2 + x1*x3 + x2*x3, x1*x2*x3 - 1},
        {x1, x2, x3}, lex);

gb := buchberger(
        {10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10},
        {x1, x2}, lex);

end;  % of file


