
% cyclic-5 system in lex

system := {
  x1 + x2 + x3 + x4 + x5,
  x1*x2 + x1*x3 + x1*x4 + x1*x5 + x2*x3 + x2*x4 + x2*x5 + x3*x4 + x3*x5 + x4*x5,
  x1*x2*x3 + x1*x2*x4 + x1*x2*x5 + x1*x3*x4 + x1*x3*x5 + x1*x4*x5 + x2*x3*x4 + x2*x3*x5 + x2*x4*x5 + x3*x4*x5,
  x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5,
  x1*x2*x3*x4*x5 - 1
};

vars := {x1, x2, x3, x4, x5};
ordering := lex;
torder(vars, ordering);

