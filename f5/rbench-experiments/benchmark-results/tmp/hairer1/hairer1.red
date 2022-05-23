% cyclic-5 system in lex

load_package groebner;

system := {
  a - f,
  b - g - h,
  c + d + ee - 1,
  b*c + a*d - 1/2,
  b^2*c + a^2*d - 1/3,
  a*c*g - 1/6
}$

vars := {a,b,c,d,ee,f,g,h}$
torder(vars, revgradlex)$

gb := groebner(system)$

end;
