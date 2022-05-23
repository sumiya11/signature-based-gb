% cassou system in revgradlex
% 0 dim

load_package f5;

system := {
  a - f,
  b - g - h,
  c + d + ee - 1,
  b*c + a*d - 1/2,
  b^2*c + a^2*d - 1/3,
  a*c*g - 1/6
}$

vars := {a,b,c,d,ee,f,g,h}$

gb := f5(system, vars, 'revgradlex)$

end;
