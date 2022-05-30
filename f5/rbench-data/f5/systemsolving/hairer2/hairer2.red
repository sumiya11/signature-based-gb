% hairer-2 system in revgradlex
% characteristic 0
% 0 dim
%
% PoSSo test suite
% https://www-sop.inria.fr/saga/POL/BASE/3.posso/ha-mo1.dir/index.html

load_package f5;

system := {
  d + ee + f + g - 1,
  c*d + b*ee + a*f - 1/2,
  c^2*d + b^2*ee + a^2*f - 1/3,
  a*ee*ii + a*d*l + b*d*m - 1/6,
  c^3*d + b^3*ee + a^3*f - 1/4,
  a*b*ee*ii + a*c*d*l + b*c*d*m - 1/8,
  a^2*ee*ii + a^2*d*l + b^2*d*m - 1/2,
  a*d*ii*m - 1/24,
  a - h,
  b - ii - j,
  c - k - l - m
}$

vars := {a,b,c,d,ee,f,g,h,ii,j,k,l,m}$
torder(vars, revgradlex)$

gb := f5(system)$

end;
