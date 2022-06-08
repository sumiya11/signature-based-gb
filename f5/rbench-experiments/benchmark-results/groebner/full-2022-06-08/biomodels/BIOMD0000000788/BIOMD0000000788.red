load_package groebner$
torder({}, revgradlex)$

k1 := 3;
k2 := 1/5;
k3 := 3/4;
k4 := 1/10;
k5 := 1/10;
k6 := 3/100;
k7 := 10;
k8 := 1/100;
k9 := 1;
k10 := 1/100;
k11 := 1;
k12 := 1/10;
k13 := 10;
k14 := 1/10;
k15 := 1/20;
k16 := 1/20;
k17 := 1/10;
k18 := 1/100;
k19 := 1/100;
k20 := 1;
k21 := 10;
k22 := 1;
k23 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k22*k7*x1*x2 + (-1)*k22*k9*x1*x3 + 1*k22*k8*x4 + 1*k22*k10*x5 + (-1)*k22*k4*x1 + (-1)*k22*k5*x1 + 1*k6*x7/k1 + 1*k2*x8*k1)/k22,
  diff(x2, t) = ((-1)*k22*k7*x1*x2 + (-1)*k22*k21*x5*x2 + 1*k22*k8*x4 + 1*k22*k19*x6 + 1*k22*k11 + (-1)*k22*k12*x2)/k22,
  diff(x3, t) = ((-1)*k22*k9*x1*x3 + (-1)*k22*k20*x4*x3 + 1*k22*k10*x5 + 1*k22*k18*x6 + 1*k22*k13 + (-1)*k22*k14*x3)/k22,
  diff(x4, t) = (1*k22*k7*x1*x2 + (-1)*k22*k20*x4*x3 + 1*k22*k19*x6 + (-1)*k22*(k8 + k15)*x4)/k22,
  diff(x5, t) = (1*k22*k9*x1*x3 + (-1)*k22*k21*x5*x2 + 1*k22*k18*x6 + (-1)*k22*(k10 + k16)*x5)/k22,
  diff(x6, t) = (1*k22*k20*x4*x3 + 1*k22*k21*x5*x2 + (-1)*k22*(k18 + k19 + k17)*x6)/k22,
  diff(x7, t) = (1*k5*x1*k1 + (-1)*k23*k6*x7)/k23,
  diff(x8, t) = (-1)*k23*k2*x8/k23 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file