load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 2;
k3 := 1;
k4 := 1;
k5 := 1/5;
k6 := 1;
k7 := 50;
k8 := 50;
k9 := 500;
k10 := 10;
k11 := 50;
k12 := 1;
k13 := 100;
k14 := 1/100;
k15 := 50;
k16 := 1/2;
k17 := 5;
k18 := 0;
k19 := 1;
k20 := 1;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = ((-1)*x2*x10*k2 + (-1)*x2*k3 + (-1)*x2*x4*k6 + 1*x5*k7 + 1*x8*k11 + (-1)*x2*x6*k12 + 1*x7*k13)/k1,
  diff(x3, t) = (1*x2*x10*k2 + (-1)*x3*k4)/k1,
  diff(x4, t) = ((-1)*x2*x4*k6 + 1*x5*k7)/k1,
  diff(x5, t) = (1*x2*x4*k6 + (-1)*x5*k7)/k1,
  diff(x6, t) = ((-1)*x2*x6*k12 + 1*x7*k13)/k1,
  diff(x7, t) = (1*x2*x6*k12 + (-1)*x7*k13)/k1,
  diff(x8, t) = (1*x4*k8 + 1*x5*k9 + (-1)*x8*k10)/k1,
  diff(x9, t) = (1*x6*k14 + 1*x7*k15 + (-1)*x9*k16)/k1,
  diff(x10, t) = ((-1)*x2*x10*k2 + 1*x3*k4 + (-1)*x10*k5 + 1*x9*k17)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file