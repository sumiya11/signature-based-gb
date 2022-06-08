load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 1/100;
k3 := 509/100000;
k4 := 47/100000;
k5 := 11/10000;
k6 := 89/12500;
k7 := 439/100000;
k8 := 9/50000;
k9 := 5567/50000;
k10 := 14359/100000;
k11 := 3/20000;
k12 := 6257/50000;
k13 := 15;
operator diff$
odes := { diff(x1, t) = ((-1)*k2*x1 + 1*k3*x2 + (-1)*k4*x1 + (-1)*k8*x1*x10)/k1,
  diff(x2, t) = (1*k2*x1 + (-1)*k3*x2 + (-1)*k5*x2 + (-1)*k6*x2 + (-1)*k11*x2*x10)/k1,
  diff(x3, t) = (1*k4*x1 + 1*k5*x2)/k1,
  diff(x4, t) = (2*k6*x2 + (-1)*k7*x4)/k1,
  diff(x5, t) = (1*k7*x4 + 1*k9*x7)/k1,
  diff(x6, t) = 1*k7*x4/k1,
  diff(x7, t) = (1*k8*x1*x10 + (-1)*k9*x7 + (-1)*k10*x7)/k1,
  diff(x8, t) = (1*k10*x7 + 1*k11*x2*x10 + (-1)*k12*x8)/k1,
  diff(x9, t) = (1*k4*x1 + 1*k5*x2)/k1,
  diff(x10, t) = ((-1)*k8*x1*x10 + 1*k9*x7 + (-1)*k11*x2*x10)/k1,
  diff(x11, t) = 1*k12*x8/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file