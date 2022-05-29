load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 83/25000000000000000000;
k3 := 1/100;
k4 := 4;
k5 := 1/2500;
k6 := 1;
k7 := 1/100000;
k8 := 1/250;
k9 := 11/100;
k10 := 10000;
k11 := 10000;
operator diff$
odes := { diff(x1, t) = (-1)*k1*(k2*x1*x2 - k3*x7)/k1,
  diff(x2, t) = ((-1)*k1*(k2*x1*x2 - k3*x7) + 1*k1*k4 + (-1)*k1*k5*x2)/k1,
  diff(x3, t) = (1*k1*k6*x5*x4 + (-1)*k1*k7*x7*x3)/k1,
  diff(x4, t) = ((-1)*k1*k6*x5*x4 + 1*k1*k7*x7*x3)/k1,
  diff(x5, t) = ((-1)*k1*k6*x5*x4 + 1*k1*k9*x6)/k1,
  diff(x6, t) = (1*k1*k7*x7*x3 + (-1)*k1*k9*x6)/k1,
  diff(x7, t) = (1*k1*(k2*x1*x2 - k3*x7) + (-1)*k1*k8*x7)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file