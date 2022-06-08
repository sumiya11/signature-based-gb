load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 100;
k4 := 100;
k5 := 1/100000;
k6 := 1/100000;
k7 := 1/10;
k8 := 1/10;
k9 := 1;
k10 := 0;
k11 := 0;
k12 := 0;
k13 := 0;
operator diff$
odes := { diff(x1, t) = (1*k9*k3*x3 + (-1)*k9*k1*x1 + (-2)*k9*k6*x1*(x1 - 1)/2*x4 + 2*k9*k8*x6)/k9,
  diff(x2, t) = (1*k9*k4*x4 + (-1)*k9*k2*x2 + (-2)*k9*k5*x2*(x2 - 1)/2*x3 + 2*k9*k7*x5)/k9,
  diff(x3, t) = ((-1)*k9*k5*x2*(x2 - 1)/2*x3 + 1*k9*k7*x5)/k9,
  diff(x4, t) = ((-1)*k9*k6*x1*(x1 - 1)/2*x4 + 1*k9*k8*x6)/k9,
  diff(x5, t) = (1*k9*k5*x2*(x2 - 1)/2*x3 + (-1)*k9*k7*x5)/k9,
  diff(x6, t) = (1*k9*k6*x1*(x1 - 1)/2*x4 + (-1)*k9*k8*x6)/k9,
  diff(x7, t) = ((-1)*k9*k3*x3 + 1*k9*k1*x1)/k9,
  diff(x8, t) = ((-1)*k9*k4*x4 + 1*k9*k2*x2)/k9 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file