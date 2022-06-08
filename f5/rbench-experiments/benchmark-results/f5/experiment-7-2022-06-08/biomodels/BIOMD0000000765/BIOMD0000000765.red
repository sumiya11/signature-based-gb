load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 10;
k3 := 1/10;
k4 := 1;
k5 := 0;
k6 := 283/500;
k7 := 0;
k8 := 1/10;
k9 := 0;
k10 := 0;
k11 := 400;
k12 := 400;
k13 := 10;
k14 := 1;
k15 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k14*k3*x3*x1 + 1*k14*k8*x4 + (-1)*k14*(k1 + k10)*x1 + 1*k9*x2/k2)/k14,
  diff(x2, t) = (1*k10*x1*k2 + (-1)*k15*k9*x2)/k15,
  diff(x3, t) = ((-1)*k14*k3*x3*x1 + 1*k14*k8*x4 + 1*k14*k7 + (-1)*k14*k6*x3)/k14,
  diff(x4, t) = (1*k14*k3*x3*x1 + (-1)*k14*(k8 + k5)*x4)/k14 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file