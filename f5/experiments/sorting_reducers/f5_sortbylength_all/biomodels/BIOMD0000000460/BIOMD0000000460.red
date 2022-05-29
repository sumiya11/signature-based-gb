load_package f5$
torder({}, revgradlex)$

k1 := 29/5000000000;
k2 := 100;
k3 := 21/250000;
k4 := 3/25000000000000;
k5 := 1/500000;
k6 := 13/250;
k7 := 17/10000000;
k8 := 1;
k9 := 100;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (1*(k9*k2 - x2*(k1 + k3*x4)) + (-1)*((-k6*x3) + k7*x2) + (-1)*(k5*x2 - k4*x4))/k8,
  diff(x3, t) = 1*((-k6*x3) + k7*x2)/k8,
  diff(x4, t) = 1*(k5*x2 - k4*x4)/k8 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file