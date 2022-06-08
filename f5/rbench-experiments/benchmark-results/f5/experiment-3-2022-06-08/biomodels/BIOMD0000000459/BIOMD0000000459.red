load_package f5$
torder({}, revgradlex)$

k1 := 41203/2500000;
k2 := 100;
k3 := 108559/100000000000000;
k4 := 932517/10000000000000;
k5 := 26923/200000000000;
k6 := 451769/500000000000;
k7 := 79443/25000000;
k8 := 1;
k9 := 100;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (1*(k9*k2 - k1*x2) + (-1)*(k6*x2 - x3*(k5 + k7*x4)) + (-1)*(k4*x2 - k3*x4))/k8,
  diff(x3, t) = 1*(k6*x2 - x3*(k5 + k7*x4))/k8,
  diff(x4, t) = 1*(k4*x2 - k3*x4)/k8 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file