load_package f5$
1;
torder({}, revgradlex)$

k1 := 289/10000;
k2 := 309/1000000;
k3 := 507/500;
k4 := 26553/1000000;
k5 := 18637/100000;
k6 := 44489/10000;
k7 := 53/500000;
k8 := 85341/100000000;
k9 := 289/10;
k10 := 5000/829;
k11 := 326/829;
k12 := 2000000000000;
k13 := 829;
k14 := 326;
k15 := 10;
k16 := 1;
operator diff$
odes := { diff(x1, t) = (-(k1 + k2))*x1 + k7*x2 + k3*k12/(k13*400000)*1/1000*x3,
  diff(x2, t) = k1*x1 - k7*x2,
  diff(x3, t) = k2*k12/(k13*400000)*1/1000*x1 - (k3 + k4)*x3 + k5*x4 - k8*(k9 - x5)*x3 + k14/k13*k6*x5,
  diff(x4, t) = k4*x3 - k5*x4,
  diff(x5, t) = k8*k14/k13*(k9 - x5)*x3 - k6*x5 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file