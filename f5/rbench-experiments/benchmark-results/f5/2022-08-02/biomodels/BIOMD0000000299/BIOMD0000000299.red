load_package f5$
1;
torder({}, revgradlex)$

k1 := 8/5;
k2 := 101/200;
k3 := 1/2;
k4 := 1;
k5 := 4;
k6 := 1/5;
k7 := 1/2;
k8 := 7/5;
k9 := 13/100;
k10 := 1/2;
k11 := 3/5;
k12 := 1;
k13 := 1;
operator diff$
odes := { diff(x1, t) = k1*k4^k5/(k4^k5 + x3^k5) - k2*x1/(k3 + x1),
  diff(x2, t) = k7*x1 + k11*x3 - (k8*x2/(k9 + x2) + k10*x2),
  diff(x3, t) = k10*x2 - k11*x3 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file