load_package f5$
torder({}, revgradlex)$

parameters := k1 := 156/17;
k2 := 160/17;
k3 := 6;
k4 := 9;
k5 := 1/2;
k6 := 1/2;
k7 := 10;
k8 := 1;
k9 := 1;
k10 := 1;
k11 := 1/10;
k12 := 1;
k13 := 57/10;
k14 := 3/10;
k15 := 30;
k16 := 1/2;
k17 := 2;
k18 := 325;
k19 := 17/10;
k20 := 23/50;
k21 := 2;
k22 := 4;
k23 := 7/10;
k24 := 10;
operator diff$
odes := { diff(x1, t) = (1*k3*(1 + k4*x3^4/(k5^4 + x3^4))*x1^2/(x1^2 + k7/(1 + x3^4/k6^4)) + (-1)*k10*x1 + 1*k11)/k8,
  diff(x2, t) = (1*k15*x3^k17/(k16^k17 + x3^k17) + (-1)*k18*x2^k21/(k19^k21 + x2^k21)*x3^k22/(k20^k22 + x3^k22) + (-1)*k23*x2)/k9,
  diff(x3, t) = (1*k12 + 1*k13*k14 + (-1)*k15*x3^k17/(k16^k17 + x3^k17) + 1*k18*x2^k21/(k19^k21 + x2^k21)*x3^k22/(k20^k22 + x3^k22) + 1*k23*x2 + (-1)*k24*x3)/k8 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file