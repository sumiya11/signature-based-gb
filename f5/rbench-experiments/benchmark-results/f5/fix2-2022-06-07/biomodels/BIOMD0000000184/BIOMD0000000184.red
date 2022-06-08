load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 1/20;
k2 := 1/2;
k3 := 40;
k4 := 3/20;
k5 := 101/50;
k6 := 3/20;
k7 := 11/5;
k8 := 1/10;
k9 := 15;
k10 := 1/10;
k11 := 1/2;
k12 := 1/20;
k13 := 3/10;
k14 := 2/25;
k15 := 1;
k16 := 1;
operator diff$
odes := { diff(x1, t) = (1*k15*k1 + (-1)*k15*k2*x1 + 1*k16*4*k3*k4^k5*x1^k5/((x1^k5 + k4^k5)*(x1^k5 + k6^k5))*x3^k7/(x3^k7 + k8^k7)*(x2 - x1) + (-1)*k15*k9*x1^2/(x1^2 + k10^2) + 1*k16*k11*(x2 - x1))/k15,
  diff(x2, t) = ((-1)*k16*4*k3*k4^k5*x1^k5/((x1^k5 + k4^k5)*(x1^k5 + k6^k5))*x3^k7/(x3^k7 + k8^k7)*(x2 - x1) + 1*k15*k9*x1^2/(x1^2 + k10^2) + (-1)*k16*k11*(x2 - x1))/k16,
  diff(x3, t) = (1*k15*k12*x1^2/(x1^2 + k13^2) + (-1)*k15*k14*x3)/k15 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file