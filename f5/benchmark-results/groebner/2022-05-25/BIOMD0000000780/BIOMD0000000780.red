load_package groebner$
torder({}, revgradlex)$

k1 := 1/50;
k2 := 1/50;
k3 := 1/5;
k4 := 2/5;
k5 := 4/5;
k6 := 1/100;
k7 := 4/5;
k8 := 1/25;
k9 := 1/10;
k10 := 1/1000;
k11 := 1/2;
k12 := 1/100;
k13 := 16/3;
k14 := 16/3;
k15 := 0;
k16 := 1;
operator diff$
odes := { diff(x1, t) = (1*k16*k1 + (-1)*k16*k2*x1 + (-1)*k16*k3*x1*x2 + (-1)*k16*k4*x1*x3)/k16,
  diff(x2, t) = (1*k16*k5*k3*x1*x2 + (-1)*k16*(k2 + k6)*x2)/k16,
  diff(x3, t) = (1*k16*k7*k4*x1*x3 + (-1)*k16*(k2 + k8)*x3 + (-1)*k16*k9*x3*x4)/k16,
  diff(x4, t) = (1*k16*k10 + 1*k16*k11*k9*x3*x4 + (-1)*k16*(k2 + k12)*x4)/k16 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file