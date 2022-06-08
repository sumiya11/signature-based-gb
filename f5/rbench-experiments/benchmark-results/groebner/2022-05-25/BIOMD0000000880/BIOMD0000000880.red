load_package groebner$
torder({}, revgradlex)$

k1 := 431/1000;
k2 := 117/500;
k3 := 17/1000;
k4 := 51/50000000000;
k5 := 1/80000;
k6 := 1/2000;
k7 := 13/1250000000;
k8 := 21/500000000;
k9 := 1/1250000000;
k10 := 1/50000;
k11 := 1/50000;
k12 := 1/50000;
k13 := 0;
k14 := 0;
k15 := 1;
operator diff$
odes := { diff(x1, t) = (1*k15*k1*x1*(1 - k4*x1) + (-1)*k15*k8*x1*x2 + (-1)*k15*k13*x1)/k15,
  diff(x2, t) = (1*k15*k10*(x3 + x4)*x2 + (-1)*k15*k9*x1*x2 + (-1)*k15*k7*x2)/k15,
  diff(x3, t) = (1*k15*k2*x3*(1 - k5*x3) + (-1)*k15*k11*x3*x2 + 1*k15*k14)/k15,
  diff(x4, t) = (1*k15*k3*x4*(1 - k6*x4) + (-1)*k15*k12*x4*x2)/k15 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file