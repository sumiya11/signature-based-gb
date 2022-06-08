load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 0;
k2 := 100;
k3 := 1309/2500000000000000;
k4 := 14399/125000000000000000;
k5 := 51051/5000000000000000;
k6 := 282741/5000000000000000000000;
k7 := 1309/62500000000000;
k8 := 14399/5000000000000000;
k9 := 43773/200000000000000;
k10 := 314159/500000000000000000;
k11 := 121737/5000000000000000;
k12 := 33301/10000000000000000;
k13 := 314159/500000000000000000;
k14 := 2;
k15 := 1;
k16 := 62496021135727/6250000000000;
operator diff$
odes := { diff(x1, t) = 1*k3*(k7*x5 - k6*x1*x6)/k3/k3,
  diff(x2, t) = (-1)*k3*(k13*x2*x8 - k12*x3)/k3/k3,
  diff(x3, t) = ((-1)*k3*(k9*x3 - k8*x4)/k3 + 1*k3*(k13*x2*x8 - k12*x3)/k3)/k3,
  diff(x4, t) = ((-1)*k3*(k5*x4 - k4*x5)/k3 + 1*k3*(k9*x3 - k8*x4)/k3)/k3,
  diff(x5, t) = (1*k3*(k5*x4 - k4*x5)/k3 + (-1)*k3*(k7*x5 - k6*x1*x6)/k3)/k3,
  diff(x6, t) = (1*k3*(k7*x5 - k6*x1*x6)/k3 + (-1)*k3*(k11*x6 - k10*x7*x8)/k3)/k3,
  diff(x7, t) = 1*k3*(k11*x6 - k10*x7*x8)/k3/k3,
  diff(x8, t) = (1*k3*(k11*x6 - k10*x7*x8)/k3 + (-1)*k3*(k13*x2*x8 - k12*x3)/k3)/k3 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file