load_package groebner$
torder({}, revgradlex)$

k1 := 21/250;
k2 := 4/125;
k3 := 2;
k4 := 1;
k5 := 7/200;
k6 := 40;
k7 := 9/10;
k8 := 0;
k9 := 2/125;
k10 := 1/50;
k11 := 2/625;
k12 := 1/100;
k13 := 1/100;
k14 := 0;
k15 := 0;
operator diff$
odes := { diff(x1, t) = ((-1)*k4*k5*x3*x1 + (-1)*k4*k7*x4*x1 + (-1)*k4*k2*x1 + 1*k4*k12)/k4,
  diff(x2, t) = (1*k4*k5*x3*x1 + 1*k4*k7*x4*x1 + (-1)*k4*k1*x2 + (-1)*k4*(k8*x2*k14 - k9*k15) + 1*k4*k10*k15)/k4,
  diff(x3, t) = ((-1)*k4*k6*x2^k3*x3 + (-1)*k4*k1*x3 + 1*k4*k11)/k4,
  diff(x4, t) = (1*k4*k6*x2^k3*x3 + (-1)*k4*k1*x4)/k4,
  diff(x5, t) = 0,
  diff(x6, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file