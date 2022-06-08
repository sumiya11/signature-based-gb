load_package groebner$
torder({}, revgradlex)$

k1 := 53/100;
k2 := 9/1250;
k3 := 5/8;
k4 := 49/20000;
k5 := 63/2000;
k6 := 4/5;
k7 := 3/400;
k8 := 71/1000;
k9 := 23/25;
k10 := 61/50000;
k11 := 87/100;
k12 := 1;
k13 := 2;
k14 := 5/2;
k15 := 5/2;
k16 := 5/2;
k17 := 3;
operator diff$
odes := { diff(x1, t) = ((-1)*k12*k1*x1*x2 + 1*k12*k2*x3 + 1*k12*k5*x4)/k12,
  diff(x2, t) = ((-1)*k12*k1*x1*x2 + 1*k12*k2*x3 + 1*k12*k11*x11)/k12,
  diff(x3, t) = (1*k12*k1*x1*x2 + (-1)*k12*k2*x3 + (-1)*k12*k3*x3*x9 + 1*k12*k4*x4)/k12,
  diff(x4, t) = (1*k12*k3*x3*x9 + (-1)*k12*k4*x4 + (-1)*k12*k5*x4)/k12,
  diff(x5, t) = (1*k12*k5*x4 + (-1)*k12*k6*x5*x7 + 1*k12*k7*x8)/k12,
  diff(x6, t) = (1*k12*k5*x4 + (-1)*k12*k9*x6*x10 + 1*k12*k10*x11)/k12,
  diff(x7, t) = ((-1)*k12*k6*x5*x7 + 1*k12*k7*x8 + 1*k12*k8*x8)/k12,
  diff(x8, t) = (1*k12*k6*x5*x7 + (-1)*k12*k7*x8 + (-1)*k12*k8*x8)/k12,
  diff(x9, t) = ((-1)*k12*k3*x3*x9 + 1*k12*k4*x4 + 1*k12*k8*x8)/k12,
  diff(x10, t) = ((-1)*k12*k9*x6*x10 + 1*k12*k10*x11 + 1*k12*k11*x11)/k12,
  diff(x11, t) = (1*k12*k9*x6*x10 + (-1)*k12*k10*x11 + (-1)*k12*k11*x11)/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file