load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 1;
k2 := 57/10000;
k3 := 39/2500;
k4 := 31/2000;
k5 := 397/5000;
k6 := 907/10000;
k7 := 137/5000;
k8 := 17/80;
k9 := 0;
k10 := 3817/2000;
k11 := 707/10000;
k12 := 1131/10000;
k13 := 1/1250;
k14 := 11/5000;
k15 := 17/5000;
k16 := 159/10000;
k17 := 67/5000;
k18 := 9;
operator diff$
odes := { diff(x1, t) = ((-1)*k2*x1 + (-1)*k3*x1 + (-1)*k4*x1)/k1,
  diff(x2, t) = (1*k2*x1 + (-1)*k5*x2 + (-1)*k11*x2 + (-1)*k12*x2)/k1,
  diff(x3, t) = (1*k3*x1 + (-1)*k8*x3 + (-1)*k17*x3)/k1,
  diff(x4, t) = (1*k4*x1 + 1*k6*x6 + 1*k9*x8 + (-1)*k15*x4*x5 + (-1)*k16*x4)/k1,
  diff(x5, t) = (1*k4*x1 + 1*k5*x2 + 1*k8*x3 + 1*k11*x2 + 1*k12*x2 + (-1)*k15*x4*x5 + 1*k17*x3)/k1,
  diff(x6, t) = (1*k5*x2 + (-1)*k6*x6 + (-1)*k7*x6 + 1*k14*x11)/k1,
  diff(x7, t) = (1*k7*x6 + 1*k16*x4)/k1,
  diff(x8, t) = (1*k8*x3 + (-1)*k9*x8 + (-1)*k10*x8)/k1,
  diff(x9, t) = (1*k10*x8 + 1*k16*x4)/k1,
  diff(x10, t) = (1*k11*x2 + (-1)*k13*x10)/k1,
  diff(x11, t) = (1*k12*x2 + 1*k13*x10 + (-1)*k14*x11)/k1,
  diff(x12, t) = 1*k16*x4/k1,
  diff(x13, t) = 1*k15*x4*x5/k1,
  diff(x14, t) = 1*k17*x3/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file