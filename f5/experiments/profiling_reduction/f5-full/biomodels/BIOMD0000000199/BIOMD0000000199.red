load_package f5$
torder({}, revgradlex)$

k1 := 13/5;
k2 := 9/10;
k3 := 26;
k4 := 13/5;
k5 := 9/10;
k6 := 26;
k7 := 5;
k8 := 13/5;
k9 := 1/10000;
k10 := 13/10000;
k11 := 1;
k12 := 40;
k13 := 140;
k14 := 0;
k15 := 0;
k16 := 0;
k17 := 0;
k18 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k11*k1*x1 + 1*k11*k7*x7 + 1*k11*k10*x8*k13)/k11,
  diff(x2, t) = (1*k11*k1*x1 + (-1)*k11*k2*x2*k13 + 1*k11*k9*x8)/k11,
  diff(x3, t) = (1*k11*k2*x2*k13 + (-1)*k11*k3*x3)/k11,
  diff(x4, t) = (1*k11*k3*x3 + (-1)*k11*k4*x4)/k11,
  diff(x5, t) = (1*k11*k4*x4 + (-1)*k11*k5*x5*k13)/k11,
  diff(x6, t) = (1*k11*k5*x5*k13 + (-1)*k11*k6*x6)/k11,
  diff(x7, t) = (1*k11*k6*x6 + (-1)*k11*k7*x7 + (-1)*k11*k8*x7)/k11,
  diff(x8, t) = (1*k11*k8*x7 + (-1)*k11*k9*x8 + (-1)*k11*k10*x8*k13)/k11,
  diff(x9, t) = 0,
  diff(x10, t) = 0,
  diff(x11, t) = 0,
  diff(x12, t) = 0,
  diff(x13, t) = 0,
  diff(x14, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file