load_package f5$
torder({}, revgradlex)$

parameters := k1 := 10;
k2 := 1/50000;
k3 := 50;
k4 := 1/100000;
k5 := 1/250000;
k6 := 3/5000000;
k7 := 1/10000000;
k8 := 500;
k9 := 1;
k10 := 1/100;
k11 := 100;
k12 := 1/2;
k13 := 1/2;
k14 := 1/20;
k15 := 2/25;
k16 := 1000;
k17 := 401/50000000000;
k18 := 12;
k19 := 1/50;
k20 := 1/10;
k21 := 1/1000;
k22 := 1;
k23 := 0;
k24 := 6000;
k25 := 1;
k26 := 11000;
operator diff$
odes := { diff(x1, t) = ((-1)*k3*x4*x1 + 1*k4*x5 + 1*k5*x5*x14 + (-1)*k8*x1*x3 + 1*k9*x2 + 1*k16*x11 + (-1)*k17*x1*x14)/k22,
  diff(x2, t) = (1*k8*x1*x3 + (-1)*k9*x2)/k22,
  diff(x3, t) = ((-1)*k8*x1*x3 + 1*k9*x2 + (-2)*(x3 - 1)*k10*x3/2 + (-1)*k11*x3*x7 + 1*k12*x6 + 2*k13*x7)/k22,
  diff(x4, t) = (1*k2*x8*x13 + (-1)*k3*x4*x1 + 1*k4*x5 + (-1)*k6*x4*x14 + (-2)*(x4 - 1)*k7*x4/2 + (-1)*k7*x4*x9)/k22,
  diff(x5, t) = (1*k3*x4*x1 + (-1)*k4*x5 + (-1)*k5*x5*x14)/k22,
  diff(x6, t) = (1*k11*x3*x7 + (-1)*k12*x6 + (-1)*k14*x10*x6 + 1*k15*x11)/k22,
  diff(x7, t) = (1*(x3 - 1)*k10*x3/2 + (-1)*k11*x3*x7 + 1*k12*x6 + (-1)*k13*x7)/k22,
  diff(x8, t) = (1*k1 + (-1)*k2*x8*x13 + 1*k5*x5*x14)/k22,
  diff(x9, t) = (1*(x4 - 1)*k7*x4/2 + 1*k7*x4*x9)/k22,
  diff(x10, t) = ((-1)*k14*x10*x6 + 1*k15*x11)/k22,
  diff(x11, t) = (1*k14*x10*x6 + (-1)*k15*x11)/k22,
  diff(x12, t) = 1*1/k22,
  diff(x13, t) = (1*k20 + (-1)*k21*x13)/k22,
  diff(x14, t) = ((-1)*k5*x5*x14 + (-1)*k6*x4*x14 + (-1)*k17*x1*x14 + 1*k18*x15 + (-1)*k19*x14)/k22,
  diff(x15, t) = (1*k5*x5*x14 + 1*k6*x4*x14 + 1*k17*x1*x14 + (-1)*k18*x15 + 1*k19*x14)/k22,
  diff(x16, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file