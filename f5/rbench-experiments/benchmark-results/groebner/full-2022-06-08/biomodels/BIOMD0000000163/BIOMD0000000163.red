load_package groebner$
torder({}, revgradlex)$

k1 := 103/10000;
k2 := 2869/100000;
k3 := 33/100;
k4 := 33/1000;
k5 := 33/100;
k6 := 1871/50000;
k7 := 1/200;
k8 := 2197;
k9 := 2609/100000;
k10 := 1/200;
k11 := 1/40;
k12 := 4/25;
k13 := 1;
k14 := 2/25;
k15 := 1/2;
k16 := 137/2000000;
k17 := 4/25;
k18 := 587/5000;
k19 := 4729/20;
k20 := 49261/100;
k21 := 1;
k22 := 7/20000;
k23 := 21/20000;
k25 := 4729/20;
k24 := 49261/100;
k27 := 13793/25;
k26 := 5747/5;
operator diff$
odes := { diff(x1, t) = ((-1)*k23*k12*x1 + 1*k22*k13*x2 + (-1)*k23*k16*x1*x3*x13)/k23,
  diff(x2, t) = (1*k23*k12*x1 + (-1)*k22*k13*x2 + 1*k22*k18*x15)/k22,
  diff(x3, t) = ((-1)*k23*k14*x3 + 1*k22*k15*x4 + (-1)*k23*k16*x1*x3*x13)/k23,
  diff(x4, t) = (1*k23*k14*x3 + (-1)*k22*k15*x4 + 1*k22*k18*x15)/k22,
  diff(x5, t) = (1*k23*k1 + (-1)*k23*k5*x5 + 1*k23*k6*x6 + (-1)*k23*k3*x5 + 1*k23*k4*x7 + (-1)*k23*k8*x16*x8*x5 + 1*k23*k6*x12 + 1*k23*k4*x13)/k23,
  diff(x6, t) = (1*k23*k5*x5 + (-1)*k23*k6*x6)/k23,
  diff(x7, t) = (1*k23*k3*x5 + (-1)*k23*k4*x7 + (-1)*k23*k10*x7)/k23,
  diff(x8, t) = (1*k23*k2 + (-1)*k23*k5*x8 + 1*k23*k6*x9 + (-1)*k23*k3*x8 + 1*k23*k4*x10 + (-1)*k23*k8*x16*x8*x5 + 1*k23*k6*x12 + 1*k23*k4*x13)/k23,
  diff(x9, t) = (1*k23*k5*x8 + (-1)*k23*k6*x9)/k23,
  diff(x10, t) = (1*k23*k3*x8 + (-1)*k23*k4*x10 + (-1)*k23*k11*x10)/k23,
  diff(x11, t) = (1*k23*k8*x16*x8*x5 + (-1)*k23*k5*x11 + (-1)*k23*k3*x11)/k23,
  diff(x12, t) = (1*k23*k5*x11 + (-1)*k23*k6*x12 + (-1)*k23*k9*x12*x15)/k23,
  diff(x13, t) = (1*k23*k3*x11 + (-1)*k23*k4*x13 + (-1)*k23*k7*x13)/k23,
  diff(x14, t) = (1*k23*k16*x1*x3*x13 + (-1)*k23*k17*x14)/k23,
  diff(x15, t) = (1*k23*k17*x14 + (-1)*k22*k18*x15)/k22,
  diff(x16, t) = (1*k23*k4*x10 + (-1)*k23*k8*x16*x8*x5 + 1*k23*k6*x12 + 1*k23*k4*x13)/k21 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file