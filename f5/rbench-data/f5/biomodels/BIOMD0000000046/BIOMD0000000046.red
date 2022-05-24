load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 3/1000000;
k3 := 18;
k4 := 3/20;
k5 := 13/2500;
k6 := 20;
k7 := 17;
k8 := 20;
k9 := 40;
k10 := 60;
k11 := 9/5;
k12 := 1/10;
k13 := 2/25;
k14 := 3/500;
k15 := 3/500;
k16 := 7/10;
k17 := 0;
k18 := 0;
k19 := 12;
k20 := 0;
k21 := 7/5;
k22 := 500;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*k2*x1*x2 + 1*k1*k13 + (-1)*k1*k16*x8*x1)/k1,
  diff(x2, t) = ((-1)*k1*k2*x1*x2 + (-1)*k1*k6*x9*x2 + 1*k1*k8*x10*x10 + (-1)*k1*k12*x12*x2 + 1*k1*k14*k19 + (-1)*k1*k15*x2)/k1,
  diff(x3, t) = (1*k1*k2*x1*x2 + (-1)*k1*k3*x3*x4 + 1*k1*k8*x10*x10)/k1,
  diff(x4, t) = ((-1)*k1*k3*x3*x4 + 1*k1*k5*x7*x6 + (-1)*k1*k7*x10*x4 + (-1)*k1*k11*x4*x9)/k1,
  diff(x5, t) = (1*k1*k3*x3*x4 + (-1)*k1*k4*x5*x6 + 1*k1*k9*x11*x9)/k1,
  diff(x6, t) = ((-1)*k1*k4*x5*x6 + (-1)*k1*k5*x7*x6 + 1*k1*k16*x8*x1)/k1,
  diff(x7, t) = (1*k1*k4*x5*x6 + (-1)*k1*k5*x7*x6)/k1,
  diff(x8, t) = (1*k1*k4*x5*x6 + 1*k1*k5*x7*x6 + (-1)*k1*k16*x8*x1)/k1,
  diff(x9, t) = ((-1)*k1*k6*x9*x2 + (-1)*k1*k9*x11*x9 + (-2)*k1*k10*x9*x9 + (-1)*k1*k11*x4*x9 + 1*k1*k16*x8*x1)/k1,
  diff(x10, t) = (1*k1*k6*x9*x2 + (-1)*k1*k7*x10*x4 + (-2)*k1*k8*x10*x10)/k1,
  diff(x11, t) = (1*k1*k7*x10*x4 + (-1)*k1*k9*x11*x9 + 1*k1*k12*x12*x2)/k1,
  diff(x12, t) = (1*k1*k11*x4*x9 + (-1)*k1*k12*x12*x2)/k1,
  diff(x13, t) = 0,
  diff(x14, t) = 0,
  diff(x15, t) = 0,
  diff(x16, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file