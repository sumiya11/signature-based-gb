load_package f5$
torder({}, revgradlex)$

k1 := 1/50;
k2 := 1;
k3 := 1/100;
k4 := 4/125;
k5 := 1;
k6 := 15;
k7 := 1/50;
k8 := 1;
k9 := 1/100;
k10 := 4/125;
k11 := 1;
k12 := 15;
k13 := 9/200;
k14 := 1;
k15 := 23/250;
k16 := 1;
k17 := 1/100;
k18 := 1/100;
k19 := 1;
k20 := 1/2;
k21 := 43/500;
k22 := 11/10000;
k23 := 1/100;
k24 := 1;
k25 := 47/100;
k26 := 7/50;
k27 := 9/5000;
k28 := 9/200;
k29 := 1;
k30 := 23/250;
k31 := 1;
k32 := 1/100;
k33 := 1;
k34 := 800;
k35 := 180;
k36 := 100;
operator diff$
odes := { diff(x1, t) = ((-1)*k33*(k1*x1*x5 - k2*x9) + (-1)*k33*(k7*x1*x5 - k8*x10) + 1*k33*(k21*x17 - k22*x1*x6) + 1*k33*(k26*x18 - k27*x1*x6))/k33,
  diff(x2, t) = (1*k33*k3*x9 + (-1)*k33*(k4*x2*x5 - k5*x7) + (-1)*k33*(k23*x2*x6 - k24*x13) + 1*k33*(k31*x14 - k32*x2*x6))/k33,
  diff(x3, t) = (1*k33*k9*x10 + (-1)*k33*(k10*x3*x5 - k11*x8) + 1*k33*(k16*x15 - k17*x3*x6) + (-1)*k33*(k18*x3*x6 - k19*x16))/k33,
  diff(x4, t) = (1*k33*k6*x7 + 1*k33*k12*x8 + (-1)*k33*(k13*x4*x6 - k14*x11) + (-1)*k33*(k28*x4*x6 - k29*x12))/k33,
  diff(x5, t) = ((-1)*k33*(k1*x1*x5 - k2*x9) + 1*k33*k3*x9 + (-1)*k33*(k4*x2*x5 - k5*x7) + 1*k33*k6*x7 + (-1)*k33*(k7*x1*x5 - k8*x10) + 1*k33*k9*x10 + (-1)*k33*(k10*x3*x5 - k11*x8) + 1*k33*k12*x8)/k33,
  diff(x6, t) = ((-1)*k33*(k13*x4*x6 - k14*x11) + 1*k33*(k16*x15 - k17*x3*x6) + (-1)*k33*(k18*x3*x6 - k19*x16) + 1*k33*(k21*x17 - k22*x1*x6) + (-1)*k33*(k23*x2*x6 - k24*x13) + 1*k33*(k26*x18 - k27*x1*x6) + (-1)*k33*(k28*x4*x6 - k29*x12) + 1*k33*(k31*x14 - k32*x2*x6))/k33,
  diff(x7, t) = (1*k33*(k4*x2*x5 - k5*x7) + (-1)*k33*k6*x7)/k33,
  diff(x8, t) = (1*k33*(k10*x3*x5 - k11*x8) + (-1)*k33*k12*x8)/k33,
  diff(x9, t) = (1*k33*(k1*x1*x5 - k2*x9) + (-1)*k33*k3*x9)/k33,
  diff(x10, t) = (1*k33*(k7*x1*x5 - k8*x10) + (-1)*k33*k9*x10)/k33,
  diff(x11, t) = (1*k33*(k13*x4*x6 - k14*x11) + (-1)*k33*k15*x11)/k33,
  diff(x12, t) = (1*k33*(k28*x4*x6 - k29*x12) + (-1)*k33*k30*x12)/k33,
  diff(x13, t) = (1*k33*(k23*x2*x6 - k24*x13) + (-1)*k33*k25*x13)/k33,
  diff(x14, t) = (1*k33*k30*x12 + (-1)*k33*(k31*x14 - k32*x2*x6))/k33,
  diff(x15, t) = (1*k33*k15*x11 + (-1)*k33*(k16*x15 - k17*x3*x6))/k33,
  diff(x16, t) = (1*k33*(k18*x3*x6 - k19*x16) + (-1)*k33*k20*x16)/k33,
  diff(x17, t) = (1*k33*k20*x16 + (-1)*k33*(k21*x17 - k22*x1*x6))/k33,
  diff(x18, t) = (1*k33*k25*x13 + (-1)*k33*(k26*x18 - k27*x1*x6))/k33 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file