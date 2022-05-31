load_package f5$
torder({}, revgradlex)$

k1 := 300000000;
k2 := 8000;
k3 := 150000000;
k4 := 16000;
k5 := 30000;
k6 := 700;
k7 := 300000000;
k8 := 216/25;
k9 := 150000000;
k10 := 432/25;
k11 := 27/50;
k12 := 10800;
k13 := 130;
k14 := 2740;
k15 := 300000000;
k16 := 4;
k17 := 150000000;
k18 := 8;
k19 := 197/10;
k20 := 187/50;
k21 := 397/20;
k22 := 87/50;
k23 := 20;
k24 := 81/100;
k25 := 300000000;
k26 := 4;
k27 := 150000000;
k28 := 8;
k29 := 1/20;
k30 := 3/2500;
k31 := 1/20;
k32 := 3/2500;
k33 := 1/20;
k34 := 3/2500;
k35 := 1/10000000000000000;
k36 := 1/10000000000000000000000;
k37 := 1/1000000000000000000000;
operator diff$
odes := { diff(x1, t) = (1*k35*(k3*x5*x13 - k4*x1) + (-1)*k35*(k5*x1 - k6*x12))/k35,
  diff(x2, t) = (1*k35*(k15*x11*x13 - k16*x2) + (-1)*k35*(k17*x2*x13 - k18*x9) + 1*k35*(k21*x3 - k22*x2) + (-1)*k35*(k31*x2 - k32*x10))/k35,
  diff(x3, t) = (1*k35*(k7*x4*x13 - k8*x3) + (-1)*k35*(k9*x3*x13 - k10*x12) + 1*k35*(k13*x5 - k14*x3) + (-1)*k35*(k21*x3 - k22*x2))/k35,
  diff(x4, t) = ((-1)*k35*(k7*x4*x13 - k8*x3) + 1*k35*(k11*x6 - k12*x4) + (-1)*k35*(k19*x4 - k20*x11))/k35,
  diff(x5, t) = (1*k35*(k1*x6*x13 - k2*x5) + (-1)*k35*(k3*x5*x13 - k4*x1) + (-1)*k35*(k13*x5 - k14*x3))/k35,
  diff(x6, t) = ((-1)*k35*(k1*x6*x13 - k2*x5) + (-1)*k35*(k11*x6 - k12*x4))/k35,
  diff(x7, t) = (1*k35*(k27*x10*x13 - k28*x7) + 1*k35*(k33*x9 - k34*x7))/k35,
  diff(x8, t) = ((-1)*k35*(k25*x8*x13 - k26*x10) + 1*k35*(k29*x11 - k30*x8))/k35,
  diff(x9, t) = (1*k35*(k17*x2*x13 - k18*x9) + 1*k35*(k23*x12 - k24*x9) + (-1)*k35*(k33*x9 - k34*x7))/k35,
  diff(x10, t) = (1*k35*(k25*x8*x13 - k26*x10) + (-1)*k35*(k27*x10*x13 - k28*x7) + 1*k35*(k31*x2 - k32*x10))/k35,
  diff(x11, t) = ((-1)*k35*(k15*x11*x13 - k16*x2) + 1*k35*(k19*x4 - k20*x11) + (-1)*k35*(k29*x11 - k30*x8))/k35,
  diff(x12, t) = (1*k35*(k5*x1 - k6*x12) + 1*k35*(k9*x3*x13 - k10*x12) + (-1)*k35*(k23*x12 - k24*x9))/k35,
  diff(x13, t) = ((-1)*k35*(k1*x6*x13 - k2*x5) + (-1)*k35*(k3*x5*x13 - k4*x1) + (-1)*k35*(k7*x4*x13 - k8*x3) + (-1)*k35*(k9*x3*x13 - k10*x12) + (-1)*k35*(k15*x11*x13 - k16*x2) + (-1)*k35*(k17*x2*x13 - k18*x9) + (-1)*k35*(k25*x8*x13 - k26*x10) + (-1)*k35*(k27*x10*x13 - k28*x7))/k35 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file