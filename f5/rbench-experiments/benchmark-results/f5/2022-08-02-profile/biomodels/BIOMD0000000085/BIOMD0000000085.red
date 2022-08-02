load_package f5$
torder({}, revgradlex)$

k1 := 0;
k2 := 0;
k3 := 1;
k4 := 529000;
k5 := 419/50000000;
k6 := 13/1000;
k7 := 903/1000000000;
k8 := 1/10000;
k9 := 623/10;
k10 := 853000;
k11 := 117/25000;
k12 := 132000000;
k13 := 32/25;
k14 := 386000;
k15 := 51/1250;
k16 := 64100;
k17 := 19/20;
k18 := 94700000;
k19 := 227/100000;
k20 := 13/1000;
k21 := 111/50000000000;
k22 := 2;
k23 := 1470000;
k24 := 25;
k25 := 61/250;
k26 := 1/10000;
k27 := 383/100;
k28 := 3960000000;
k29 := 543/10000000;
k30 := 1620000;
k31 := 7/800;
k32 := 6300000;
k33 := 239/500;
k34 := 25;
k35 := 297/100000;
k36 := 11/4;
k37 := 2940;
k38 := 0;
k39 := 1/100000000;
k40 := 617/1000000;
k41 := 1217/250000;
k42 := 1/100000;
operator diff$
odes := { diff(x1, t) = ((-1)*k3*(k14*x6*x1 - k15*x12) + (-1)*k3*(k16*x7*x1 - k17*x13) + (-1)*k3*(k32*x11*x1 - k33*x16))/k3,
  diff(x2, t) = ((-1)*k3*(k4*x2*x4 - k5*x6) + 1*k3*(k8*x7 - k9*x2*x9))/k3,
  diff(x3, t) = (1*k3*(k26*x13 - k27*x3*x9) + (-1)*k3*(k28*x3*x5 - k29*x15))/k3,
  diff(x4, t) = ((-1)*k3*(k4*x2*x4 - k5*x6) + (-1)*k3*(k10*x10*x4 - k11*x11) + (-1)*k3*(k30*x15*x4 - k31*x16))/k3,
  diff(x5, t) = ((-1)*k3*(k12*x6*x5 - k13*x11) + (-1)*k3*(k18*x7*x5 - k19*x14) + (-1)*k3*(k28*x3*x5 - k29*x15))/k3,
  diff(x6, t) = (1*k3*(k4*x2*x4 - k5*x6) + (-1)*k3*(k6*x6 - k7*x7*x8) + (-1)*k3*(k12*x6*x5 - k13*x11) + (-1)*k3*(k14*x6*x1 - k15*x12))/k3,
  diff(x7, t) = (1*k3*(k6*x6 - k7*x7*x8) + (-1)*k3*(k8*x7 - k9*x2*x9) + (-1)*k3*(k16*x7*x1 - k17*x13) + (-1)*k3*(k18*x7*x5 - k19*x14))/k3,
  diff(x8, t) = (1*k3*(k6*x6 - k7*x7*x8) + 1*k3*(k20*x11 - k21*x14*x8) + 1*k3*(k24*x12 - k25*x13*x8) + 1*k3*(k34*x16 - k35*x17*x8))/k3,
  diff(x9, t) = (1*k3*(k8*x7 - k9*x2*x9) + 1*k3*(k22*x14 - k23*x10*x9) + 1*k3*(k26*x13 - k27*x3*x9) + 1*k3*(k36*x17 - k37*x15*x9))/k3,
  diff(x10, t) = ((-1)*k3*(k10*x10*x4 - k11*x11) + 1*k3*(k22*x14 - k23*x10*x9))/k3,
  diff(x11, t) = (1*k3*(k10*x10*x4 - k11*x11) + 1*k3*(k12*x6*x5 - k13*x11) + (-1)*k3*(k20*x11 - k21*x14*x8) + (-1)*k3*(k32*x11*x1 - k33*x16))/k3,
  diff(x12, t) = (1*k3*(k14*x6*x1 - k15*x12) + (-1)*k3*(k24*x12 - k25*x13*x8))/k3,
  diff(x13, t) = (1*k3*(k16*x7*x1 - k17*x13) + 1*k3*(k24*x12 - k25*x13*x8) + (-1)*k3*(k26*x13 - k27*x3*x9))/k3,
  diff(x14, t) = (1*k3*(k18*x7*x5 - k19*x14) + 1*k3*(k20*x11 - k21*x14*x8) + (-1)*k3*(k22*x14 - k23*x10*x9))/k3,
  diff(x15, t) = (1*k3*(k28*x3*x5 - k29*x15) + (-1)*k3*(k30*x15*x4 - k31*x16) + 1*k3*(k36*x17 - k37*x15*x9))/k3,
  diff(x16, t) = (1*k3*(k30*x15*x4 - k31*x16) + 1*k3*(k32*x11*x1 - k33*x16) + (-1)*k3*(k34*x16 - k35*x17*x8))/k3,
  diff(x17, t) = (1*k3*(k34*x16 - k35*x17*x8) + (-1)*k3*(k36*x17 - k37*x15*x9))/k3 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file