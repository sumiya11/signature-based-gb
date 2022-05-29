load_package groebner$
torder({}, revgradlex)$

parameters := k1 := 0;
k2 := 0;
k3 := 1;
k4 := 8780000;
k5 := 8;
k6 := 529000;
k7 := 419/50000000;
k8 := 636000000;
k9 := 179/10000;
k10 := 13/1000;
k11 := 903/1000000000;
k12 := 1/10000;
k13 := 623/10;
k14 := 853000;
k15 := 117/25000;
k16 := 132000000;
k17 := 32/25;
k18 := 386000;
k19 := 51/1250;
k20 := 64100;
k21 := 19/20;
k22 := 94700000;
k23 := 227/100000;
k24 := 13/1000;
k25 := 111/50000000000;
k26 := 2;
k27 := 1470000;
k28 := 44700;
k29 := 13/156250000;
k30 := 25;
k31 := 61/250;
k32 := 1/10000;
k33 := 383/100;
k34 := 74300;
k35 := 143/25000;
k36 := 22800000;
k37 := 543/10000000;
k38 := 1620000;
k39 := 7/800;
k40 := 6200000;
k41 := 433/10000;
k42 := 6300000;
k43 := 239/500;
k44 := 25;
k45 := 297/100000;
k46 := 13000;
k47 := 137/200;
k48 := 49400000;
k49 := 421/100000;
k50 := 11/4;
k51 := 2940;
k52 := 0;
k53 := 1/100000000;
k54 := 617/1000000;
k55 := 1217/250000;
k56 := 1/1000000;
operator diff$
odes := { diff(x1, t) = ((-1)*k3*(k4*x2*x1 - k5*x3) + (-1)*k3*(k18*x6*x1 - k19*x12) + (-1)*k3*(k20*x7*x1 - k21*x13) + (-1)*k3*(k34*x10*x1 - k35*x15) + (-1)*k3*(k42*x11*x1 - k43*x16) + (-1)*k3*(k46*x14*x1 - k47*x17))/k3,
  diff(x2, t) = ((-1)*k3*(k4*x2*x1 - k5*x3) + (-1)*k3*(k6*x2*x4 - k7*x6) + (-1)*k3*(k8*x2*x5 - k9*x10) + 1*k3*(k12*x7 - k13*x2*x9))/k3,
  diff(x3, t) = (1*k3*(k4*x2*x1 - k5*x3) + (-1)*k3*(k28*x3*x4 - k29*x12) + 1*k3*(k32*x13 - k33*x3*x9) + (-1)*k3*(k36*x3*x5 - k37*x15))/k3,
  diff(x4, t) = ((-1)*k3*(k6*x2*x4 - k7*x6) + (-1)*k3*(k14*x10*x4 - k15*x11) + (-1)*k3*(k28*x3*x4 - k29*x12) + (-1)*k3*(k38*x15*x4 - k39*x16))/k3,
  diff(x5, t) = ((-1)*k3*(k8*x2*x5 - k9*x10) + (-1)*k3*(k16*x6*x5 - k17*x11) + (-1)*k3*(k22*x7*x5 - k23*x14) + (-1)*k3*(k36*x3*x5 - k37*x15) + (-1)*k3*(k40*x12*x5 - k41*x16) + (-1)*k3*(k48*x13*x5 - k49*x17))/k3,
  diff(x6, t) = (1*k3*(k6*x2*x4 - k7*x6) + (-1)*k3*(k10*x6 - k11*x7*x8) + (-1)*k3*(k16*x6*x5 - k17*x11) + (-1)*k3*(k18*x6*x1 - k19*x12))/k3,
  diff(x7, t) = (1*k3*(k10*x6 - k11*x7*x8) + (-1)*k3*(k12*x7 - k13*x2*x9) + (-1)*k3*(k20*x7*x1 - k21*x13) + (-1)*k3*(k22*x7*x5 - k23*x14))/k3,
  diff(x8, t) = (1*k3*(k10*x6 - k11*x7*x8) + 1*k3*(k24*x11 - k25*x14*x8) + 1*k3*(k30*x12 - k31*x13*x8) + 1*k3*(k44*x16 - k45*x17*x8))/k3,
  diff(x9, t) = (1*k3*(k12*x7 - k13*x2*x9) + 1*k3*(k26*x14 - k27*x10*x9) + 1*k3*(k32*x13 - k33*x3*x9) + 1*k3*(k50*x17 - k51*x15*x9))/k3,
  diff(x10, t) = (1*k3*(k8*x2*x5 - k9*x10) + (-1)*k3*(k14*x10*x4 - k15*x11) + 1*k3*(k26*x14 - k27*x10*x9) + (-1)*k3*(k34*x10*x1 - k35*x15))/k3,
  diff(x11, t) = (1*k3*(k14*x10*x4 - k15*x11) + 1*k3*(k16*x6*x5 - k17*x11) + (-1)*k3*(k24*x11 - k25*x14*x8) + (-1)*k3*(k42*x11*x1 - k43*x16))/k3,
  diff(x12, t) = (1*k3*(k18*x6*x1 - k19*x12) + 1*k3*(k28*x3*x4 - k29*x12) + (-1)*k3*(k30*x12 - k31*x13*x8) + (-1)*k3*(k40*x12*x5 - k41*x16))/k3,
  diff(x13, t) = (1*k3*(k20*x7*x1 - k21*x13) + 1*k3*(k30*x12 - k31*x13*x8) + (-1)*k3*(k32*x13 - k33*x3*x9) + (-1)*k3*(k48*x13*x5 - k49*x17))/k3,
  diff(x14, t) = (1*k3*(k22*x7*x5 - k23*x14) + 1*k3*(k24*x11 - k25*x14*x8) + (-1)*k3*(k26*x14 - k27*x10*x9) + (-1)*k3*(k46*x14*x1 - k47*x17))/k3,
  diff(x15, t) = (1*k3*(k34*x10*x1 - k35*x15) + 1*k3*(k36*x3*x5 - k37*x15) + (-1)*k3*(k38*x15*x4 - k39*x16) + 1*k3*(k50*x17 - k51*x15*x9))/k3,
  diff(x16, t) = (1*k3*(k38*x15*x4 - k39*x16) + 1*k3*(k40*x12*x5 - k41*x16) + 1*k3*(k42*x11*x1 - k43*x16) + (-1)*k3*(k44*x16 - k45*x17*x8))/k3,
  diff(x17, t) = (1*k3*(k44*x16 - k45*x17*x8) + 1*k3*(k46*x14*x1 - k47*x17) + 1*k3*(k48*x13*x5 - k49*x17) + (-1)*k3*(k50*x17 - k51*x15*x9))/k3 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file