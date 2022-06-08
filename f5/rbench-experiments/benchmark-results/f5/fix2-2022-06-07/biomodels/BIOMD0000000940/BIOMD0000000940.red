load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 1/2;
k2 := 1/200;
k3 := 1/10;
k4 := 1/2;
k5 := 1/5;
k6 := 3/5;
k7 := 2;
k8 := 5;
k9 := 1/5;
k10 := 3/5;
k11 := 3/5;
k12 := 1/10;
k13 := 2;
k14 := 3;
k15 := 2;
k16 := 7/5;
k17 := 2;
k18 := 6;
k19 := 1/5;
k20 := 1;
k21 := 5;
k22 := 1/25;
k23 := 2;
k24 := 10;
k25 := 2;
k26 := 30;
k27 := 2;
k28 := 3/50;
k29 := 1/5;
k30 := 1/2;
k31 := 1/2;
k32 := 1/200;
k33 := 3;
k34 := 9/2;
k35 := 2;
k36 := 20;
k37 := 9/10;
k38 := 9/125;
k39 := 1/2;
k40 := 9/20;
k41 := 1;
k42 := 133/10000;
k43 := 67/10000;
operator diff$
odes := { diff(x1, t) = (1*k41*k1 + (-1)*k41*k2*x1)/k41,
  diff(x2, t) = (1*k41*k3*x1 + (-1)*k41*k4*x2)/k41,
  diff(x3, t) = (1*k41*k5*x1 + (-1)*k41*k6*x3 + 1*k41*k5*x2)/k41,
  diff(x4, t) = (1*k41*k7*x3 + (-1)*k41*k8*x4 + 1*k41*k7*x5 + 1*k41*k7*x19 + 1*k41*k7*x13)/k41,
  diff(x5, t) = (1*k41*k9*x3 + (-1)*k41*k10*x5 + 1*k41*k9*x2)/k41,
  diff(x6, t) = (1*k41*k11*x4 + (-1)*k41*k13*x6 + 1*k41*k11*x10 + 1*k41*k11*x11 + 1*k41*k11*x15 + 1*k41*k11*x7 + 1*k41*k11*x17 + (-1)*k41*k13*x6*x16*k43)/k41,
  diff(x7, t) = (1*k41*k12*x4 + (-1)*k41*k14*x7 + 1*k41*k12*x15 + 1*k41*k12*x16)/k41,
  diff(x8, t) = (1*k41*k15*x5 + (-1)*k41*k16*x8)/k41,
  diff(x9, t) = (1*k41*k17*x8 + (-1)*k41*k18*x9 + 1*k41*k17*x16)/k41,
  diff(x10, t) = (1*k41*k19*x9 + (-1)*k41*k20*x10 + 1*k41*k19*x20)/k41,
  diff(x11, t) = (1*k41*k21 + (-1)*k41*k22*x11 + (-1)*k41*k22*x11*x18*k42)/k41,
  diff(x12, t) = (1*k41*k29*x6 + (-1)*k41*k30*x12)/k41,
  diff(x13, t) = (1*k41*k27*x19 + (-1)*k41*k28*x13*x12)/k41,
  diff(x14, t) = ((-1)*k41*k24*x14*x12 + 1*k41*k23*x18 + 1*k41*k23*x13)/k41,
  diff(x15, t) = (1*k41*k31 + (-1)*k41*k32*x15)/k41,
  diff(x16, t) = (1*k41*k33*x15 + (-1)*k41*k34*x16)/k41,
  diff(x17, t) = (1*k41*k35*x15 + (-1)*k41*k36*x17 + 1*k41*k35*x7)/k41,
  diff(x18, t) = (1*k41*k37*x2 + (-1)*k41*k38*x18)/k41,
  diff(x19, t) = (1*k41*k39*x2 + (-1)*k41*k40*x19 + 1*k41*k39*x14)/k41,
  diff(x20, t) = (1*k41*k25*x14 + 1*k41*k25*x13 + (-1)*k41*k26*x20)/k41 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file