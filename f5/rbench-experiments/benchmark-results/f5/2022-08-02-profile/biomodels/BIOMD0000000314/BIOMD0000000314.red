load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 174087/100000000;
k3 := 907541/100000;
k4 := 81033/250000;
k5 := 208769/500000;
k6 := 129843/500000;
k7 := 39243/10000000;
k8 := 300019/1000000;
k9 := 981611/10000000;
k10 := 426767/100000000;
k11 := 116389/10000000;
k12 := 14491/1250000;
k13 := 100;
k14 := 24;
k15 := 209;
operator diff$
odes := { diff(x1, t) = ((-1)*k2*x12*x1*k13 + (-1)*k6*x1*k13 + 1*k7*x2*k13)/k13,
  diff(x2, t) = (1*k6*x1*k13 + (-1)*k7*x2*k13)/k13,
  diff(x3, t) = (1*k2*x12*x1*k13 + (-1)*k3*x3*x7*k13)/k13,
  diff(x4, t) = (1*k3*x3*x7*k13 + (-1)*k4*x4*k13)/k13,
  diff(x5, t) = (1*k4*x4*k13 + (-1)*k5*x5*k13)/k13,
  diff(x6, t) = ((-1)*k8*x6*x3*k13 + (-1)*k8*x6*x4*k13 + 1*k9*x7*x8*k13)/k13,
  diff(x7, t) = (1*k8*x6*x3*k13 + 1*k8*x6*x4*k13 + (-1)*k9*x7*x8*k13)/k13,
  diff(x8, t) = 0/k13,
  diff(x9, t) = ((-1)*k10*x9*x7*k13 + 1*k11*x10*x8*k13)/k13,
  diff(x10, t) = (1*k10*x9*x7*k13 + (-1)*k11*x10*x8*k13)/k13,
  diff(x11, t) = 1*x10*k12*k13/k13 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file