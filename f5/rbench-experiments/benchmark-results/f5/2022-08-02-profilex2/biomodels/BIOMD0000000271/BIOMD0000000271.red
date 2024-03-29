load_package f5$
torder({}, revgradlex)$

k1 := 164683/5000000;
k2 := 516;
k3 := 41/390625;
k4 := 34427/2000000;
k5 := 748267/10000000;
k6 := 198761/20000000;
k7 := 317871/100000000;
k8 := 82021/5000000;
k9 := 203019/100;
k10 := 0;
k11 := 1;
k12 := 1;
k13 := 1;
k14 := 203019/100;
k15 := 0;
k16 := 0;
operator diff$
odes := { diff(x1, t) = (1*k1*k2*k13 + (-1)*k1*x1*k13 + (-1)*k3*x2*x1*k13 + 1*k4*x3*k13 + 1*k6*x4*k13)/k12,
  diff(x2, t) = ((-1)*k3*x2*x1*k13 + 1*k4*x3*k13 + 1*k6*x4*k13)/k11,
  diff(x3, t) = (1*k3*x2*x1*k13 + (-1)*k4*x3*k13 + (-1)*k5*x3*k13)/k12,
  diff(x4, t) = (1*k5*x3*k13 + (-1)*k6*x4*k13 + (-1)*k7*x4*k13 + (-1)*k8*x4*k13)/k13,
  diff(x5, t) = 1*k7*x4*k13/k13,
  diff(x6, t) = 1*k8*x4*k13/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file