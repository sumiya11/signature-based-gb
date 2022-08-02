load_package f5$
torder({}, revgradlex)$

k1 := 175/11;
k2 := 1/10;
k3 := 1/5;
k4 := 1;
k5 := 1;
k6 := 10;
k7 := 1/10;
k8 := 1/10;
k9 := 1/2;
k10 := 1/5;
k11 := 1/10;
k12 := 1/100;
k13 := 1/1000000;
k14 := 31/2;
k15 := 1;
k16 := 1/10;
k17 := 1;
k18 := 12;
k19 := 1/10;
k20 := 1/100;
k21 := 1;
operator diff$
odes := { diff(x1, t) = (1*k21*k6 + (-1)*k21*k7*x1 + (-1)*k21*k8*x5*x1 + (-1)*k21*k9*x5*x1)/k21,
  diff(x2, t) = (1*k21*k8*x5*x1 + (-1)*k21*k2*x2 + (-1)*k21*k12*x2 + 1*k21*k11*x4 + (-1)*k21*k13*x6*x2 + (-1)*k21*k17*x16*x2)/k21,
  diff(x3, t) = (1*k21*k12*x2 + (-1)*k21*k3*x3 + (-1)*k21*k13*x6*x3 + (-1)*k21*k17*x16*x3)/k21,
  diff(x4, t) = (1*k21*k9*x5*x1 + (-1)*k21*k11*x4 + (-1)*k21*k7*x4)/k21,
  diff(x5, t) = (1*k21*k4*x3 + (-1)*k21*k5*x5)/k21,
  diff(x6, t) = (1*k21*k10*x15 + 1*k21*k14*(x2 + x3)*x6 + (-1)*k21*k16*x6)/k21,
  diff(x7, t) = (-1)*k21*k15*x7/k21,
  diff(x8, t) = (2*k21*k15*x7 + (-1)*k21*k15*x8)/k21,
  diff(x9, t) = (2*k21*k15*x8 + (-1)*k21*k15*x9)/k21,
  diff(x10, t) = (2*k21*k15*x9 + (-1)*k21*k15*x10)/k21,
  diff(x11, t) = (2*k21*k15*x10 + (-1)*k21*k15*x11)/k21,
  diff(x12, t) = (2*k21*k15*x11 + (-1)*k21*k15*x12)/k21,
  diff(x13, t) = (2*k21*k15*x12 + (-1)*k21*k15*x13)/k21,
  diff(x14, t) = (2*k21*k15*x13 + (-1)*k21*k15*x14)/k21,
  diff(x15, t) = (2*k21*k15*x14 + (-1)*k21*k15*x15)/k21,
  diff(x16, t) = (1*k21*k20 + 1*k21*k18*(x2 + x3)*x16 + (-1)*k21*k19*x16)/k21 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file