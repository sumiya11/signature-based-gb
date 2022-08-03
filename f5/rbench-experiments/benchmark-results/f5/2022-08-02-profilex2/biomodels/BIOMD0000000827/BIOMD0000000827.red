load_package f5$
torder({}, revgradlex)$

k1 := 31/25;
k2 := 23000000000;
k3 := 27/25;
k4 := 1900000000000;
k5 := 1;
k6 := 1;
k7 := 240000000000;
k8 := 2600000000000;
k9 := 9/200;
k10 := 10000000000;
k11 := 10000000000;
k12 := 7/250;
k13 := 7/250;
k14 := 7/250;
k15 := 1;
k16 := 9/2500000000000;
k17 := 11/200000000000000;
k18 := 37/5000000000000;
operator diff$
odes := { diff(x1, t) = ((-2)*k15*(1/2*k2*x1^2 - k1*x2) + (-1)*k15*(k7*x1*x3 - k5*x5))/k15,
  diff(x2, t) = 1*k15*(1/2*k2*x1^2 - k1*x2)/k15,
  diff(x3, t) = ((-2)*k15*(1/2*k4*x3^2 - k3*x4) + (-1)*k15*(k7*x1*x3 - k5*x5))/k15,
  diff(x4, t) = (1*k15*(1/2*k4*x3^2 - k3*x4) + (-1)*k15*k11*x4*x8 + 1*k15*k14*x10)/k15,
  diff(x5, t) = (1*k15*(k7*x1*x3 - k5*x5) + (-1)*k15*k10*x5*x8 + 1*k15*k13*x9)/k15,
  diff(x6, t) = (-2)*k15*(1/2*k8*x6^2 - k6*x7)/k15,
  diff(x7, t) = (1*k15*(1/2*k8*x6^2 - k6*x7) + (-1)*k15*(k9*x7 - k12*x8))/k15,
  diff(x8, t) = 1*k15*(k9*x7 - k12*x8)/k15,
  diff(x9, t) = (1*k15*k10*x5*x8 + (-1)*k15*k13*x9)/k15,
  diff(x10, t) = (1*k15*k11*x4*x8 + (-1)*k15*k14*x10)/k15 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file