load_package f5$
torder({}, revgradlex)$

k1 := 0;
k2 := 0;
k3 := 1;
k4 := 1/5;
k5 := 1/20;
k6 := 19/50000;
k7 := 19/50000;
k8 := 1/62500;
k9 := 3/12500;
k10 := 19/50000;
k11 := 1/20;
k12 := 1/200;
k13 := 3/12500;
k14 := 1/200;
k15 := 3/12500;
k16 := 5/2;
k17 := 10;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (1*k3*k4*k16 + (-1)*k3*k7*x2 + (-2)*k3*(k8*x2^2 - k9*x5) + (-1)*k3*(k12*x2*x4 - k13*x7))/k3,
  diff(x3, t) = (1*k3*k5*x7 + (-1)*k3*k6*x3)/k3,
  diff(x4, t) = (1*k3*k5*x7 + 1*k3*k11*x6 + (-1)*k3*(k12*x2*x4 - k13*x7) + (-1)*k3*(k14*x5*x4 - k15*x6))/k3,
  diff(x5, t) = (1*k3*(k8*x2^2 - k9*x5) + (-1)*k3*k10*x5 + (-1)*k3*(k14*x5*x4 - k15*x6))/k3,
  diff(x6, t) = ((-1)*k3*k11*x6 + 1*k3*(k14*x5*x4 - k15*x6))/k3,
  diff(x7, t) = ((-1)*k3*k5*x7 + 1*k3*(k12*x2*x4 - k13*x7))/k3 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file