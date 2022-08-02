load_package f5$
torder({}, revgradlex)$

k1 := 1/200;
k2 := 50000;
k3 := 9/500000000;
k4 := 17123/50000000;
k5 := 1/2500;
k6 := 30;
k7 := 3/500000000;
k8 := 391/10000000;
k9 := 1/10;
k10 := 689/2000;
k11 := 50000;
k12 := 800;
k13 := 1;
k14 := 1;
operator diff$
odes := { diff(x1, t) = (1*k14*k1*x1*(1 - x1/k11)*(x1/k12 - 1) + (-1)*k14*k3*x1*x2)/k14,
  diff(x2, t) = (1*k14*k3*x1*x2 + (-1)*k14*k4*x2 + (-1)*k14*k5*x2)/k14,
  diff(x3, t) = (1*k13*k6 + (-1)*k7*x2*x3 + (-1)*k13*k8*x3)/k13,
  diff(x4, t) = (1*k7*x2*x3 + (-1)*k13*k8*x4 + (-1)*k13*k10*x4 + (-1)*k13*k9*x4)/k13,
  diff(x5, t) = (1*k13*k9*x4 + (-1)*k13*k8*x5)/k13 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file