load_package f5$
torder({}, revgradlex)$

k1 := 1/100;
k2 := 2;
k3 := 1/250;
k4 := 33/100;
k5 := 1/250;
k6 := 2;
k7 := 50;
k8 := 2;
k9 := 2000;
k10 := 2;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (1*k11*k2 + (-1)*k11*k1*x1 + (-1)*k11*k3*x1*x2)/k11,
  diff(x2, t) = (1*k11*k7*x3 + (-1)*k11*k8*x2)/k11,
  diff(x3, t) = (1*k11*k3*x1*x2 + (-1)*k11*k4*x3 + (-1)*k11*k5*x4*x3)/k11,
  diff(x4, t) = (1*k11*k9*x5 + (-1)*k11*k10*x4)/k11,
  diff(x5, t) = (1*k11*k5*x4*x3 + (-1)*k11*k6*x5)/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file