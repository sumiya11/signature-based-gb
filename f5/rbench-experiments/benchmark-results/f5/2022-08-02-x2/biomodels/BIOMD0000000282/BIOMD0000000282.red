load_package f5$
torder({}, revgradlex)$

k1 := 11;
k2 := 0;
k3 := 83/5;
k4 := 18/25;
k5 := 1;
k6 := 0;
k7 := 34/25;
k8 := 2;
operator diff$
odes := { diff(x1, t) = ((-1)*k5*(k1*x1*x2 - k2*x3) + 1*k5*k3*x3*x2 + 1*k5*k4*x3*k6)/k5,
  diff(x2, t) = ((-1)*k5*(k1*x1*x2 - k2*x3) + (-1)*k5*k3*x3*x2)/k5,
  diff(x3, t) = (1*k5*(k1*x1*x2 - k2*x3) + (-1)*k5*k3*x3*x2 + (-1)*k5*k4*x3*k6)/k5,
  diff(x4, t) = 1*k5*k3*x3*x2/k5,
  diff(x5, t) = 0,
  diff(x6, t) = 1*k5*k4*x3*k6/k5 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file