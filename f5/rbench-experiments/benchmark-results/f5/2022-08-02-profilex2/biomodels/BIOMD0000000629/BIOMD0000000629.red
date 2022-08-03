load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 3/5;
k3 := 1/10;
k4 := 7/500;
k5 := 1/5;
k6 := 1/2000;
k7 := 7/2000;
k8 := 30;
operator diff$
odes := { diff(x1, t) = (-1)*k1*(k2*x1*x3 - k3*x2)/k1,
  diff(x2, t) = (1*k1*(k2*x1*x3 - k3*x2) + (-1)*k1*(k4*x2*x4 - k5*x5))/k1,
  diff(x3, t) = (-1)*k1*(k2*x1*x3 - k3*x2)/k1,
  diff(x4, t) = (-1)*k1*(k4*x2*x4 - k5*x5)/k1,
  diff(x5, t) = 1*k1*(k4*x2*x4 - k5*x5)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file