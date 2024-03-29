load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 10;
k3 := 1/10;
k4 := 1/5;
k5 := 1/5;
k6 := 1;
k7 := 1/2;
k8 := 1;
operator diff$
odes := { diff(x1, t) = (1*k8*k1*x3*x1 + (-1)*k8*k1*x1*x3*(x1 + x2)/k2 + (-1)*k8*k3*x1 + (-1)*k8*k4*x3*x1)/k8,
  diff(x2, t) = (1*k8*k4*x3*x1 + 1*k8*k1*x3*x2 + (-1)*k8*k1*x2*x3*(x1 + x2)/k2 + (-1)*k8*k5*x2)/k8,
  diff(x3, t) = (1*k8*k6*x2 + (-1)*k8*k7*x3)/k8 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file