load_package f5$
torder({}, revgradlex)$

k1 := 18/25;
k2 := 36/25;
k3 := 864;
k4 := 216/5;
k5 := 20000;
k6 := 432;
k7 := 3/50;
k8 := 21/25000;
k9 := 3/1250000;
k10 := 1;
operator diff$
odes := { diff(x1, t) = k3 - (k2 + k1*x2)*x1,
  diff(x2, t) = x3*k4*x1^2/(k5 + x1^2) - k6*x2,
  diff(x3, t) = ((-k7) + k8*x1 - k9*x1^2)*x3 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file