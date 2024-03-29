load_package f5$
torder({}, revgradlex)$

k1 := 1/10;
k2 := 109/500;
k3 := 1;
k4 := 292408052354609/100000000000000;
k5 := 499999999999999/5000000000000000;
k6 := 263/1000;
k7 := 547/1000;
k8 := 1;
k9 := 292408052354609/10000000000000;
k10 := 239254806051979/1000000000000000;
k11 := 183/100;
k12 := 1460500000017999/20000000000000;
k13 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*x1 + 1*k5*x1)/k13,
  diff(x2, t) = (1*(k2 + k3*x1/(x1 + k4))*x1 + (-1)*k6*x2 + 1*k10*x2)/k13,
  diff(x3, t) = (1*(k7 + k8*x2/(x2 + k9))*x2 + (-1)*k11*x3)/k13 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file