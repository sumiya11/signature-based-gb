load_package f5$
torder({}, revgradlex)$

k1 := 1/10;
k2 := 1/10;
k3 := 1/5;
k4 := 1/100;
k5 := 1/100;
k6 := 1/20;
k7 := 1/20;
k8 := 1/100;
k9 := 1/100;
k10 := 10;
k11 := 20;
k12 := 5;
k13 := 1;
k14 := 1/10;
k15 := 1/100;
k16 := 1/100;
k17 := 1/100;
k18 := 9/10;
k19 := 1;
k20 := 1023/1000;
k21 := 0;
k22 := 0;
k23 := 10;
k24 := 10;
k25 := 1;
k26 := 1/4;
k27 := 1;
k28 := 1;
operator diff$
odes := { diff(x1, t) = (1*k28*k1*x1*(1 - x1/k10) + (-1)*k28*k5*x1 + (-1)*k28*k16*x3*x1)/k28,
  diff(x2, t) = (1*k28*k2*x2*(1 - x2/k11) + (-1)*k28*k6*x2 + (-1)*k28*k15*x3*x2)/k28,
  diff(x3, t) = (1*k28*k3*x5*x3*(1 - x3/k12) + (-1)*k28*k7*x3 + (-1)*k28*k17*x1*x3)/k28,
  diff(x4, t) = (1*k28*k15*x3*x2 + (-1)*k28*k8*x4 + (-1)*k28*k14*x4)/k28,
  diff(x5, t) = (1*k28*k14*x4 + 1*k28*k4*x3*x5*(1 - x5/k13) + (-1)*k28*k9*x5)/k28 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file