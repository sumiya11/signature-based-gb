load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 110;
k4 := 1500;
k5 := 1000;
k6 := 70;
k7 := 321/2;
k8 := 847;
k9 := 420;
k10 := 360;
k11 := 70;
k12 := 847;
k13 := 321/2;
k14 := 420;
k15 := 847;
k16 := 360;
k17 := 133/100;
k18 := 16;
k19 := 13/1000;
k20 := 90;
k21 := 330;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = (1*k2*k3*x3/(k4*(1 + x2/k5) + x3) + (-1)*k2*(k13/(1 + k21/k11)*x2/(k16*(1 + x3/k11) + x2) + k14/(1 + k21/k12)*x2/(k15*(1 + x3/k12) + x2)) + (-1)*k2*k19*x2/(k20 + x2))/k2,
  diff(x3, t) = ((-1)*k2*k3*x3/(k4*(1 + x2/k5) + x3) + 1*k1*(k21/(k6 + k21)*k7/(1 + x2/k10 + x3/k6) + k21/(k8 + k21)*k9/(1 + x2/k10 + x3/k8)) + (-1)*k2*k17*x3/(k18 + x3))/k2 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file