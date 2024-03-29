load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 127;
k2 := 31/2000;
k3 := 231/2000;
k4 := 6879/2500;
k5 := 2659/1250;
k6 := 9/20;
k7 := 75669/50000;
k8 := 329/100;
k9 := 2587/2000;
k10 := 12607/10000;
k11 := 4712/5283;
k12 := 565/10566;
k13 := 2/1761;
k14 := 0;
k15 := 5283/5000;
k16 := 6879/2500;
k17 := 9/20;
k18 := 1;
k19 := 10001/10000;
operator diff$
odes := { diff(x1, t) = ((-1)*k18*k4*x4*x1 + (-1)*k18*k2*(1 - x3/k3))/k18,
  diff(x2, t) = (1*k18*k2*(1 - x3/k3) + (-1)*k18*k16*(1 - k17)*x4*x2)/k18,
  diff(x3, t) = 1*k18*k2*(1 - x3/k3)/k18,
  diff(x4, t) = (1*k18*k4*x4*x1 + 1*k18*k16*(1 - k17)*x4*x2 + (-1)*k18*k5*x4)/k18,
  diff(x5, t) = 1*k18*k5*x4/k18 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file