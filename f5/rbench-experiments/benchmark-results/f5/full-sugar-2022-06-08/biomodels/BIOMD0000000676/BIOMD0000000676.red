load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 1/10;
k2 := 1/10;
k3 := 189/100;
k4 := 57/5;
k5 := 0;
k6 := 91/100;
k7 := 539/10;
k8 := 33/1000;
k9 := 1/10;
k10 := 1/10;
k11 := 91/100;
k12 := 57/5;
k13 := 189/100;
k14 := 129/50;
k15 := 98;
k16 := 63/5;
k17 := 91/100;
k18 := 333/100;
k19 := 899/10;
k20 := 147/5;
k21 := 91/100;
k22 := 1;
k23 := 100;
k24 := 150;
k25 := 3/200;
operator diff$
odes := { diff(x1, t) = 0,
  diff(x2, t) = ((-1)*k22*(k1*k23*x2 - k2*x3) + (-1)*k22*k6*x2 + 1*k22*k7*x6 + 1*k22*k8*k24*x8 + 1*k22*(k9*x9 - k10*x2*x10))/k22,
  diff(x3, t) = (1*k22*(k1*k23*x2 - k2*x3) + (-1)*k22*k11*x3)/k22,
  diff(x4, t) = ((-1)*k22*(k3*k23*x4 - k4*x5) + 1*k22*k6*x2 + 1*k22*(k12*x11 - k13*x4*x10))/k22,
  diff(x5, t) = (1*k22*(k3*k23*x4 - k4*x5) + 1*k22*k11*x3 + (-1)*k22*k14*k24*x5 + 1*k22*k15*x13)/k22,
  diff(x6, t) = ((-1)*k22*k7*x6 + 1*k22*k20*x14 + (-1)*k22*k21*x6)/k22,
  diff(x7, t) = 1*k22*k7*x6/k22,
  diff(x8, t) = ((-1)*k22*k8*k24*x8 + 1*k22*k21*x6)/k22,
  diff(x9, t) = ((-1)*k22*(k9*x9 - k10*x2*x10) + 1*k22*k16*x13 + (-1)*k22*k17*x9)/k22,
  diff(x10, t) = (1*k22*(k9*x9 - k10*x2*x10) + 1*k22*(k12*x11 - k13*x4*x10))/k22,
  diff(x11, t) = ((-1)*k22*(k12*x11 - k13*x4*x10) + 1*k22*k17*x9 + (-1)*k22*k18*k24*x11 + 1*k22*k19*x14)/k22,
  diff(x12, t) = 0,
  diff(x13, t) = (1*k22*k14*k24*x5 + (-1)*k22*k15*x13 + (-1)*k22*k16*x13)/k22,
  diff(x14, t) = (1*k22*k18*k24*x11 + (-1)*k22*k20*x14 + (-1)*k22*k19*x14)/k22,
  diff(x15, t) = 1*k22*k20*x14/k22 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file