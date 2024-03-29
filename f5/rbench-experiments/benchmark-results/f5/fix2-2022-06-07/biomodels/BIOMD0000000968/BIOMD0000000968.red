load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 1;
k2 := 1/10;
k3 := 1/10;
k4 := 1;
k5 := 1/10;
k6 := 1;
k7 := 1/10;
k8 := 1;
k9 := 1/10;
k10 := 1/10;
k11 := 0;
k12 := 1/1000000000000;
k13 := 1;
k14 := 9963234242563/4000000000000;
k15 := 41513476010679/25000000000000;
operator diff$
odes := { diff(x1, t) = (-1)*k12*(k1*x1*x7 - k2*x6)/k12,
  diff(x2, t) = 0/k12,
  diff(x3, t) = ((-1)*k12*(k4*x8*x3 - k5*x9) + 1*k12*k10*x4)/k12,
  diff(x4, t) = ((-1)*k12*k10*x4 + 1*k12*k6*x9)/k12,
  diff(x5, t) = (1*k12*k7*x4 + (-1)*k12*k8*x8*x5)/k12,
  diff(x6, t) = (1*k12*(k1*x1*x7 - k2*x6) + (-1)*k12*(k3*x6 - k9*x8))/k12,
  diff(x7, t) = (-1)*k12*(k1*x1*x7 - k2*x6)/k12,
  diff(x8, t) = (1*k12*(k3*x6 - k9*x8) + (-1)*k12*(k4*x8*x3 - k5*x9) + (-1)*k12*k8*x8*x5 + 1*k12*k6*x9)/k12,
  diff(x9, t) = (1*k12*(k4*x8*x3 - k5*x9) + (-1)*k12*k6*x9)/k12,
  diff(x10, t) = 1*k12*k8*x8*x5/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file