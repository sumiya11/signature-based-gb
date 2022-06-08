load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 1/25;
k2 := 1/1000;
k3 := 5;
k4 := 1/10;
k5 := 3/100;
k6 := 100;
k7 := 1;
k8 := 1/100;
k9 := 1/100;
k10 := 3;
k11 := 1/100;
k12 := 1;
k13 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k12*k2*x1*x3 + 1*k12*k1*x2 + 1*k12*k6 + (-1)*k12*k8*x1 + (-1)*k12*k7*x1*x5)/k12,
  diff(x2, t) = (1*k12*k2*x1*x3 + (-1)*k12*k1*x2)/k12,
  diff(x3, t) = ((-1)*k12*k2*x1*x3 + 1*k12*k1*x2)/k12,
  diff(x4, t) = (1*k12*k4*x3 + 1*k12*k3*x2 + (-1)*k12*k5*x4)/k12,
  diff(x5, t) = (1*k12*k10*x4 + (-1)*k12*k9*x5 + (-1)*k12*k7*x1*x5)/k12,
  diff(x6, t) = (1*k12*k7*x1*x5 + (-1)*k12*k11*x6)/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file