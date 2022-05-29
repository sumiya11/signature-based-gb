load_package groebner$
torder({}, revgradlex)$

k1 := 9/25;
k2 := 1/2;
k3 := 9/25;
k4 := 639/5000;
k5 := 12/25;
k6 := 4/25;
k7 := 1/5;
k8 := 3/5;
k9 := 29/100;
k10 := 4/25;
k11 := 2;
k12 := 1;
operator diff$
odes := { diff(x1, t) = (1*k12*k1*x1*(1 - (x1 + x2)) + (-1)*k12*k2*x1*x3 + (-1)*k12*(k3*x1*x4 + k4*x1))/k12,
  diff(x2, t) = (1*k12*k2*x1*x3 + (-1)*k12*(k5*x2*x4 + x2))/k12,
  diff(x3, t) = ((-1)*k12*k2*x1*x3 + 1*k12*k11*x2 + (-1)*k12*(k6*x3*x4 + k7*x3))/k12,
  diff(x4, t) = (1*k12*(k8*x2*x4 + k9*x1*x4) + (-1)*k12*k10*x4)/k12 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file