load_package groebner$
torder({}, revgradlex)$

k1 := 9/50;
k2 := 5000000;
k3 := 1101/10000000000;
k4 := 103/2500;
k5 := 1711/5000000000000;
k6 := 49/2000;
k7 := 10000000;
k8 := 1/20;
k9 := 0;
k10 := 1;
operator diff$
odes := { diff(x1, t) = (1*k10*k1*x1*(1 - x1/k2) + (-1)*k10*k3*x1*x2 + (-1)*k10*k8*x1)/k10,
  diff(x2, t) = (1*k10*k9*x2*x3 + (-1)*k10*k4*x2 + (-1)*k10*k5*x1*x2)/k10,
  diff(x3, t) = (1*k10*k6*x3*(1 - x3/k7) + (-1)*k10*k9*x2*x3)/k10 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file