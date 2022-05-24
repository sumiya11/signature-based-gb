load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 67/50;
k4 := 1600000000;
k5 := 8000;
k6 := 40000000;
k7 := 1;
k8 := 3/50;
k9 := 0;
operator diff$
odes := { diff(x1, t) = ((-1)*x1*k8*k3*k2 + (-1)*x1*x4*k4*k2 + 1*x3*k7*k2)/k2,
  diff(x2, t) = 0,
  diff(x3, t) = (1*k8*x4*k5*k2 + (-1)*x3*k7*k2)/k2,
  diff(x4, t) = (1*x1*k8*k3*k2 + (-1)*x1*x4*k4*k2 + 1*k8*x4*k5*k2 + (-2)*x4^2*k6*k2)/k2,
  diff(x5, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file