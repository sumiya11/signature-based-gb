load_package groebner$
torder({}, revgradlex)$

k1 := 1000;
k2 := 1/10;
k3 := 3/100;
k4 := 7/5;
k5 := 7200;
k6 := 5000;
k7 := 3/5;
k8 := 1/5;
k9 := 11;
k10 := 1;
operator diff$
odes := { diff(x1, t) = k1 - k6*x1*x3 - k2*x1 + (k5 + k8)*x4,
  diff(x2, t) = k3*x1^2 - k7*x2,
  diff(x3, t) = k4*x2 - k6*x1*x3 + (k5 + k9)*x4 - k8*x3,
  diff(x4, t) = k6*x1*x3 - (k5 + k9)*x4 - k8*x4 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file