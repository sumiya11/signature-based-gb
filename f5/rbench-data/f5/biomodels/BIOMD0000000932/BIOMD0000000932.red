load_package f5$
torder({}, revgradlex)$

k1 := 1713/5000;
k2 := 53/10;
k3 := 4;
k4 := 2;
k5 := 23/10;
k6 := 30;
k7 := 1/10;
k8 := 1;
operator diff$
odes := { diff(x1, t) = (1*k8*(k1*k6*x1 - k2*x3*x1) + (-1)*k8*k4*x1)/k8,
  diff(x2, t) = (1*k8*k4*x1 + (-1)*k8*k5*x2)/k8,
  diff(x3, t) = (1*k8*k5*x2 + (-1)*k8*k3*x3)/k8,
  diff(x4, t) = 1*k8*k7*x3*x1*x4/k8 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file