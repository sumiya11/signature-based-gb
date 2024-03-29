load_package f5$
1;
torder({}, revgradlex)$

k1 := 1;
k2 := 3/10;
k3 := 1;
k4 := 21/10;
k5 := 1/100;
k6 := 6/5;
k7 := 3/50;
k8 := 20;
k9 := 1;
operator diff$
odes := { diff(x1, t) = (1*k9*((1 + k3*(1 - k4*k2) + 1/2*k4*k6*k2)*x1 + k4*(1 + k3*k2)*x3) + (-1)*k9*(x1*x1 + x1*x2))/k9,
  diff(x2, t) = (1*k9*x2*x3 + (-1)*k9*(k7 - k5)*x2)/k9,
  diff(x3, t) = (1*k9*(k3*(1 - k4*k2) + k4*k6*k2)*x1 + (-1)*k9*(k8 - (k3*(1 + k4*k2) + k4*(1 + 1/2*k6*k2)))*x3)/k9 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file