load_package f5$
torder({}, revgradlex)$

k1 := 1/200;
k2 := 19/2000;
k3 := 1/10000;
k4 := 1/10000;
k5 := 1/20;
k6 := 1;
operator diff$
odes := { diff(x1, t) = (1*k6*k5 + (-1)*k6*(k2*x1*x2 + k3*x1))/k6,
  diff(x2, t) = (1*k6*(k2*x1*x2 + k4*x3) + (-1)*k6*k1*x1*x2)/k6,
  diff(x3, t) = (1*k6*(k1*x1*x2 + k3*x1) + (-1)*k6*k4*x3)/k6 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file