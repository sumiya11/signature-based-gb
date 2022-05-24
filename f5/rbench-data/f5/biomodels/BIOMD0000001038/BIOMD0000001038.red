load_package f5$
torder({}, revgradlex)$

k1 := 539/1250;
k2 := 299/100000000;
k3 := 9817/10000;
k4 := 2213/5000;
k5 := 2/5;
k6 := 2291/10000;
k7 := 561/625;
k8 := 9611/10000;
k9 := 443/2000;
k10 := 199/400;
k11 := 1;
operator diff$
odes := { diff(x1, t) = (1*k11*(k1*x1*(1 - k2*x1) + k9*x1*x3) + (-1)*k11*k3*x1*x2)/k11,
  diff(x2, t) = (1*k11*(k4*x2*(1 - k5*x2) + k6*x1*x2) + (-1)*k11*k10*x2*x3)/k11,
  diff(x3, t) = (1*k11*k7 + (-1)*k11*k8*x3)/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file