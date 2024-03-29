load_package f5$
1;
torder({}, revgradlex)$

k1 := 1/2000000000;
k2 := 1/2000000000;
k3 := 1;
k4 := 4;
k5 := 2;
k6 := 2400;
k7 := 9/10;
k8 := 14;
k9 := 1;
k10 := 400000000;
k11 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k11*k1*x5*x1 + (-1)*k11*(k2*x6*x1 - k3*x4))/k11,
  diff(x2, t) = (1*k11*k1*x5*x1 + (-1)*k11*(k5*x2 + k2*x6*x2) + (-1)*k11*k4*x2)/k11,
  diff(x3, t) = (1*k11*(k5*x2 + k2*x6*x2) + (-1)*k11*k4*x3)/k11,
  diff(x4, t) = 1*k11*(k2*x6*x1 - k3*x4)/k11,
  diff(x5, t) = (1*k11*(k6*x2 + (1 - k7)*k6*x3) + (-1)*k11*k8*x5)/k11,
  diff(x6, t) = (1*k11*k9*k6*(x2 + x3) + (-1)*k11*k8*x6)/k11 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file