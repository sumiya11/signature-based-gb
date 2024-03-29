load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 60000;
k2 := 100;
k3 := 20000000000;
k4 := 100;
k5 := 20000000000;
k6 := 20000;
k7 := 0;
k8 := 1;
k9 := 1;
k10 := 1/10000000000;
k11 := 1/10000000000;
operator diff$
odes := { diff(x1, t) = ((-1)*k9*(k5*x2*x1 - k6*x3)/k9 + (-1)*k9*(k3*x2*x1 - k4*x4)/k9 + (-1)*k9*(k5*x2*x1 - k6*x5)/k9 + (-1)*k9*(k3*x3*x1 - k4*x6)/k9 + (-1)*k9*(k5*x3*x1 - k6*x7)/k9 + (-1)*k9*(k5*x4*x1 - k6*x6)/k9 + (-1)*k9*(k5*x4*x1 - k6*x8)/k9 + (-1)*k9*(k5*x5*x1 - k6*x7)/k9 + (-1)*k9*(k3*x5*x1 - k4*x8)/k9 + (-1)*k9*(k5*x6*x1 - k6*x9)/k9 + (-1)*k9*(k3*x7*x1 - k4*x9)/k9 + (-1)*k9*(k5*x8*x1 - k6*x9)/k9)/k9,
  diff(x2, t) = ((-1)*k9*(k5*x2*x1 - k6*x3)/k9 + (-1)*k9*(k3*x2*x1 - k4*x4)/k9 + (-1)*k9*(k5*x2*x1 - k6*x5)/k9)/k9,
  diff(x3, t) = (1*k9*(k5*x2*x1 - k6*x3)/k9 + (-1)*k9*(k3*x3*x1 - k4*x6)/k9 + (-1)*k9*(k5*x3*x1 - k6*x7)/k9)/k9,
  diff(x4, t) = (1*k9*(k3*x2*x1 - k4*x4)/k9 + (-1)*k9*(k5*x4*x1 - k6*x6)/k9 + (-1)*k9*(k5*x4*x1 - k6*x8)/k9)/k9,
  diff(x5, t) = (1*k9*(k5*x2*x1 - k6*x5)/k9 + (-1)*k9*(k5*x5*x1 - k6*x7)/k9 + (-1)*k9*(k3*x5*x1 - k4*x8)/k9)/k9,
  diff(x6, t) = (1*k9*(k3*x3*x1 - k4*x6)/k9 + 1*k9*(k5*x4*x1 - k6*x6)/k9 + (-1)*k9*(k5*x6*x1 - k6*x9)/k9 + (-1)*k9*(k1*x6 - k2*x10)/k9)/k9,
  diff(x7, t) = (1*k9*(k5*x3*x1 - k6*x7)/k9 + 1*k9*(k5*x5*x1 - k6*x7)/k9 + (-1)*k9*(k3*x7*x1 - k4*x9)/k9)/k9,
  diff(x8, t) = (1*k9*(k5*x4*x1 - k6*x8)/k9 + 1*k9*(k3*x5*x1 - k4*x8)/k9 + (-1)*k9*(k5*x8*x1 - k6*x9)/k9 + (-1)*k9*(k1*x8 - k2*x11)/k9)/k9,
  diff(x9, t) = (1*k9*(k5*x6*x1 - k6*x9)/k9 + 1*k9*(k3*x7*x1 - k4*x9)/k9 + 1*k9*(k5*x8*x1 - k6*x9)/k9 + (-1)*k9*(k1*x9 - k2*x13)/k9 + (-1)*k9*(k1*x9 - k2*x12)/k9)/k9,
  diff(x10, t) = 1*k9*(k1*x6 - k2*x10)/k9/k9,
  diff(x11, t) = 1*k9*(k1*x8 - k2*x11)/k9/k9,
  diff(x12, t) = (1*k9*(k1*x9 - k2*x12)/k9 + (-1)*k9*(k1*x12 - k2*x14)/k9)/k9,
  diff(x13, t) = (1*k9*(k1*x9 - k2*x13)/k9 + (-1)*k9*(k1*x13 - k2*x14)/k9)/k9,
  diff(x14, t) = (1*k9*(k1*x13 - k2*x14)/k9 + 1*k9*(k1*x12 - k2*x14)/k9)/k9 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file