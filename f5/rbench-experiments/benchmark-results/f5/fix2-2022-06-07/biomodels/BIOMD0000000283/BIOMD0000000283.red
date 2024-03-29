load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 0;
k2 := 1/2;
k3 := 1;
k4 := 8;
k5 := 1;
operator diff$
odes := { diff(x1, t) = (-1)*k3*(x2*x1 - k1*x3)/k3,
  diff(x2, t) = ((-1)*k3*(x2*x1 - k1*x3) + 1*k3*k2*x3)/k3,
  diff(x3, t) = (1*k3*(x2*x1 - k1*x3) + (-1)*k3*k2*x3)/k3,
  diff(x4, t) = 1*k3*k2*x3/k3 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file