load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 1;
k2 := 1/10;
k3 := 1/10;
k4 := 1/5;
k5 := 1;
k6 := 1/100;
k7 := 3/2;
k8 := 1/100;
k9 := 6;
k10 := 1/100;
k11 := 5/2;
k12 := 1/100;
operator diff$
odes := { diff(x1, t) = (1*k1*k2*x2 + (-1)*k3*x1/(k4 + x1))/k1,
  diff(x2, t) = (1*k5*(1 - x2)/(k6 + 1 - x2) + (-1)*k7*x3*x2/(k8 + x2))/k1,
  diff(x3, t) = (1*x1*k9*(1 - x3)/(k10 + 1 - x3) + (-1)*k11*x3/(k12 + x3))/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file