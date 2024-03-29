load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 59/100;
k2 := 84/125;
k3 := 339/500;
k4 := 49/1250;
k5 := 277/500;
k6 := 2;
k7 := 0;
k8 := 1;
operator diff$
odes := { diff(x1, t) = (-k1)*(k7 + x8)*x1 - k2*x7*x1,
  diff(x2, t) = k2*x7*x1 - k3*x7^k6*x2,
  diff(x3, t) = k3*x7^k6*x2 - 4*k4*x3,
  diff(x4, t) = 4*k4*x3 - 4*k4*x4,
  diff(x5, t) = 4*k4*x4 - 4*k4*x5,
  diff(x6, t) = 4*k4*x5 - 4*k4*x6,
  diff(x7, t) = 4*k4*x6 - k5*x7,
  diff(x8, t) = k5*x7,
  diff(x9, t) = k1*(k7 + x8)*x1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file