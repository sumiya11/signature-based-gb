load_package groebner$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 50;
k4 := 1000;
k5 := 1/20;
k6 := 57/2500;
k7 := 19/50000;
k8 := 1/62500;
k9 := 3/12500;
k10 := 19/50000;
k11 := 1/20;
k12 := 1/200;
k13 := 3/12500;
k14 := 1/200;
k15 := 3/12500;
k16 := 50;
k17 := 42;
k18 := 3/3125;
k19 := 9/6250;
k20 := 19/50000;
k21 := 57/2500;
k22 := 10;
operator diff$
odes := { diff(x1, t) = (1*k2*k4*k1/(k3 + k1) + (-1)*k2*k7*x1 + (-2)*k2*(k8*x1^2 - k9*x4) + (-1)*k2*(k12*x1*x3 - k13*x6))/k2,
  diff(x2, t) = (1*k2*k5*x6 + (-1)*k2*k6*x2 + (-1)*k2*(k18*x7*x2 - k19*x8))/k2,
  diff(x3, t) = (1*k2*k5*x6 + 1*k2*k11*x5 + (-1)*k2*(k12*x1*x3 - k13*x6) + (-1)*k2*(k14*x4*x3 - k15*x5))/k2,
  diff(x4, t) = (1*k2*(k8*x1^2 - k9*x4) + (-1)*k2*k10*x4 + (-1)*k2*(k14*x4*x3 - k15*x5))/k2,
  diff(x5, t) = ((-1)*k2*k11*x5 + 1*k2*(k14*x4*x3 - k15*x5))/k2,
  diff(x6, t) = ((-1)*k2*k5*x6 + 1*k2*(k12*x1*x3 - k13*x6))/k2,
  diff(x7, t) = (1*k2*k17*k1/(k16 + k1) + (-1)*k2*(k18*x7*x2 - k19*x8) + (-1)*k2*k21*x7)/k2,
  diff(x8, t) = (1*k2*(k18*x7*x2 - k19*x8) + (-1)*k2*k20*x8)/k2 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file