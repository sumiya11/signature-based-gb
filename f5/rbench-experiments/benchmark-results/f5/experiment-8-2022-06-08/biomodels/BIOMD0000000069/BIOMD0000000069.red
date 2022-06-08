load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 4/5;
k3 := 1;
k4 := 10;
k5 := 1;
k6 := 1;
k7 := 1/20;
k8 := 3/20;
k9 := 7/200;
k10 := 1/10000;
k11 := 0;
k12 := 1/10;
k13 := 1/100;
k14 := 1;
k15 := 0;
k16 := 1;
k17 := 1;
k18 := 0;
k19 := 0;
k20 := 1/10000;
k21 := 1;
k22 := 1;
k23 := 1;
k24 := 1;
k25 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*(k2*(k11 + k17*x8)*x1 - k1*x5*x2)*k21 + 1*k21*k4*k7*x4)/k21,
  diff(x2, t) = (1*(k2*(k11 + k17*x8)*x1 - k1*x5*x2)*k21 + (-1)*(k3*(k15*x2 + k14*x3 + k10 + k16*x4)*x2 - k7*x3)*k21)/k21,
  diff(x3, t) = (1*(k3*(k15*x2 + k14*x3 + k10 + k16*x4)*x2 - k7*x3)*k21 + (-1)*(k1*x5*x3 - k2*(k11 + k17*x8)*x4)*k21)/k21,
  diff(x4, t) = (1*(k1*x5*x3 - k2*(k11 + k17*x8)*x4)*k21 + (-1)*k21*k4*k7*x4)/k21,
  diff(x5, t) = 1*(x10*k12*x6 - k13*x5)*k21/k21,
  diff(x6, t) = (-1)*(x10*k12*x6 - k13*x5)*k21/k21,
  diff(x7, t) = (-1)*k21*((k5*(k15*x2 + k14*x3 + k10 + k16*x4) + k9)*x7 - k8*x8)/k21,
  diff(x8, t) = 1*k21*((k5*(k15*x2 + k14*x3 + k10 + k16*x4) + k9)*x7 - k8*x8)/k21,
  diff(x9, t) = (-1)*k6*(k15*x2 + k14*x3 + k10 + k16*x4)*x9*k21/k21,
  diff(x10, t) = ((-1)*(x10*k12*x6 - k13*x5)*k21 + 1*k6*(k15*x2 + k14*x3 + k10 + k16*x4)*x9*k21)/k21 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file