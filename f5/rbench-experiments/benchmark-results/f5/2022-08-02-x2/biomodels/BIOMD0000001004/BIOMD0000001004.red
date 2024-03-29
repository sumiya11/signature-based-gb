load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 627/10000;
k3 := 5111/5000;
k4 := 37037/2500;
k5 := 22047/10000;
k6 := 423/5000;
k7 := 163647/10000;
k8 := 221/10000;
k9 := 27893/10000;
k10 := 193743/10000;
k11 := 386891/500;
k12 := 479/10000;
k13 := 0;
k14 := 367/5000;
k15 := 74237/10000;
k16 := 2023/2500;
k17 := 53/10000;
k18 := 2641/2000;
k19 := 4967/1000;
k20 := 1;
k21 := 1;
operator diff$
odes := { diff(x1, t) = (-1)*k1*k2*x1/k1,
  diff(x2, t) = 1*k1*k2*x1/k1,
  diff(x3, t) = (1*k1*k3 + 1*k1*k4*x4 + (-1)*k1*k5*x3)/k1,
  diff(x4, t) = (1*k1*k7*x5*x2 + (-1)*k1*k9*x4)/k1,
  diff(x5, t) = (1*k1*k6*x3 + (-1)*k1*k7*x5*x2 + (-1)*k1*k8*x5)/k1,
  diff(x6, t) = (-1)*k1*k10*x6/k1,
  diff(x7, t) = 1*k1*k10*x6/k1,
  diff(x8, t) = ((-1)*k1*k11*x8*x9 + (-1)*k1*k12*x8 + 1*k1*k19*x7*x4)/k1,
  diff(x9, t) = (1*k1*k17*x10 + (-1)*k1*k18*x9)/k1,
  diff(x10, t) = (1*k1*k13 + 1*k1*k14*x7 + (-1)*k1*k15*x10*x4 + (-1)*k1*k16*x10)/k1 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file