load_package f5$
torder({}, revgradlex)$

k1 := 307000;
k2 := 1940;
k3 := 5490000;
k4 := 209/10000;
k5 := 77/200000;
k6 := 991/1000;
k7 := 3/8;
k8 := 104/5;
k9 := 145;
k10 := 357/100;
k11 := 53/50;
k12 := 409/100;
k13 := 35/2;
k14 := 14;
k15 := 168;
k16 := 5;
k17 := 123/125;
k18 := 43/500;
k19 := 611/1000;
k20 := 1;
k21 := 1;
k22 := 0;
k23 := 0;
k24 := 5490000;
k25 := 77/200000;
k26 := 209/10000;
k27 := 1;
operator diff$
odes := { diff(x7, t) = (1*k27*k18*x5 + (-1)*k27*k18*x7)/k27,
  diff(x8, t) = (1*k27*k18*x7 + (-1)*k27*k18*x8)/k27,
  diff(x9, t) = (1*k27*k18*x8 + (-1)*k27*k18*x9)/k27,
  diff(x10, t) = (1*k27*k18*x9 + (-1)*k27*k18*x10)/k27,
  diff(x11, t) = (1*k27*k19*x6 + (-1)*k27*k19*x11)/k27,
  diff(x12, t) = (1*k27*k19*x11 + (-1)*k27*k19*x12)/k27,
  diff(x13, t) = (1*k27*k19*x12 + (-1)*k27*k19*x13)/k27,
  diff(x14, t) = (1*k27*k19*x13 + (-1)*k27*k19*x14)/k27,
  diff(x1, t) = (1 - x3)*(1 - x4)*k26*x1*(1 - x1/k24) - (1 + x10)*(1 + x14)*k25*x1,
  diff(x2, t) = (1 + x10)*(1 + x14)*k25*x1 - k25*x2 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file