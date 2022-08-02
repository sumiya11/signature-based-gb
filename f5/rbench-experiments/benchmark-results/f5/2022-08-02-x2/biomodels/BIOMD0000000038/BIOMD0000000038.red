load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 1960;
k3 := 480000;
k4 := 108000;
k5 := 294;
k6 := 14000;
k7 := 14000;
k8 := 84000;
k9 := 3360;
k10 := 21960;
k11 := 21960;
k12 := 4392;
k13 := 3384;
k14 := 880;
k15 := 880;
k16 := 2640;
k17 := 960;
k18 := 260;
k19 := 389;
k20 := 4800;
k21 := 27/5000;
k22 := 2800;
k23 := 900;
k24 := 50;
k25 := 500;
k26 := 5;
k27 := 50;
k28 := 40;
k29 := 10;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*(k2*k22*x1 - k3*x2) + 1*k1*(k8*x5 - k9*x1*x6))/k1,
  diff(x2, t) = (1*k1*(k2*k22*x1 - k3*x2) + (-1)*k1*(k4*x2 - k5*k23*x3))/k1,
  diff(x3, t) = (1*k1*(k4*x2 - k5*k23*x3) + (-1)*k1*(k6*x3*x4 - k7*x5))/k1,
  diff(x4, t) = ((-1)*k1*(k6*x3*x4 - k7*x5) + 1*k1*(k12*x8 - k13*x4*x9))/k1,
  diff(x5, t) = (1*k1*(k6*x3*x4 - k7*x5) + (-1)*k1*(k8*x5 - k9*x1*x6))/k1,
  diff(x6, t) = (1*k1*(k8*x5 - k9*x1*x6) + (-1)*k1*(k10*x6*x7 - k11*x8))/k1,
  diff(x7, t) = ((-1)*k1*(k10*x6*x7 - k11*x8) + 1*k1*(k16*x11 - k17*x7*x12))/k1,
  diff(x8, t) = (1*k1*(k10*x6*x7 - k11*x8) + (-1)*k1*(k12*x8 - k13*x4*x9))/k1,
  diff(x9, t) = (1*k1*(k12*x8 - k13*x4*x9) + (-1)*k1*(k14*x9*x10 - k15*x11))/k1,
  diff(x10, t) = ((-1)*k1*(k14*x9*x10 - k15*x11) + 1*k1*(k20*x13 - k21*x10*k24))/k1,
  diff(x11, t) = (1*k1*(k14*x9*x10 - k15*x11) + (-1)*k1*(k16*x11 - k17*x7*x12))/k1,
  diff(x12, t) = (1*k1*(k16*x11 - k17*x7*x12) + (-1)*k1*(k18*x12*k25 - k19*x13))/k1,
  diff(x13, t) = (1*k1*(k18*x12*k25 - k19*x13) + (-1)*k1*(k20*x13 - k21*x10*k24))/k1,
  diff(x14, t) = 0,
  diff(x15, t) = 0,
  diff(x16, t) = 0,
  diff(x17, t) = 0 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file