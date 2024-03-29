load_package f5$
torder({}, revgradlex)$

k1 := 1;
k2 := 1;
k3 := 1267031539/100000000;
k4 := 213430191/200000000;
k5 := 12222573/200000000;
k6 := 727706671/100000000;
k7 := 20105521/200000000;
k8 := 522645841/200000000;
k9 := 137763703/1000000000;
k10 := 746666581/500000000;
k11 := 19305821/125000000;
k12 := 10725099/250000000;
k13 := 179809059/200000000;
k14 := 572065273/500000000;
k15 := 355490081/1000000000;
k16 := 21879693/500000000;
k17 := 134371419/1000000000;
k18 := 15336713/200000000;
k19 := 304695409/1000000000;
k20 := 192119917/1000000000;
k21 := 445547231/1000000000;
k22 := 40272103/200000000;
k23 := 54570911/1000000000;
k24 := 23306949/250000000;
k25 := 121370929/1000000000;
k26 := 11186909/250000000;
k27 := 125873837/1000000000;
k28 := 1269699/40000000;
k29 := 30471301/500000000;
k30 := 21/50;
k31 := 37/100;
k33 := 0;
k32 := 100;
operator diff$
odes := { diff(x1, t) = ((-1)*x1*k3 + 1*x4*k6 + (-1)*x1*k8 + 1*x5*k9 + (-1)*x1*k10 + 1*x6*k11 + (-1)*x1*k12 + (-1)*x1*k13 + (-1)*x1*k14 + (-1)*x1*k17 + 1*x11*k18 + (-1)*x1*k19 + 1*x12*k20 + (-1)*x1*k21 + 1*x13*k22 + (-1)*x1*k23 + 1*x14*k24 + (-1)*x1*k25 + (-1)*x1*k26 + 1*x16*k27 + (-1)*x1*k28 + 1*x17*k29 + 1*k30*x7)/k2,
  diff(x2, t) = (1*x1*k3 + (-1)*x2*k4 + (-1)*x2*k7)/k2,
  diff(x3, t) = (1*x2*k4 + (-1)*x3*k5)/k2,
  diff(x4, t) = (1*x3*k5 + (-1)*x4*k6 + 1*x2*k7)/k2,
  diff(x5, t) = (1*x1*k8 + (-1)*x5*k9)/k2,
  diff(x6, t) = (1*x1*k10 + (-1)*x6*k11)/k2,
  diff(x7, t) = (1*x1*k12 + (-1)*k30*x7)/k2,
  diff(x8, t) = (1*x1*k14 + (-1)*x8*k16)/k2,
  diff(x9, t) = (1*x1*k13 + (-1)*x9*k15)/k2,
  diff(x10, t) = (1*x9*k15 + 1*x8*k16 + 1*k31*x15)/k1,
  diff(x11, t) = (1*x1*k17 + (-1)*x11*k18)/k2,
  diff(x12, t) = (1*x1*k19 + (-1)*x12*k20)/k2,
  diff(x13, t) = (1*x1*k21 + (-1)*x13*k22)/k2,
  diff(x14, t) = (1*x1*k23 + (-1)*x14*k24)/k2,
  diff(x15, t) = (1*x1*k25 + (-1)*k31*x15)/k2,
  diff(x16, t) = (1*x1*k26 + (-1)*x16*k27)/k2,
  diff(x17, t) = (1*x1*k28 + (-1)*x17*k29)/k2 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$
in "~/signature-based-gb/f5/sed.red";

end; % of file