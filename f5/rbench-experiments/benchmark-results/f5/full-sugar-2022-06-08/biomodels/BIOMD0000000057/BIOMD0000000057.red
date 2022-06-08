load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 0;
k2 := 2221/7265;
k3 := 16/25;
k4 := 3/25;
k5 := 17/10;
k6 := 4/5;
k7 := 10;
k8 := 1/40;
k9 := 10761/7265;
k10 := 187/5;
k11 := 17/10;
k12 := 72204/3235;
k13 := 1/25;
k14 := 7/5;
k15 := 149/5;
k16 := 5/2;
k17 := 547/10;
k18 := 6017/64700;
k19 := 11/100;
k20 := 492580/647;
k21 := 4;
k22 := 4707;
k23 := 1791/12650;
k24 := 27/50;
k25 := 57/5;
k26 := 2221/1265;
k27 := 10;
k28 := 1;
k29 := 10;
k30 := 1/25;
k31 := 4/5;
k32 := 149/5;
k33 := 1/25;
k34 := 4/5;
k35 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k28*((k10*k8 + k11*k7)/(k8 + k7*(1 + k8/k4))*k29*x1 - (k14 + k16*k7)/(1 + k7/k17)*x2) + (-1)*k28*((k3*k4 + k5)*k7/(k4 + k7*(1 + k4/k8))*x1 - (k30 + k31)*x3))/k28,
  diff(x2, t) = (1*k28*((k10*k8 + k11*k7)/(k8 + k7*(1 + k8/k4))*k29*x1 - (k14 + k16*k7)/(1 + k7/k17)*x2) + (-1)*k28*(k19*k17/(k17 + k7)*x2 - k32*x4) + (-1)*k28*((k21*k17 + k22)*k7/(k17 + k7)*x2 - k4*(k24 + k25)/(k4 + k7)*x5))/k28,
  diff(x3, t) = 1*k28*((k3*k4 + k5)*k7/(k4 + k7*(1 + k4/k8))*x1 - (k30 + k31)*x3)/k28,
  diff(x4, t) = 1*k28*(k19*k17/(k17 + k7)*x2 - k32*x4)/k28,
  diff(x5, t) = (1*k28*((k21*k17 + k22)*k7/(k17 + k7)*x2 - k4*(k24 + k25)/(k4 + k7)*x5) + (-1)*k28*((k3*k4 + k5)*k7/(k4 + k7)*x5 - (k33 + k34)*x6))/k28,
  diff(x6, t) = 1*k28*((k3*k4 + k5)*k7/(k4 + k7)*x5 - (k33 + k34)*x6)/k28 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file