load_package f5$
on f5sugar;
torder({}, revgradlex)$

k1 := 6987/5000000000;
k2 := 85551/1000000;
k3 := 2211/100000;
k4 := 15783/2000000;
k5 := 7617/20000000;
k6 := 27761/500000;
k7 := 15639/100000000000;
k8 := 22493/250;
k9 := 829;
k10 := 6018000000000000/200003009;
k11 := 400000000000000000000/200003009;
k12 := 326;
k13 := 3009/200000000;
k14 := 72592125/200003009;
k15 := 326/829;
k16 := 10;
k17 := 15599/100000000;
k18 := 29553/100000000;
k19 := 24125;
k20 := 12913/100000000;
k21 := 15789/50000000;
k22 := 1;
operator diff$
odes := { diff(x1, t) = (-(k17 + k1))*x1 + k18*x2 + k2*k13*x3,
  diff(x2, t) = k17*x1 - (k18 + k1)*x2 + k2*k13*x4,
  diff(x3, t) = k1/k13*x1 - (k2 + k17 + k3)*x3 + k18*x4 + k4/(k19*k13/(1 + k13))*x5,
  diff(x4, t) = k1/k13*x2 + k17*x3 - (k18 + k2)*x4,
  diff(x5, t) = k3*k19*k13/(1 + k13)*x3 - (k4 + k20)*x5 + k21*x6 + k6*k12/k9*x7 - k5*(k8 - x7)*x5,
  diff(x6, t) = k20*x5 - k21*x6 + k7*k12/k9*x7,
  diff(x7, t) = k5/(k12/k9)*(k8 - x7)*x5 - (k6 + k7)*x7 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file