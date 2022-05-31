load_package f5$

                   load_package groebner$ torder({}, revgradlex)$

                   parameters := k1 = 1/200;
k2 = 7/10;
k3 = 7/10;
k4 = 21/10000;
k5 = 7/10000;
k6 = 1/20;
k7 = 0;
k8 = 0;
k9 = 0;
k10 = 10;
k11 = 1/100;
k12 = 10;
k13 = 29/50000;
k14 = 17/1000;
k15 = 1/50;
k16 = 3;
k17 = 189/1000;
k18 = 3000000;
k19 = 7/20;
k20 = 200000;
k21 = 86;
k22 = 1000;
k23 = 250;
k24 = 11627/59127;
k25 = 7/200;
k26 = 92390375/1544401932;
k27 = 5/258;
k28 = 0;
k29 = 125/43;
k30 = 150;
k31 = 1;$

                   operator diff$

                   odes := { diff(x1, t) = k5*(x3 + k6*k1)/(x3 + k1) - k6*k3*x1/((x3 + k6*k1)/(x3 + k1)),
  diff(x2, t) = k6*k3*x1/((x3 + k6*k1)/(x3 + k1)) - k17*x2,
  diff(x3, t) = k4*k13/k14*k18*(k9/k21 + k23/k21)/(k9/k21 + k16/k15)*x2/(1 + k13*k10/k14 + k11/(k12*k19)*(k8 + k20*x1/((k9/k21 + k23/k21)/(k9/k21 + k16/k15))))*(1 + k7/k22) - k2*(x3 + k6*k1)/(x3 + k1)*x3 }$
                   odes := for each o in odes collect part(o, 2)$ 

                   f5(odes);

                   end;