load_package f5$
1;
torder({}, revgradlex)$

k1 := 1;
k2 := 0;
k3 := 0;
k4 := 1;
k5 := 0;
k6 := 0;
k7 := 11773/2500;
k8 := 63/1000;
k9 := 1/20;
k10 := 0;
k11 := 741699/2500000;
k12 := 1/20;
k13 := 9010000;
k14 := 9010000;
k15 := 48000000;
k16 := 1335000000;
k17 := 11773/2500;
k18 := 63/1000;
k19 := 1/20;
k20 := 303/400;
k21 := 3917/10000;
k22 := 643/10000;
k23 := 4797/10000;
k24 := 1237/2000;
k25 := 161/500;
k26 := 2967/500;
k27 := 1/20;
k28 := 1/20;
k29 := 6079/10000;
k30 := 61/125;
k31 := 957/5000;
k32 := 15283/10000;
k33 := 1941/10000;
k34 := 1/20;
k35 := 5753/10000;
k36 := 5157/10000;
k37 := 2189/10000;
k38 := 11773/2500;
k39 := 4797/10000;
k40 := 5753/10000;
k41 := 6079/10000;
k42 := 303/400;
k43 := 15283/10000;
k44 := 2967/500;
k45 := 11773/2500;
k46 := 1335000000;
k47 := 48000000;
k48 := 9010000;
k49 := 0;
k50 := 0;
k51 := 0;
k52 := 0;
k53 := 1;
k54 := 1;
k55 := 63/1000;
k56 := 1237/2000;
k57 := 5157/10000;
k58 := 61/125;
k59 := 3917/10000;
k60 := 1941/10000;
k61 := 1/20;
k62 := 63/1000;
k63 := 1/20;
k64 := 161/500;
k65 := 2189/10000;
k66 := 957/5000;
k67 := 643/10000;
k68 := 1/20;
k69 := 1/20;
k70 := 1/20;
k71 := 0;
k72 := 1;
k73 := 9010258;
operator diff$
odes := { diff(x1, t) = (-1)*k72*k38*k55*x1*x2/(k54*k48 + k50*k47 + k49*k46)/k72,
  diff(x2, t) = (1*k72*k38*k55*x1*x2/(k54*k48 + k50*k47 + k49*k46) + (-1)*k72*(k54*(k53*k62 + k52*k59 + k51*k56) + k50*(k53*k61 + k52*k58) + k49*(k53*k60 + k52*k57))*x2)/k72,
  diff(x3, t) = (1*k72*(k54*(k53*k62 + k52*k59 + k51*k56) + k50*(k53*k61 + k52*k58) + k49*(k53*k60 + k52*k57))*x2 + (-1)*k72*(k63 + (1 - k63)*k71)*x3)/k72,
  diff(x4, t) = 1*k72*(k63 + (1 - k63)*k71)*x3/k72 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file