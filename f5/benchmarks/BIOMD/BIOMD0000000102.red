load_package f5$

                   load_package groebner$ torder({}, revgradlex)$

                   parameters := k1 = 1/500;
k2 = 1/10;
k3 = 1/200000;
k4 = 7/20000;
k5 = 1/5000;
k6 = 1/5000;
k7 = 1/20000;
k8 = 7/2000;
k9 = 1/500;
k10 = 1/10;
k11 = 1/1000;
k12 = 1/1000;
k13 = 1/1000;
k14 = 1/1000;
k15 = 1/1000;
k16 = 1/1000;
k17 = 1/1000;
k18 = 1/1000;
k19 = 1/500;
k20 = 1/10;
k21 = 1/500;
k22 = 1/10;
k23 = 3/1000;
k24 = 1/1000;
k25 = 1/1000;
k26 = 1/50;
k27 = 1/1000;
k28 = 1/50;
k29 = 1/1000;
k30 = 1/25;
k31 = 1/1000;
k32 = 1/1000;
k33 = 1/1000;
k34 = 1/1000;
k35 = 1/5;
k36 = 1/1000;
k37 = 1/1000;
k38 = 1/1000;
k39 = 1/1000;
k40 = 1/1000;
k41 = 1/1000;
k42 = 1;$

                   operator diff$

                   odes := { diff(x1, t) = ((-1)*k42*(k1*x1*x2 - k2*x6) + (-1)*k42*(k19*x3*x1 - k20*x5) + (-1)*k42*(k9*x11*x1 - k10*x12) + (-1)*k42*(k21*x10*x1 - k22*x13) + 1*k42*(k26 - k25*x1))/k42,
  diff(x2, t) = ((-1)*k42*(k1*x1*x2 - k2*x6) + (-1)*k42*(k11*x2*x4 - k12*x3) + (-1)*k42*k5*x2*x8 + 1*k42*(k28 - k27*x2))/k42,
  diff(x3, t) = (1*k42*(k11*x2*x4 - k12*x3) + (-1)*k42*(k19*x3*x1 - k20*x5) + (-1)*k42*k31*x3)/k42,
  diff(x4, t) = ((-1)*k42*(k11*x2*x4 - k12*x3) + (-1)*k42*(k13*x6*x4 - k14*x5) + (-1)*k42*(k23*x8*x4 - k24*x9) + (-1)*k42*(k15*x11*x4 - k16*x10) + (-1)*k42*(k17*x12*x4 - k18*x13) + 1*k42*(k30 - k29*x4))/k42,
  diff(x5, t) = (1*k42*(k13*x6*x4 - k14*x5) + 1*k42*(k19*x3*x1 - k20*x5) + (-1)*k42*k32*x5)/k42,
  diff(x6, t) = (1*k42*(k1*x1*x2 - k2*x6) + (-1)*k42*(k13*x6*x4 - k14*x5) + (-1)*k42*k6*x6*x8 + (-1)*k42*k33*x6)/k42,
  diff(x7, t) = ((-1)*k42*k3*x7*x2 + (-1)*k42*k4*x7*x6 + (-1)*k42*k7*x7*x11 + (-1)*k42*k8*x7*x12 + 1*k42*(k35 - k34*x7))/k42,
  diff(x8, t) = (1*k42*k3*x7*x2 + 1*k42*k4*x7*x6 + (-1)*k42*(k23*x8*x4 - k24*x9) + 1*k42*k7*x7*x11 + 1*k42*k8*x7*x12 + (-1)*k42*k36*x8)/k42,
  diff(x9, t) = (1*k42*(k23*x8*x4 - k24*x9) + (-1)*k42*k37*x9)/k42,
  diff(x10, t) = (1*k42*(k15*x11*x4 - k16*x10) + (-1)*k42*(k21*x10*x1 - k22*x13) + (-1)*k42*k38*x10)/k42,
  diff(x11, t) = (1*k42*k5*x2*x8 + (-1)*k42*(k9*x11*x1 - k10*x12) + (-1)*k42*(k15*x11*x4 - k16*x10) + (-1)*k42*k39*x11)/k42,
  diff(x12, t) = (1*k42*k6*x6*x8 + 1*k42*(k9*x11*x1 - k10*x12) + (-1)*k42*(k17*x12*x4 - k18*x13) + (-1)*k42*k40*x12)/k42,
  diff(x13, t) = (1*k42*(k17*x12*x4 - k18*x13) + 1*k42*(k21*x10*x1 - k22*x13) + (-1)*k42*k41*x13)/k42 }$
                   odes := for each o in odes collect part(o, 2)$ 

                   f5(odes);

                   end;