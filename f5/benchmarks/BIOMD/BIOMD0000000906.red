load_package f5$

                   load_package groebner$ torder({}, revgradlex)$

                   parameters := k1 = 9/10;
k2 = 1/20;
k3 = 1/5;
k4 = 1/25;
k5 = 4/5;
k6 = 3/10;
k7 = 1/10;
k8 = 12/5;
k9 = 1/10;
k10 = 1;$

                   operator diff$

                   odes := { diff(x1, t) = (1*k10*k1*x1 + (-1)*k10*(k2*x2*x1 + k3*x1*x1))/k10,
  diff(x2, t) = (1*k10*(k4 + k6*x2*x1) + (-1)*k10*(k5*x2 + k7*k2*x2*x1))/k10,
  diff(x3, t) = (1*k10*k8*x1 + (-1)*k10*k9*x3)/k10 }$
                   odes := for each o in odes collect part(o, 2)$ 

                   f5(odes);

                   end;