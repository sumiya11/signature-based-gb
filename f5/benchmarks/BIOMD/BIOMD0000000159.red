load_package f5$

                   load_package groebner$ torder({}, revgradlex)$

                   parameters := k1 = 3/10;
k2 = 1;
k3 = 0;
k4 = 16/5;
k5 = 2/5;
k6 = 1/10;
k7 = 1/10;
k8 = 1;$

                   operator diff$

                   odes := { diff(x1, t) = (1*k8*k1*k2 + (-1)*k8*k3*x1 + (-1)*k8*k4*x2*x1)/k8,
  diff(x2, t) = (1*k8*k7*x3 + (-1)*k8*k6*x2)/k8,
  diff(x3, t) = (1*k8*k5*x1*k2 + (-1)*k8*k7*x3)/k8 }$
                   odes := for each o in odes collect part(o, 2)$ 

                   f5(odes);

                   end;