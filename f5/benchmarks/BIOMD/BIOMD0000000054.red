load_package f5$

                   load_package groebner$ torder({}, revgradlex)$

                   parameters := k1 = 1;
k2 = 1/100;
k3 = 121/1000;
k4 = 100;
k5 = 1/5;
k6 = 337/25;
k7 = 1/50;
k8 = 1/100;
k9 = 6/5;
k10 = -1;
k11 = 1;$

                   operator diff$

                   odes := { diff(x1, t) = (1*k11*k3*k4 + (-3)*k11*k5*x1*(x3 + 3*x2 - (6*x3*x2 - 3*x2^2 + x3^2)^(1/2))/6)/k11,
  diff(x2, t) = ((-1)*k11*k5*x1*(x3 + 3*x2 - (6*x3*x2 - 3*x2^2 + x3^2)^(1/2))/6 + 1*k11*k6*((x3 + 3*x2 - (6*x3*x2 - 3*x2^2 + x3^2)^(1/2))/6)^(13/25)*((7*x3 - 3*x2 - (6*x3*x2 - 3*x2^2 + x3^2)^(1/2))/6)^(41/100) + (-1)*k11*2*k7)/k11,
  diff(x3, t) = 1*k11*k7*(1 - k8*((x3 + 3*x2 - (6*x3*x2 - 3*x2^2 + x3^2)^(1/2))/6)^k9*((7*x3 - 3*x2 - (6*x3*x2 - 3*x2^2 + x3^2)^(1/2))/6)^k10)/k11 }$
                   odes := for each o in odes collect part(o, 2)$ 

                   f5(odes);

                   end;