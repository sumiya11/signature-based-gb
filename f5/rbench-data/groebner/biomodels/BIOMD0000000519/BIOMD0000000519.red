load_package groebner$
torder({}, revgradlex)$

parameters := k1 := 1/10;
k2 := 420467092599869/1000000000000000;
k3 := 55069267386123/50000000000000;
k4 := 1495588966300553/20000000000000;
k5 := 120;
k6 := 36779737044933/20000000000000;
k7 := 389066458140967532339483328595180417255401/320000000000000000000000000000000000000000;
k8 := 9018759018759/62500000000000;
k9 := 53481240981241/62500000000000;
k10 := 147002558002553/25000000000000;
k11 := 121454266376232554106270728286370779423/2000000000000000000000000000000000000;
k12 := 172622515189057/1000000000000000;
k13 := 827377484810943/1000000000000000;
k14 := 1;
operator diff$
odes := { diff(x1, t) = ((-1)*k1*x1 + 1*(1 - k9)*k6*x1*(1 - (x1 + x2 + x3)/k5))/k14,
  diff(x2, t) = (1*k9*k6*x1*(1 - (x1 + x2 + x3)/k5) + (-1)*k2*x2 + 1*(1 - k13)*k10*x2*(1 - (x1 + x2 + x3)/k5))/k14,
  diff(x3, t) = (1*k13*k10*x2*(1 - (x1 + x2 + x3)/k5) + (-1)*k3*x3)/k14 }$
odes := for each o in odes collect part(o, 2)$

gb := groebner(odes)$

end; % of file