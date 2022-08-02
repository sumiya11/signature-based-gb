load_package f5$
1;
torder({}, revgradlex)$

k1 := 281441;
k2 := 171977;
k3 := 412/73;
k4 := 40758549/3650000;
k5 := 1/730;
k6 := 1/25258;
k7 := 1/6159375000;
k8 := 9137/2635182;
k9 := 9137/5270364;
k10 := 180000;
k11 := 120000;
k12 := 850;
k13 := 1/16425000000;
k14 := 22841/6587955;
k15 := 0;
k16 := 0;
k17 := 3/1000;
k18 := 0;
k19 := 0;
k20 := 171977;
k21 := 1;
k22 := 1;
operator diff$
odes := { diff(x1, t) = (1*k21*2060/365 + (-1)*k21*1/(2*365)*x1 + (-1)*k21*45685/(439197*30)*x1*x2)/k21,
  diff(x2, t) = (1*k21*45685/(439197*30)*x1*x2 + (-1)*k21*45682/(439197*30)*x2 + (-1)*k21*1/(2*365)*x2)/k21,
  diff(x3, t) = (1*k22*237/10*k20/365000 + (-1)*k22*1/(346/5*365)*x3 + (-1)*16/(108*1000000*1825/2)*x3*x2)/k22,
  diff(x4, t) = (1*16/(108*1000000*1825/2)*x3*x2 + (-1)*k22*6/(108*1000000*1825/2)*x4 + (-1)*k22*1/(346/5*365)*x4)/k22 }$
odes := for each o in odes collect part(o, 2)$

gb := f5(odes)$

end; % of file