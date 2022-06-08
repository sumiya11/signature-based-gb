% symmetric nbody-4 system in revgradlex
% characteristic 0
% 0 dim
%
% PoSSo test suite
% https://www-sop.inria.fr/saga/POL/BASE/2.multipol/centralpos.html

load_package f5;

system := {
2*Theta-ss*s-2*ff+2,
8*ff*p-4*p*ss-2*ff*s^2+ss*s^2+4*Theta-2*ss*s, -2*s-4*p+s^2+f+1,
s*Theta^2-p*s*pp-p*ss*Theta-2, p^3*pp^2-1, ff^2*f^3-1
}$

vars := {s,p,ss,pp,TH,ff,f}$
torder(vars, revgradlex)$

gb := f5(system)$

end;
