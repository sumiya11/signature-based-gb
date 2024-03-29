% nbody-4 system in revgradlex
% characteristic 0
% 0 dim
%
% PoSSo test suite
% https://www-sop.inria.fr/saga/POL/BASE/2.multipol/centralpos.html

load_package f5;
on f5sugar;

system := {
	(b-d)*(bb-dd)-2*ff+2,
	(b-d)*(bb+dd-2*ff)+2*(bb-dd),
	(b-d)^2-2*(b+d)+f+1,
	bb^2*b^3-1,
	dd^2*d^3-1,
	ff^2*f^3-1
}$

vars := {bb,dd,ff,b,d,f}$
torder(vars, revgradlex)$

gb := f5(system)$

end;
