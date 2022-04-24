
module f5lp;

% The lp module provides the Labeled Polynomial interface --
% a special polynomial type to be used in the f5-style algorithms.
% The LabeledPolynomial object `p` is stored internally as a 4-item list:
%   {'lp, evaluation of `p`, signature index of `p`, signature monomial of `p`}
%
% The interface provides, in particular,
% functions `lp_eval` and `lp_sgn`, that return
% second and {third, fouth} items of the internal list respectively

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instantiates LabeledPolynomial from the given Polynomial `poly` and
% places garbage in the signature position
%
% Currently, used only in computing plain polynomial normal form
asserted inline procedure lp_LabeledPolynomial0(
                              poly: Polynomial): LabeledPolynomial;
  lp_LabeledPolynomial2(poly, {0, poly_zeroExp()});

% instantiates LabeledPolynomial from Polynomial and the leading index
asserted inline procedure lp_LabeledPolynomial1(
                              poly: Polynomial,
                              idx: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly, {idx, poly_zeroExp()});

% instantiates LabeledPolynomial from Polynomial and its signature
asserted inline procedure lp_LabeledPolynomial2(
                              poly: Polynomial,
                              sgn: List): LabeledPolynomial;
  'lbl . poly . sgn;

asserted inline procedure lp_getPoly(lp: LabeledPolynomial): Polynomial;
  cadr lp;

asserted inline procedure lp_getLabel(lp: LabeledPolynomial): List;
  cddr lp;

% Checks if LP is zero
%
% Zero LabeledPolynomial is represented as
%   ('lp, zero Polynomial, smth, smth)
%
% Same as lp_isSyzygy!?, in fact
asserted inline procedure lp_iszero!?(lp);
  poly_iszero!?(lp_eval(lp));

% return the evaluation of module element f
asserted inline procedure lp_setEval(f, ev): Polynomial;
  cadr f := ev;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% SIGNATURE MANIPULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return the signature {i, m}
asserted inline procedure lp_sgnInit(i: Integer, m: List): List;
  {i, m};

% return the signature of module element f
asserted inline procedure lp_sgn(f: LabeledPolynomial): List;
  lp_getLabel(f);

% sgn1 == sgn2 ?
asserted procedure lp_eqSignature(sgn1: List, sgn2: List);
  car sgn1 #= car sgn2 and lp_eqList(cdr sgn1, cdr sgn2);

% return the index of the signature of module element f
asserted inline procedure lp_index(f: LabeledPolynomial): Integer;
  car lp_getLabel(f);

% return the monom of the signature of module element f
asserted inline procedure lp_monom(f: LabeledPolynomial): Integer;
  cadr lp_getLabel(f);

% return the evaluation of module element f
asserted inline procedure lp_eval(f: LabeledPolynomial): Polynomial;
   lp_getPoly(f);

% multiply signature sgn by monomial ev
asserted procedure lp_multSignature(sgn: List, ev: List);
  begin scalar  term;
        integer idx;
    idx  := car sgn;
    term := cadr sgn;
    return lp_sgnInit(idx, poly_sumExp(term, ev))
  end;

% multiply signature sgn by monomial ev
asserted procedure lp_sgnMonom(sgn);
  cadr sgn;

% multiply signature sgn by monomial ev
asserted procedure lp_sgnIndex(sgn);
  car sgn;

% multiply signature sgn by monomial ev
asserted procedure lp_multSgn(sgn: List, ev: List);
  begin scalar  term;
        integer idx;
    idx  := car sgn;
    term := cadr sgn;
    return {idx, poly_sumExp(term, ev)}
  end;

% compare signatures sgn1 vs. sgn2 with the
% position over term extension
asserted procedure lp_sgnCmp(sgn1, sgn2);
  begin integer idx1, idx2;
        scalar  ev1, ev2;
    idx1 := car sgn1;
    ev1  := cadr sgn1;
    idx2 := car sgn2;
    ev2  := cadr sgn2;
    return if idx1 #= idx2 then poly_cmpExp(ev1, ev2)
      else (idx1 #< idx2)
  end;

% checks if lp represents a syzygy
asserted inline procedure lp_isSyzygy(lp: LabeledPolynomial): Boolean;
  poly_iszero!?(lp_eval(lp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% COEFFICIENT MANIPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mainly, these fall back to simple polynomial coefficient operations

% asserted procedure lp_normalizeInplace(f: LabeledPolynomial);
%  poly_normalizeInplace(lp_eval(f));

asserted procedure lp_normalize(f: LabeledPolynomial);
  lp_LabeledPolynomial2(poly_normalize(lp_eval(f)), lp_sgn(f));

asserted procedure lp_scaleDenominators(f: LabeledPolynomial);
  lp_LabeledPolynomial2(poly_scaleDenominators(lp_eval(f)), lp_sgn(f));

asserted procedure lp_scaleDenominatorsInplace(f: LabeledPolynomial);
  lp_LabeledPolynomial2(poly_scaleDenominatorsInplace(lp_eval(f)), lp_sgn(f));

asserted procedure lp_reduceCoeffs(f: LabeledPolynomial, prime);
  lp_LabeledPolynomial2(poly_reduceCoeffs(lp_eval(f), prime), lp_sgn(f));

asserted procedure lp_reconstructCoeffs(f: LabeledPolynomial, prime);
  lp_LabeledPolynomial2(poly_reconstructCoeffs(lp_eval(f), prime), lp_sgn(f));

asserted procedure lp_crtCoeffs(polyaccum, modulo, polycomp, prime);
  begin scalar underlyingPoly;
    underlyingPoly := poly_crtCoeffs(lp_eval(polyaccum), modulo, lp_eval(polycomp), prime);
    return lp_LabeledPolynomial2(underlyingPoly, lp_sgn(polyaccum))
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% l1 = l2 elementwise ?
asserted procedure lp_eqList(l1, l2);
  if null l1 then
    t
  else
    car l1 #= car l2 and lp_eqList(cdr l1, cdr l2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LabeledPolynomial sorting ad-hoc

% return true if lead(poly1) < lead(poly2)
asserted procedure lp_cmpLPLead(lp1, lp2);
  poly_cmpPolyLead(lp_eval(lp1), lp_eval(lp2));

% return true if lead(poly1) > lead(poly2)
asserted procedure lp_cmpLPLeadReverse(lp1, lp2);
  poly_cmpPolyLead(lp_eval(lp2), lp_eval(lp1));

% return true if lead(poly2) < lead(poly1)
% asserted procedure lp_leadTotalDegreeCmp(lp1, lp2);
%  poly_cmpPolyLead(lp_eval(lp2), lp_eval(lp1));

% return true if total_degree(lead(lp1)) < total_degree(lead(lp2))
asserted procedure lp_leadTotalDegreeCmp(lp1, lp2);
  poly_tdegCmp(poly_leadExp(lp_eval(lp1)), poly_leadExp(lp_eval(lp2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trst lp_LabeledPolynomial;
% trst lp_spoly;
% trst lp_spolyMultSignatures;
% trst lp_spolyCofactors;
% trst lp_multSignature;

% trst lp_principalSyzygy;
% trst lp_tryReduce1;
% trst lp_tryReduce;
% trst lp_normalForm;
% trst lp_isSingularlyTopReducible;

% trst lp_principalSyzygy;

endmodule;


end;  % of file
