module f5lp;

% The lp module provides the Labeled Polynomial interface --
% a special polynomial type to be used in the f5-style algorithms.
% The LabeledPolynomial object `p` is stored internally as a 4-item list:
%   {'lp, evaluation of `p`, signature index of `p`, signature monomial of `p`}
%
% The interface provides, in particular,
% functions `lp_eval` and `lp_sgn`, that return
% second and {third, fouth} items of the internal list respectively

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
asserted inline procedure lp_sgnIndex(f: LabeledPolynomial): Integer;
  car lp_getLabel(f);

% return the monom of the signature of module element f
asserted inline procedure lp_sgnMonom(f: LabeledPolynomial): Integer;
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

% compare signatures sgn1 vs. sgn2 with the
% position over term extension
asserted procedure lp_potCmpSignature(sgn1: List, sgn2: List);
  begin integer idx1, idx2;
        scalar  ev1, ev2;
    idx1 := car sgn1;
    ev1  := cadr sgn1;
    idx2 := car sgn2;
    ev2  := cadr sgn2;
    return if idx1 equal idx2 then poly_cmpExp(ev1, ev2)
      else (idx1 #> idx2)
  end;

% minimal signature based on the pot comparison strategy
asserted inline procedure lp_potMinSignature(sgn1: List, sgn2: List): List;
  if lp_potCmpSignature(sgn1, sgn2) then sgn1 else sgn2;

% maximal signature based on pot comparison strategy
asserted inline procedure lp_potMaxSignature(sgn1: List, sgn2: List): List;
  if lp_potCmpSignature(sgn1, sgn2) then sgn2 else sgn1;

% compare LabeledPolynomials as module elems lp1 vs. lp2 with
% the position over term extension
asserted procedure lp_potCmpLP(lp1: LabeledPolynomial, lp2: LabeledPolynomial): Boolean;
  begin scalar sgn1, sgn2;
    sgn1 := lp_sgn(lp1);
    sgn2 := lp_sgn(lp2);
    return lp_potCmpSignature(sgn1, sgn2)
  end;

% checks if lp represents a syzygy
asserted inline procedure lp_isSyzygy(lp: LabeledPolynomial): Boolean;
  poly_iszero!?(lp_eval(lp));

% checks if signature sgn1 divides signature sgn2
asserted procedure lp_sgnDivides(sgn1, sgn2): Boolean;
  (car sgn1 equal car sgn2) and poly_divExp!?(cadr sgn1, cadr sgn2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constructs the principal syzygy of LabeledPolynomials (fi, fj)
% as a LabeledPolynomial itself.
% Assumes i != j.
% Principal syzygy - a syz of form ej*fi - ei*fj,
% with the evaluation equal to fj*fi - fi*fj == 0,
% and the signature = { min(j, i), lt(fi) if j < i else lt(fj) };
asserted procedure lp_principalSyzygy(fi: LabeledPolynomial,
                                      fj: LabeledPolynomial): LabeledPolynomial;
  begin scalar i, j, ei, ej, sgn;
      i := lp_sgnIndex(fi);
      j := lp_sgnIndex(fj);
      ei := poly_leadExp(lp_eval(fi));
      ej := poly_leadExp(lp_eval(fj));
      sgn := if i < j then {i, ej} else {j, ei};

      % TODO: omit poly_zero ?
      return lp_LabeledPolynomial2(poly_zero(), sgn)
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% COEFFICIENT MANIPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mainly, these fall back to simple polynomial coefficient operations

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

% try to reduce module elem f with module elem g.
% Returns success flag with the result of reduction
asserted procedure lp_tryReduce1(f: LabeledPolynomial,
                                 g: LabeledPolynomial,
                                 isregular);
  begin scalar fexps, fcoeffs, fe, fc, zeroexp, gexplead,
                gcoefflead, flag, fevl, gevl, fsgn, gsgn,
                ans, zeroexp, gmult;

    % note that isregular=nil breaks the signature invariant
    % (storage of signature may become garbage as well)
    fevl := lp_eval(f);
    gevl := lp_eval(g);
    fsgn := lp_sgn(f);
    gsgn := lp_sgn(g);

    zeroexp    := poly_zeroExp();
    gexplead   := poly_leadExp(gevl);
    gcoefflead := poly_leadCoeff(gevl);

    fexps   := poly_getExps(fevl);
    fcoeffs := poly_getCoeffs(fevl);

    ans   := fevl;

    while fexps do <<
      fe . fexps   := fexps;
      fc . fcoeffs := fcoeffs;
      if not flag then <<
        if poly_divExp!?(gexplead, fe) then <<
          gmult := poly_subExp(fe, gexplead);
          if (not isregular) or not (fsgn equal lp_multSignature(gsgn, gmult)) then <<
            flag := t;
            ans := poly_paircomb(fevl, zeroexp, fc, gevl, gmult, gcoefflead);
          >>
        >>
      >>
    >>;

    return flag . lp_LabeledPolynomial2(ans, fsgn)
  end;

% try to reduce module elem f once with each element of module G,
% returns success flag with the result of reduction
asserted procedure lp_tryReduce(f: LabeledPolynomial, basis: List, isregular);
  begin scalar flag1, flag2, reducer;
    return if null basis then
      nil . f
    else <<
      reducer := car basis;
      flag1 . f := lp_tryReduce1(f, reducer, isregular);
      flag2 . f := lp_tryReduce(f, cdr basis, isregular);
      (flag1 or flag2) . f
    >>
  end;

% compute the normal form of f w.r.t. basis and return it
% If isregular is set, performs only regular reductions,
% otherwise reductions may corrupt the signature of f
%
% Contract: use isregular=nil only in the final interreduction
%           and in computing non-significant normal forms,
%           when corrupted signatures are not a problem
asserted procedure lp_normalForm(
                              f: LabeledPolynomial,
                              basis: List,
                              isregular): LabeledPolynomial;
  begin scalar state;
    state := t;
    if lp_iszero!?(f) then
      state := nil;
    while state do <<
      state . f := lp_tryReduce(f, basis, isregular);
      state := state and (not lp_iszero!?(f))
    >>;
    return f
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% l1 = l2 elementwise ?
asserted procedure lp_eqList(l1, l2);
  if null l1 then
    t
  else
    car l1 #= car l2 and lp_eqList(cdr l1, cdr l2);

% Checks if f is is singularly-top-reducible by some element of G
asserted procedure lp_isSingularlyTopReducible(
                                      f: LabeledPolynomial,
                                      basis: List);
  begin scalar g, fevl, gevl, ft, gt, flead, glead, m, ans;
    fevl  := lp_eval(f);
    ft    := lp_sgn(f);
    flead := poly_leadExp(fevl);
    for each g in basis do <<
      gevl  := lp_eval(g);
      gt    := lp_sgn(g);
      glead := poly_leadExp(gevl);
      if poly_divExp!?(glead, flead) then <<
        m := poly_subExp(flead, glead);
        if lp_eqSignature(ft, lp_multSignature(gt, m)) then
          ans := t
      >>
    >>;
    return ans
  end;

% compute cofactors a, b of Spoly of f and g: a*f - b*g
asserted procedure lp_spolyCofactors(
                          f: LabeledPolynomial,
                          g: LabeledPolynomial);
  begin scalar p1, p2, e1, e2, elcm, mult1, mult2;
    p1 := lp_eval(f);
    p2 := lp_eval(g);
    e1 := poly_leadExp(p1);
    e2 := poly_leadExp(p2);

    elcm := poly_lcmExp(e1, e2);

    mult1 := poly_subExp(elcm, e2);
    mult2 := poly_subExp(elcm, e1);

    return mult1 . mult2
  end;

% compute signatures of p1 and p2
% after multiplication by cofactors of Spolynomial of p1 and p2
asserted procedure lp_spolyMultSignatures(
                                p1: LabeledPolynomial,
                                p2: LabeledPolynomial): List;
  begin scalar mi, mj, si, sj;
    mi . mj := lp_spolyCofactors(p1, p2);
    si := lp_sgn(p1);
    sj := lp_sgn(p2);

    return lp_multSignature(si, mj) . lp_multSignature(sj, mi)
  end;

% compute Spoly of f and g, where f and g are tuples in internal representation
asserted procedure lp_spoly(
                    f: LabeledPolynomial,
                    g: LabeledPolynomial): LabeledPolynomial;
  begin scalar p1, p2, e1, e2, c1, c2, elcm, mult1, mult2, ans,
                sgn1, sgn2, sgn11, sgn22, anssgn;

    % here, signatures do not cancel each other,
    % Maybe write an assert?..
    % so the updated signature can be computed easily

    p1 := lp_eval(f);
    p2 := lp_eval(g);
    e1 := poly_leadExp(p1);
    e2 := poly_leadExp(p2);

    elcm := poly_lcmExp(e1, e2);

    mult1 := poly_subExp(elcm, e2);
    mult2 := poly_subExp(elcm, e1);

    ans := poly_paircomb(p1, mult2, poly_leadCoeff(p1), p2, mult1, poly_leadCoeff(p2));

    sgn1 := lp_sgn(f);
    sgn2 := lp_sgn(g);

    sgn11 := lp_multSignature(sgn1, mult2);
    sgn22 := lp_multSignature(sgn2, mult1);

    anssgn := lp_potMaxSignature(sgn11, sgn22);

    return lp_LabeledPolynomial2(ans, anssgn)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LabeledPolynomial sorting ad-hoc

% return true if lead(poly1) < lead(poly2)
asserted procedure lp_cmpLPLead(lp1, lp2);
  poly_cmpPolyLead(lp_eval(lp1), lp_eval(lp2));

% return true if lead(poly2) < lead(poly1)
asserted procedure lp_cmpLPLeadRev(lp1, lp2);
  poly_cmpPolyLead(lp_eval(lp2), lp_eval(lp1));

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
