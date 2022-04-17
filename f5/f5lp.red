
module f5lp;

% The lp module provides Labeled Polynomial interface --
% a tuple of polynomials for usage in f5-style algorithms.
% LabeledPolynomial `p` is stored internally as a 4-item list:
%     {'lp, evaluation of `p`, signature index of `p`, signature monomial of `p`}
%
% The interface provides, in particular,
% functions lp_evaluation and lp_signature, that return
% second and {third, fouth} items of internal list respectively

% instantiate LabeledPolynomial from Polynomial and leading index
asserted inline procedure lp_LabeledPolynomial1(
                              poly: Polynomial,
                              idx: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly, {idx, poly_zeroExp()});

% instantiate LabeledPolynomial from Polynomial and signature
asserted inline procedure lp_LabeledPolynomial2(
                              poly: Polynomial,
                              sgn: List): LabeledPolynomial;
  'lbl . poly . sgn;

asserted inline procedure lp_getPoly(lp: LabeledPolynomial): Polynomial;
  cadr lp;

asserted inline procedure lp_getLabel(lp: LabeledPolynomial): List;
  cddr lp;

% if lp is zero (represented as ('lp, (p, nil, nil), smth, smth) )
% Same as lp_isSyzygy, in fact
asserted inline procedure lp_iszero!?(lp);
  poly_iszero!?(lp_evaluation(lp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return the signature of module element f
asserted inline procedure lp_signature(f: LabeledPolynomial): List;
  lp_getLabel(f);

% return the index of signature of module element f
asserted inline procedure lp_index(f: LabeledPolynomial): Integer;
  car lp_getLabel(f);

% return the evaluation of module element f
asserted inline procedure lp_evaluation(f: LabeledPolynomial): Polynomial;
   lp_getPoly(f);

% multiply signature sgn by monomial exponent vector ev
asserted inline procedure lp_multSignature(sgn: List, ev: List);
  begin scalar idx, term;
    idx . term := sgn;
    term := car term;
    return {idx, poly_sumExp(term, ev)}
  end;

% compare signatures sgn1 vs. sgn2 with the
% position over term extension
asserted procedure lp_potCmpSignature(sgn1: List, sgn2: List): Boolean;
  begin integer idx1, idx2;
        scalar ev1, ev2, sgn1, sgn2;
    idx1 . ev1 := sgn1;
    idx2 . ev2 := sgn2;
    return if idx1 equal idx2 then poly_cmpExp(ev1, ev2)
      else (idx1 > idx2)
  end;

% minimal signature based on pot comparison strategy
asserted inline procedure lp_potMinSignature(sgn1: List, sgn2: List): List;
  if lp_potCmpSignature(sgn1, sgn2) then sgn1 else sgn2;

% maximal signature based on pot comparison strategy
asserted inline procedure lp_potMaxSignature(sgn1: List, sgn2: List): List;
  if lp_potCmpSignature(sgn1, sgn2) then sgn2 else sgn1;

% compare LabeledPolynomials as module elems lp1 vs. lp2 with
% position over term extension
asserted procedure lp_potCmpLP(lp1: LabeledPolynomial, lp2: LabeledPolynomial): Boolean;
  begin scalar sgn1, sgn2;
    sgn1 := lp_signature(lp1);
    sgn2 := lp_signature(lp2);
    return lp_potCmpSignature(sgn1, sgn2)
  end;

% checks if lp represents a syzygy
asserted inline procedure lp_isSyzygy(lp: LabeledPolynomial): Boolean;
  poly_iszero!?(lp_evaluation(lp));

% checks if signature sgn1 divides signature sgn2
asserted procedure lp_signatureDivides(sgn1, sgn2): Boolean;
  (car sgn1 equal car sgn2) and poly_divExp!?(cadr sgn1, cadr sgn2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct a principal syzygy of LabeledPolynomials (fi, fj)
% as a LabeledPolynomial itself.
% Assumes i != j.
% A syzygy of form ej*fi - ei*fj, with evaluation equal to fj*fi - fi*fj == 0
% and signature = { min(j, i), lt(fi) if j < i else lt(fj) };
asserted procedure lp_principalSyzygy(fi: LabeledPolynomial,
                                      fj: LabeledPolynomial): LabeledPolynomial;
  begin scalar i, j, ei, ej, sgn;
      i := lp_index(fi);
      j := lp_index(fj);
      ei := poly_leadExp(lp_evaluation(fi));
      ej := poly_leadExp(lp_evaluation(fj));
      sgn := if i < j then {i, ej} else {j, ej};

      % TODO: poly_zero ?
      return lp_LabeledPolynomial2(poly_zero(), sgn)
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficient manipulation

asserted procedure lp_normalize(f: LabeledPolynomial);
  lp_LabeledPolynomial2(poly_normalize(lp_evaluation(f)), lp_signature(f));

asserted procedure lp_scaleDenominators(f: LabeledPolynomial);
  lp_LabeledPolynomial2(poly_scaleDenominators(lp_evaluation(f)), lp_signature(f));

asserted procedure lp_reduceCoeffs(f: LabeledPolynomial, prime);
  lp_LabeledPolynomial2(poly_reduceCoeffs(lp_evaluation(f), prime), lp_signature(f));

asserted procedure lp_reconstructCoeffs(f: LabeledPolynomial, prime);
  lp_LabeledPolynomial2(poly_reconstructCoeffs(lp_evaluation(f), prime), lp_signature(f));

asserted procedure lp_crtCoeffs(polyaccum, polycomp, modulo, prime);
  begin underlyingpoly;
    underlyingpoly := poly_crtCoeffs(lp_evaluation(polyaccum), lp_evaluation(polycomp), modulo, prime);
    return lp_LabeledPolynomial2(underlyingpoly, lp_signature(polyaccum))
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try to reduce module elem f with module elem g
% returns success flag with the result of reduction
asserted procedure lp_tryReduce1(f: LabeledPolynomial,
                                 g: LabeledPolynomial,
                                 isregular);
  begin scalar fexps, fcoeffs, fe, fc, zeroexp, gexplead,
                gcoefflead, flag, fevl, gevl, fsgn, gsgn,
                ans, zeroexp, gmult;

    % note that isregular=nil breaks the signature invariant
    % (storage of signature may become garbage as well)
    fevl := lp_evaluation(f);
    gevl := lp_evaluation(g);
    fsgn := lp_signature(f);
    gsgn := lp_signature(g);

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

% compute the normal form of LP f w.r.t. basis and return it
% If isregular is set, performs only regular reductions,
% otherwise the reduction may corrupt the signature of f
%
% Contract: use isregular=nil only in the final interreduction
%           when corrupted signature is not a problem
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

asserted procedure lp_eqList(l1, l2);
  if null l1 then
    t
  else
    car l1 #= car l2 and lp_eqList(cdr l1, cdr l2);

asserted procedure lp_eqSignature(sgn1, sgn2);
  car sgn1 #= car sgn2 and lp_eqList(cdr sgn1, cdr sgn2);

% Checks if f is is singularly-top-reducible by some element of G
asserted procedure lp_isSingularlyTopReducible(f: LabeledPolynomial, basis: List): Boolean;
  begin scalar g, fevl, gevl, ft, gt, flead, glead, m, ans;
    fevl  := lp_evaluation(f);
    ft    := lp_signature(f);
    flead := poly_leadExp(fevl);
    for each g in basis do <<
      gevl  := lp_evaluation(g);
      gt    := lp_signature(g);
      glead := poly_leadExp(gevl);
      if poly_divExp!?(glead, flead) then <<
        m := poly_subExp(flead, glead);
        if lp_eqSignature(ft, lp_multSignature(gt, m)) then
          ans := t
      >>
    >>;
    return ans
  end;

% compute cofactors of Spoly of p1 and p2
asserted procedure lp_spolyCofactors(f: LabeledPolynomial, g: LabeledPolynomial);
  begin scalar p1, p2, e1, e2, elcm, mult1, mult2;
    p1 := lp_evaluation(f);
    p2 := lp_evaluation(g);
    e1 := poly_leadExp(p1);
    e2 := poly_leadExp(p2);

    elcm := poly_lcmExp(e1, e2);

    mult1 := poly_subExp(elcm, e2);
    mult2 := poly_subExp(elcm, e1);

    return mult1 . mult2
  end;

% compute signatures of p1 and p2
% after multiplication by cofactors of Spolynomial of p1 and p2
asserted procedure lp_spolyMultSignatures(p1: LabeledPolynomial, p2: LabeledPolynomial): List;
  begin scalar mi, mj, si, sj;
    mi . mj := lp_spolyCofactors(p1, p2);
    si := lp_signature(p1);
    sj := lp_signature(p2);

    return lp_multSignature(si, mj) . lp_multSignature(sj, mi)
  end;

% compute Spoly of f and g, where f and g are tuples in internal representation
asserted procedure lp_spoly(f: LabeledPolynomial, g: LabeledPolynomial): LabeledPolynomial;
  begin scalar p1, p2, e1, e2, c1, c2, elcm, mult1, mult2, ans,
                sgn1, sgn2, sgn11, sgn22, anssgn;

    % here, signatures do not cancel each other,
    % Maybe write an assert?..
    % so the updated signature can be computed easily

    p1 := lp_evaluation(f);
    p2 := lp_evaluation(g);
    e1 := poly_leadExp(p1);
    e2 := poly_leadExp(p2);

    elcm := poly_lcmExp(e1, e2);

    mult1 := poly_subExp(elcm, e2);
    mult2 := poly_subExp(elcm, e1);

    ans := poly_paircomb(p1, mult2, poly_leadCoeff(p1), p2, mult1, poly_leadCoeff(p2));

    sgn1 := lp_signature(f);
    sgn2 := lp_signature(g);

    sgn11 := lp_multSignature(sgn1, mult2);
    sgn22 := lp_multSignature(sgn2, mult1);

    anssgn := lp_potMaxSignature(sgn11, sgn22);

    return lp_LabeledPolynomial2(ans, anssgn)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LabeledPolynomial sorting ad-hoc

% return true if lead(poly1) < lead(poly2)
asserted procedure lp_cmpLPLead(lp1, lp2);
  poly_cmpPolyLead(lp_evaluation(lp1), lp_evaluation(lp2));

% return false if lead(poly1) < lead(poly2)
asserted procedure lp_cmpLPLeadRev(lp1, lp2);
  not poly_cmpPolyLead(lp_evaluation(lp1), lp_evaluation(lp2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trst lp_LabeledPolynomial;
% trst lp_spoly;
% trst lp_spolyMultSignatures;
% trst lp_spolyCofactors;
% trst lp_multSignature;

% trst lp_tryReduce1;
% trst lp_tryReduce;
% trst lp_normalForm;
% trst lp_isSingularlyTopReducible;

% trst lp_principalSyzygy;

endmodule;


end;  % of file
