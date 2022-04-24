
module f5core;

% Core module that contains f5 algorithm implementation.
%
% BasisKeeper is a structure that stores all emerging polynomials in a vector.
% More precicely, the structure of BasisKeeper is the following
%   {vector, filled, capacity}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure core_BasisKeeper(capacity);
  begin scalar v;
        integer filled;
    filled := 0;
    v := mkvect(capacity);
    return {v, filled, capacity}
  end;

asserted procedure core_addPoly(r, f);
  begin;
    putv(car r, cadr r, f);
    cadr r := cadr r #+ 1
  end;

asserted procedure core_setPoly(r, i, f);
  putv(car r, i, f);

asserted procedure core_getPoly(r, i);
  getv(car r, i);

asserted procedure core_getBasisIdx(r);
  (cadr r) #- 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure core_pairTotalDegreeCmp(p1, p2);
  poly_totalDeg(core_getPairLcm(p1)) #< poly_totalDeg(core_getPairLcm(p2));

asserted procedure core_pairLcmCmp(p1, p2);
  poly_cmpExp(core_getPairLcm(p1), core_getPairLcm(p2));

asserted procedure core_assocSgnCmp(pr1, pr2);
  lp_sgnCmp(lp_sgn(cdr pr1), lp_sgn(cdr pr2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure core_getPairLcm(p);
  car p;

asserted procedure core_getPairFirst(p);
  (cadr p) . (caddr p);

asserted procedure core_getPairSecond(p);
  (cadddr p) . (car cddddr p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the given set of polynomials `inputBasis` construct
% the standard basis of unit vectors corresponding to F
% , i.e.,
% {f, g}  --> {(1, 0), (0, 1)}
asserted procedure core_constructModule(inputBasis: List): List;
  begin scalar poly, elem, outputModule;
        integer i;

    inputBasis := sort(inputBasis, 'poly_leadTotalDegreeCmp);

    i := 0;
    while inputBasis do <<
      poly . inputBasis := inputBasis;
      elem := lp_LabeledPolynomial1(poly, i);
      outputModule := elem . outputModule;
      i := i + 1
    >>;
    return reversip(outputModule)
  end;

asserted procedure core_findReductor(k, Gprev, newGcurr, r, Rule);
  begin scalar tt, rk, rkev, rksgn, rj, ri, u, rjsgn, urjsgn, tj;
        integer j, ans;

    ans := 0;
    rk := core_getPoly(r, k);
    rkev := lp_eval(rk);
    rksgn := lp_sgn(rk);

    tt := poly_leadExp(lp_eval(rk));
    for each j in newGcurr do <<
      rj := core_getPoly(r, j);
      tj := poly_leadExp(lp_eval(rj));
      if (ans = 0) and poly_divExp!?(tj, tt) then <<
        u := poly_subExp(tt, tj);
        rjsgn := lp_sgn(rj);
        urjsgn := lp_multSgn(rjsgn, u);
        if not lp_eqSignature(urjsgn, rksgn) then
          if not core_isRewritable(u, j, r, Rule) then
            if not core_isTopReducibleMonom(lp_sgnMonom(urjsgn), Gprev, r) then
              ans := j
      >>
    >>;

    return ans
  end;

asserted procedure core_topReductionF5(k, Gprev, newGcurr, r, Rule);
  begin scalar rk, p, q, u, newr, flag, pev;
        integer j;

    rk := core_getPoly(r, k);
    if lp_iszero!?(rk) then <<
      prin2t "Reduction to zero!";
      return nil . nil
    >>;

    j := core_findReductor(k, Gprev, newGcurr, r, Rule);
    % COEFF
    if j #= 0 then <<
      core_setPoly(r, k, lp_normalize(rk));
      return k . nil
    >>;

    p := rk;
    q := core_getPoly(r, j);
    u := poly_subExp(poly_leadExp(lp_eval(p)), poly_leadExp(lp_eval(q)));
    % c := mod_div(poly_leadCoeff(lp_eval(p)), poly_leadCoeff(lp_eval(q)));

    flag . pev := core_reducePolyByTop(lp_eval(p), lp_eval(q));

    lp_setEval(p, pev);
    core_setPoly(r, k, p);

    % COEFF
    if not poly_iszero!?(pev) then
      p := lp_normalize(p);

    return if lp_sgnCmp(lp_multSgn(lp_sgn(q), u), lp_sgn(p)) then
      nil . k . nil
    else <<
      newr := lp_LabeledPolynomial2(lp_eval(p), lp_multSgn(lp_sgn(q), u));
      core_addPoly(r, newr);
      core_addRule(Rule, core_getBasisIdx(r), lp_sgn(newr), r);
      nil . k . core_getBasisIdx(r) . nil
    >>
  end;

asserted procedure core_reducePolyByTop(f, g);
  begin scalar newf, updated, fexps, fcoeffs,
                glead, gc, fc, flag, fex, fmult, gmult;
    newf := f;
    updated := nil;

    fe := poly_leadExp(f);
    fc := poly_leadCoeff(f);

    ge := poly_leadExp(g);
    gc := poly_leadCoeff(g);

    flag := poly_divExp!?(ge, fe);

    if flag then <<
      fmult := poly_zeroExp();
      gmult := poly_subExp(fe, ge);
      % f*mult - g*mult
      % COEFF
      newf := poly_paircomb(f, fmult, fc, g, gmult, gc);
      updated := t;
    >>;

    return updated . newf
  end;

asserted procedure core_reducePolyBy(f, g);
  begin scalar newf, updated, fexps, fcoeffs,
                glead, gc, fc, flag, fex, fmult, gmult;
    newf := f;

    updated := nil;

    fexps := poly_getExps(f);
    fcoeffs := poly_getCoeffs(f);
    glead := poly_leadExp(g);
    gc := poly_leadCoeff(g);

    while (not updated) and fexps do <<
      fex := pop(fexps);
      fc  := pop(fcoeffs);
      flag := poly_divExp!?(glead, fex);
      if flag then <<
        fmult := poly_zeroExp();
        gmult := poly_subExp(fex, glead);
        % f*mult - g*mult
        % COEFF
        newf := poly_paircomb(f, fmult, fc, g, gmult, gc);
        updated := t;
      >>;
    >>;

    return updated . newf
  end;

asserted procedure core_normalForm(f, Gprev, r);
  begin scalar newf, updated, reducer, modified,
                updatedToreturn;
    newf := f;
    updatedToreturn := nil;

  start:
    updated := nil;
    for each g in Gprev do <<
      reducer := lp_eval(core_getPoly(r, g));
      if not poly_iszero!?(reducer) and not poly_iszero!?(newf) then <<
        modified . newf := core_reducePolyBy(newf, reducer);
        updated := updated or modified;
        updatedToreturn := updated or updatedToreturn
      >>
    >>;

    if updated then
      goto start;

    return updatedToreturn . newf
  end;

asserted procedure core_normalFormTop(f, Gprev, r);
  begin scalar newf, updated, reducer, modified,
                updatedToreturn;
    newf := f;
    updatedToreturn := nil;

  start:
    updated := nil;
    for each g in Gprev do <<
      reducer := lp_eval(core_getPoly(r, g));
      if not poly_iszero!?(reducer) and not poly_iszero!?(newf) then <<
        modified . newf := core_reducePolyByTop(newf, reducer);
        updated := updated or modified;
        updatedToreturn := updated or updatedToreturn
      >>
    >>;

    if updated then
      goto start;

    return updatedToreturn . newf
  end;

asserted procedure core_insertSorted(todo, j, r);
  begin scalar tmp, prev, sj, s;
    if null todo then
      return j . nil;

    sj  := lp_sgn(core_getPoly(r, j));
    s   := lp_sgn(core_getPoly(r, car todo));
    if lp_sgnCmp(sj, s) then
      return j . todo;

    tmp := todo;
    while cdr todo and lp_sgnCmp(sj, lp_sgn(core_getPoly(r, cadr todo))) do <<
      todo := cdr todo
    >>;

    cdr todo := j . cdr todo;

    return tmp
  end;

asserted procedure core_reduction(S, Gprev, Gcurr, r, Rule);
  begin scalar todo, S, completed, newGcurr, rk, rknfeval,
                newcompleted, redo, flag;
        integer k, j;

    todo := S;
    completed := nil;

    newGcurr := copy(Gcurr);

    while todo do <<
      k := pop(todo);
      rk := core_getPoly(r, k);

      flag . rknfeval := if !*f5fullreduce then
        core_normalForm(lp_eval(rk), Gprev, r)
      else
        core_normalFormTop(lp_eval(rk), Gprev, r);

      lp_setEval(rk, rknfeval);
      core_setPoly(r, k, rk);

      newcompleted . redo := core_topReductionF5(k, Gprev, newGcurr, r, Rule);

      if newcompleted then <<
        completed := newcompleted . completed;
        newGcurr := newcompleted . newGcurr
      >>;

      for each j in redo do
        todo := core_insertSorted(todo, j, r);
    >>;

    return completed
  end;

asserted procedure core_isTopReducible(f, Gprev, r);
  begin scalar tf, g, gi, ans;
    ans := nil;
    tf := poly_leadExp(f);
    while Gprev and (not ans)  do <<
      gi := pop(Gprev);
      g  := core_getPoly(r, gi);
      if poly_divExp!?(poly_leadExp(lp_eval(g)), tf) then
        ans := t
    >>;
    return ans;
  end;

asserted procedure core_isTopReducibleMonom(m, Gprev, r);
  begin scalar tf, g, gi, ans;
    ans := nil;
    while Gprev and (not ans) do <<
      gi := pop(Gprev);
      g  := core_getPoly(r, gi);
      if poly_divExp!?(poly_leadExp(lp_eval(g)), m) then
        ans := t
    >>;
    return ans;
  end;

asserted procedure core_criticalPair(i, k, l, Gprev, r);
  begin scalar rk, rl, tk, tl, tt, u1, u2, sgn1, sgn2,
                usgn1, usgn2;
        integer i;

    rk := core_getPoly(r, k);
    rl := core_getPoly(r, l);

    tk := poly_leadExp(lp_eval(rk));
    tl := poly_leadExp(lp_eval(rl));

    tt := poly_lcmExp(tk, tl);

    u1 := poly_subExp(tt, tk);
    u2 := poly_subExp(tt, tl);

    sgn1 := lp_sgn(rk);
    sgn2 := lp_sgn(rl);

    usgn1 := lp_multSgn(sgn1, u1);
    usgn2 := lp_multSgn(sgn2, u2);

    if (lp_sgnIndex(sgn1) #= i) and core_isTopReducibleMonom(lp_sgnMonom(usgn1), Gprev, r) then
      return nil;

    if (lp_sgnIndex(sgn2) #= i) and core_isTopReducibleMonom(lp_sgnMonom(usgn2), Gprev, r) then
      return nil;

    if lp_sgnCmp(usgn1, usgn2) then <<
      u1 . u2 := u2 . u1;
      k . l := l . k
    >>;

    return {tt, k, u1, l, u2}
  end;

asserted procedure core_RewriteRule(idx, monom);
  {idx, monom};

asserted procedure core_ruleMonom(rr);
  cadr rr;

asserted procedure core_ruleIndex(rr);
  car rr;

asserted procedure core_addRule(Rule, k, sgn, r);
  begin scalar s, srule;
    s := lp_sgn(core_getPoly(r, k));
    srule := getv(Rule, lp_sgnIndex(s));
    srule := core_RewriteRule(k, lp_sgnMonom(s)) . srule;
    putv(Rule, lp_sgnIndex(s), srule)
  end;

asserted procedure core_findRewriting(u, k, r, Rule);
  begin scalar sgn, flag, krules, krule;
        integer ans;

    sgn := lp_sgn(core_getPoly(r, k));
    krules := getv(Rule, lp_sgnIndex(sgn));
    ans := k;

    % be careful with order in krules
    while krules and (not flag) do <<
      krule := pop(krules);
      flag := poly_divExp!?(core_ruleMonom(krule), lp_sgnMonom(lp_multSgn(sgn, u)));
      if flag then
        ans := core_ruleIndex(krule)
    >>;

    return ans
  end;

asserted procedure core_isRewritable(u, k, r, Rule);
  not (core_findRewriting(u, k, r, Rule) #= k);

asserted procedure core_spoly(f, g);
  begin scalar e1, e2, elcm, mult1, mult2;
    e1 := poly_leadExp(f);
    e2 := poly_leadExp(g);

    elcm := poly_lcmExp(e1, e2);

    mult1 := poly_subExp(elcm, e2);
    mult2 := poly_subExp(elcm, e1);

    return poly_paircomb(f, mult2, poly_leadCoeff(f), g, mult1, poly_leadCoeff(g))
  end;

asserted procedure core_computeSpolys(pairs, r, Rule);
  begin scalar S, pairs, p, u, v, lpk, lpl, alS,
                evals, sgns, lpnew, localSgnCmp;
        integer l, k;

    S := nil;
    pairs := sort(pairs, 'core_pairLcmCmp);

    while pairs do <<
      p := pop(pairs);
      k . u := core_getPairFirst(p);
      l . v := core_getPairSecond(p);
      if (not core_isRewritable(u, k, r, Rule)) and (not core_isRewritable(v, l, r, Rule)) then <<
        lpk := core_getPoly(r, k);
        lpl := core_getPoly(r, l);

        evals := core_spoly(lp_eval(lpk), lp_eval(lpl));
        sgns := lp_multSgn(lp_sgn(lpk), u);

        lpnew := lp_LabeledPolynomial2(evals, sgns);

        core_addPoly(r, lpnew);

        core_addRule(Rule, core_getBasisIdx(r), sgns, r);

        if not poly_iszero!?(evals) then
          S := core_getBasisIdx(r) . S
      >>
    >>;

    % Sort indices S w.r.t. values in r
    alS := for each i in S collect i . core_getPoly(r, i);
    alS := sort(alS, 'core_assocSgnCmp);

    S := for each pr in alS collect car pr;

    return S
  end;

asserted procedure core_incrementalBasis(i, Gprev, r, Rule);
  begin scalar Gcurr, pairs, p, S, dpairs, reduced, k, tmp;
        integer i, j, d, currIdx;

    currIdx := core_getBasisIdx(r);

    Gcurr := copy(Gprev);
    Gcurr := currIdx . Gcurr;

    pairs = nil;
    for each j in Gprev do <<
      p := core_criticalPair(i, currIdx, j, Gprev, r);
      if p then
        pairs := p . pairs
    >>;

    while pairs do <<
      pairs := sort(pairs, 'core_pairTotalDegreeCmp);

      p := pop(pairs);
      d := poly_totalDeg(core_getPairLcm(p));
      dpairs := p . nil;
      while pairs and (poly_totalDeg(core_getPairLcm(car pairs)) #= d) do <<
        dpairs := p . dpairs;
        p := pop(pairs)
      >>;

      S := core_computeSpolys(dpairs, r, Rule);

      reduced := core_reduction(S, Gprev, Gcurr, r, Rule);

      % prin2t r;

      while reduced do <<
        k := pop(reduced);
        tmp := Gcurr;
        while tmp do <<
          j := pop(tmp);
          p := core_criticalPair(i, j, k, Gprev, r);
          if p then
            pairs := p . pairs
        >>;
        Gcurr := k . Gcurr
      >>
    >>;

    return Gcurr
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure core_getReducers(i, G);
  begin scalar reducers, poly;
        integer j;
    j := 1;
    while G do <<
      poly := pop(G);
      if not (j #= i) then
        reducers := poly . reducers;
      j := j #+ 1
    >>;
    return reducers
  end;

% Given the set of generators `basis` apply the autoreduction algorithm
% to interreduce generators w.r.t. each other
asserted procedure core_interreduceBasis(Gprev, r);
  begin scalar updated, Gfull, newGprev, f, reducers, p;
        integer m, idx, full;

  full := nil;

  start:
    m := length(Gprev);
    Gfull := Gprev;
    newGprev := nil;
    for i := 1:m do <<
      idx := pop(Gprev);
      f := core_getPoly(r, idx);
      reducers := core_getReducers(i, Gfull);

      updated . p := if full then
        core_normalForm(lp_eval(f), reducers, r)
      else
        core_normalFormTop(lp_eval(f), reducers, r);

      lp_setEval(f, p);
      core_setPoly(r, idx, f);

      if not poly_iszero!?(p) then
        newGprev := idx . newGprev
    >>;
    Gprev := newGprev;

    if not full and updated then
      goto start;

    if not full and !*f5fullreduce then <<
      full := t;
      goto start
    >>;

    return Gprev
  end;

% normalize each generator in the given `basis` by the leading coefficient
asserted procedure core_normalizeBasis(basis: List): List;
  for each x in basis collect lp_normalize(x);

% given a Groebner `basis`
% returns the unique Groebner basis of the corresponding ideal <basis>
%
% Output invariants:
%  . the basis is interreduced
%  . the basis is normalized by leading coefficients
%  . the basis is sorted by leading terms increasingly
asserted procedure core_standardizeOutput(basis: List): List;
  begin;
    basis := core_normalizeBasis(basis);
    return sort(basis, 'lp_cmpLPLeadReverse)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The heart of the package, --
% the function to take the list of ideal generators `basis`
% and to return a Groebner basis of this idel
asserted procedure core_groebner1(basis: List): List;
  begin scalar f1, r, Gprev, Rule, fi;
        integer m, i, initialBasisSize;

    m  := length(basis);
    f1 := pop(basis);

    % COEFF
    f1 := lp_normalize(f1);

    initialBasisSize := 1000;
    r := core_BasisKeeper(initialBasisSize);
    core_addPoly(r, f1);

    Gprev := {0};
    i := 1;

    Rule := mkvect(m);

    % % WARN: added !
    % % core_addRule(Rule, 0, lp_sgn(f1), r);

    while i #< m do <<
      fi := pop(basis);

      % COEFF
      fi := lp_normalize(fi);
      core_addPoly(r, fi);
      putv(Rule, i, nil);

      % % WARN: added !
      % % core_addRule(Rule, core_getBasisIdx(), lp_sgn(fi), r);

      Gprev := core_incrementalBasis(i, Gprev, r, Rule);

      % prin2t {i, getv(Rule, i)};

      i := i #+ 1
    >>;

    if !*f5fullreduce then
      Gprev := core_interreduceBasis(Gprev, r);

    basis := nil;
    while Gprev do <<
      i := pop(Gprev);
      basis := core_getPoly(r, i) . basis
    >>;

    basis := core_standardizeOutput(basis);

    return basis
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if <polys> contains in <basis>,
% assuming `basis` is a Groebner basis
asserted procedure core_checkIdealInclusion1(basis: List, polys: List);
  begin;
    return t
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if `basis` is a Groebner basis by checking all S-polynomials
asserted procedure core_isGroebner1(basis: List);
  begin;
    return t
  end;


endmodule;

% trst core_criticalPair;
% trst core_groebner1;
% trst core_incrementalBasis;
% trst core_computeSpolys;
% trst core_addRule;
% trst core_reduction;
% trst core_constructModule;
% trst core_isTopReducibleMonom;
% trst core_isRewritable;
% trst core_findRewriting;
% trst core_isTopReducible;
% trst core_findReductor;
% trst core_topReduction;
% trst core_interreduceBasis;
% trst core_normalForm2;
% trst core_insertSorted;

end; % of file
