module f5core;
% The main module with the f5 algorithm implementation.

% fluid '(NREDUCTIONSF5 NNORMALFORMS);
% NREDUCTIONSF5 := 0;
% NNORMALFORMS := 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CriticalPair

% The struct implements the usual Critical Pair interface:
% the pair stores information about two polynomials to construct
% an S-polynomial from it later.
% Assuming there exist two polynomials `f1` and `f2` and their S-polynomial is
%   u*f1 - v*f2
% Then CriticalPair objects are represented with a 6-item list
%   {'cp, tt, k, u, l, v}
%
% Where
%   tt is a Term for the lcm(leadTerm(f1), leadTerm(f2)),
%   k is an Integer for the index of `f1` in the current Basistracker,
%   u is a Term for the multiplier of `f1`,
%   l is an Integer for the index of `f2` in the current Basistracker,
%   v is a Term for the multiplier of `f2`
%
% It is safe to assume that (k, u) >= (l, v) as Signatures
%             TODO: or even (k, u) > (l, v)

asserted inline procedure core_CriticalPair(tt: Term, k: Integer, u: Term,
                                    l: Integer, v: Term): CriticalPair;
  {'cp, tt, k, u, l, v};

% Return tt
asserted inline procedure core_getPairLcm(p: CriticalPair): Term;
  car p;

% Return k . u1
asserted inline procedure core_getPairFirst(p: CriticalPair): DottedPair;
  (cadr p) . (caddr p);

% Return l . u2
asserted inline procedure core_getPairSecond(p: CriticalPair): DottedPair;
  (cadddr p) . (car cddddr p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RewriteRule

% The struct implements Rewrite Rule interface to keep track of signatures
% that were already produced during execution before.
%
% The RewriteRule object is a list
%   {'rr, index, term}
% Here, index is an Integer, an index of polynomial from the Basistracker structure,
% and term is a Term, a multiplier of some signature.

asserted inline procedure core_RewriteRule(idx: Integer, tt: Term): RewriteRule;
  {'rr, idx, tt};

% Return index
asserted inline procedure core_getRuleIndex(r: RewriteRule): Integer;
  car r;

% Return multiplier
asserted inline procedure core_getRuleTerm(r: RewriteRule): Term;
  cadr r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basistracker

% Implements the Basistracker struct. Basistracker instance is an object
% that stores all polynomials produced during f5 execution in a single vector.
% More precicely, the structure of Basistracker instance is represented as
%   {'bt, polys, filled, capacity}
% Where
%   polys is a Vector of `Polynomial`s,
%   filled is an Integer equal to the number of elements filled in `polys`,
%   capacity is an Integer equal to the length of `polys`
%
% It safe to assume that polynomials added to the current Basistracker
% are normalized.

fluid '(core_initialBasisSize!*);
% We take a big number from the start,
% so that there is no need to copy the storage vector very often
% TODO: scale this with the basis growing
core_initialBasisSize!* := 10;

asserted procedure core_Basistracker(capacity: Integer): Basistracker;
  {'bt, mkvect(capacity), 0, capacity};

% Adds LabeledPolynomial f to the basis
asserted procedure core_addPoly(r: Basistracker, f: LabeledPolynomial);
  <<
    putv(cadr r, caddr r, f);
    caddr r := caddr r #+ 1
  >>;

% Sets the i'th polynomial in the basis to f
asserted inline procedure core_setPoly(r: Basistracker, i: Integer, f: Polynomial);
  putv(cadr r, i, f);

% Returns the i'th polynomial from the basis
asserted inline procedure core_getPoly(r: Basistracker, i: Integer): LabeledPolynomial;
  getv(cadr r, i);

% Returns the index of the last polynomial added to the basis
asserted inline procedure core_getBasisIdx(r: Basistracker): Integer;
  (caddr r) #- 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional comparators for sorting

% TODO: add the current term order cmp as a tie-breaker
asserted procedure core_pairTotalDegreeCmp(p1: CriticalPair, p2: CriticalPair);
  poly_totalDegTerm(core_getPairLcm(p1)) #< poly_totalDegTerm(core_getPairLcm(p2));

% compare critical pairs by lcm term according to the current term order
asserted procedure core_pairLcmCmp(p1: CriticalPair, p2: CriticalPair);
  poly_cmpTerm(core_getPairLcm(p1), core_getPairLcm(p2));

% compare associative list elements by their signature
asserted procedure core_assocSgnCmp(pr1: DottedPair, pr2: DottedPair);
  lp_cmpSgn(lp_sgn(cdr pr1), lp_sgn(cdr pr2));

% compare associative list elements by their leading term
% in the current term order
asserted procedure core_assocLeadCmp(pr1: DottedPair, pr2: DottedPair);
  poly_leadTotalDegreeCmp(lp_eval(cdr pr1), lp_eval(cdr pr2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the given set of polynomials `inputBasis` first
%   . interreduce input polynomials
%   . sort them by the leading term increasing
% and then construct the standard basis of unit vectors corresponding to F, i.e.,
%   {f, g}  --> {(1, 0), (0, 1)}
% where each basis vector is represented as a LabeledPolynomial
asserted procedure core_constructModule(inputBasis: List): List;
  begin scalar outputModule;
        integer i;
    inputBasis := core_interreduceInput(inputBasis);
    inputBasis := sort(inputBasis, 'poly_leadTotalDegreeCmp);
    i := 0;
    while inputBasis do <<
      push(lp_LabeledPolynomial1(pop(inputBasis), i), outputModule);
      i := i + 1
    >>;
    return reversip(outputModule)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the normal form of f w.r.t. polynomials r_i where i ∈ Gprev
%
% Two cases are possible:
%  1. At least one reduction step happened,
%     then we return the true flag together with the result of all reductions.
%  2. No reductions where performed, then nil flag and f itself are returned
%
% During reductions, the topReduce flag controls desired reduction type.
% If topReduce is set, only top-reductions happen.
asserted procedure core_normalForm(f: Polynomial, Gprev: List,
                                      r: Basistracker, topReduce): DottedPair;
  begin scalar updated, reducer, modified, updatedToreturn;
    updated := t;
    while updated do <<
      updated := nil;
      for each g in Gprev do <<
        reducer := lp_eval(core_getPoly(r, g));
        if not poly_iszero!?(reducer) and not poly_iszero!?(f) then <<
          modified . f := if topReduce then
            poly_tryTopReductionStep(f, reducer)
          else
            poly_tryReductionStep(f, reducer);
          updated := updated or modified;
          updatedToreturn := updated or updatedToreturn
        >>
      >>
    >>;
    return updatedToreturn . newf
  end;

% Same as the above, but all possible reducers are
% already stored in `reducers` as polynomials.
asserted procedure core_normalFormReducers(f: Polynomial, reducers: List,
                                              topReduce): DottedPair;
  begin scalar updated, reducer, modified, updatedToreturn;
    updated := t;
    while updated do <<
      updated := nil;
      for each reducer in reducers do <<
        if not poly_iszero!?(reducer) and not poly_iszero!?(f) then <<
          modified . f := if topReduce then
            poly_tryTopReductionStep(f, reducer)
          else
            poly_tryReductionStep(f, reducer);
          updated := updated or modified;
          updatedToreturn := updated or updatedToreturn
        >>
      >>
    >>;
    return updatedToreturn . f
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return G without the i-th element, e.g.,
% for G = {8, 4, 5} and i = 2 return {8, 5}
asserted procedure core_getReducers(i: Integer, G: List): List;
  begin scalar reducers, poly;
        integer j;
    j := 1;
    while G do <<
      poly := pop(G);
      if not (j #= i) then
        push(poly, reducers);
      j := j #+ 1
    >>;
    return reducers
  end;

% Given the set of basis generators `input`
% apply passes of the Autoreduction algorithm
% until no further reductions happen
asserted procedure core_interreduceInput(input: List): List;
  begin scalar reducers, updated, reduced, p, f,
                newInput, i;
    updated := t;
    while updated do <<
      updated := nil;
      newInput := nil;
      while input do <<
        f := pop(input);
        reducers := append(newInput, input);
        reduced . p := core_normalFormReducers(f, reducers, nil);
        updated := updated or reduced;
        if not poly_iszero!?(p) then
          push(p, newInput)
      >>;
      input := newInput
    >>;
    return newInput
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part of the main F5 criterion:
%   checks if m is top reducible by the previous Groebner basis;
%   If it is, returns t;
%   Otherwise, return nil
asserted procedure core_isTopReducibleTerm(m: Term, Gprev: List,
                                            r: Basistracker);
  begin scalar tf, glead, gi, isReducible;
    while Gprev and (not isReducible) do <<
      gi := pop(Gprev);
      glead := oly_leadTerm(lp_eval(core_getPoly(r, gi)));
      if poly_divTerm!?(glead, m) then
        isReducible := t
    >>;
    return isReducible;
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add rewrite rule with index k and term term(sgn) to the set
% of rewrite rules at position index(sgn)
asserted procedure core_addRule(Rule: Vector, sgn: Signature, k: Integer);
  begin scalar ruleS;
    ruleS := getv(Rule, lp_indexSgn(sgn));
    push(core_RewriteRule(k, lp_termSgn(sgn)), ruleS);
    putv(Rule, lp_indexSgn(sgn), ruleS)
  end;

% Searches for the first rewriter for signature u*sgn(r_k)
% and returns the index of it.
% If there are no proper rewriters other than self, returns k
asserted procedure core_findRewriting(u: Term, k: Integer,
                                r: Basistracker, Rule: Vector): Integer;
  begin scalar sgn, usgnt, foundRewriter, rulesAtK, ruleK;
        integer rewriter;
    rewriter := k;
    sgn      := lp_sgn(core_getPoly(r, k));
    usgnt    := lp_termSgn(lp_mulSgn(sgn, u));
    rulesAtK := getv(Rule, lp_sgnIndex(sgn));
    % Be careful with order in rulesAtK:
    % we want to traverse rewrite rules at index k
    % from the last to the first added ones,
    % and terminate as soon as the first matching rule is found
    while rulesAtK and (not foundRewriter) do <<
      ruleK := core_getRuleTerm(pop(rulesAtK));
      foundRewriter := poly_divTerm!?(ruleK, usgnt);
      if foundRewriter then
        rewriter := core_getRuleIndex(ruleK)
    >>;
    return rewriter
  end;

% Main Rewrite Criterion: Is signature u*sgn(r_k) rewritable?
asserted inline procedure core_isRewritable(u: Term, k: Integer,
                                        r: Basistracker, Rule: Vector);
  not (core_findRewriting(u, k, r, Rule) #= k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
asserted procedure core_filterRedundant(Gprev: List, r: Basistracker): List;
  begin scalar Gnew, alG, Gsort, reducers, glead, g;
        integer gi;
    % Sort indices S w.r.t. leading terms in r
    alG := for each i in Gprev collect i . core_getPoly(r, i);
    alG := sort(alG, 'core_assocLeadCmp);
    Gsort := for each pr in alG collect car pr;
    while Gsort do <<
      gi := pop(Gsort);
      g := lp_eval(core_getPoly(r, gi));
      glead := poly_leadTerm(g);
      if not core_isTopReducibleTerm(glead, Gnew, r) then
        push(gi, Gnew)
    >>;
    return Gnew
  end;

% normalize each generator in the given `basis`
% by dividing it by the leading coefficient
asserted procedure core_normalizeBasis(basis: List): List;
  for each x in basis collect lp_normalize(x);

% given a Groebner `basis`
% returns the unique Groebner basis of the corresponding ideal <basis>
%
% Output invariants:
%  . the basis is interreduced (if f5fullreduce if off, only top-reductions are performed)
%  . the basis is normalized by leading coefficients
%  . the basis is sorted by leading terms increasingly
asserted procedure core_standardizeOutput(basis: List): List;
  begin scalar normalizedBasis;
    normalizedBasis := core_normalizeBasis(basis);
    return sort(normalizedBasis, 'lp_cmpLPLeadReverse)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if <polys> contains in <basis> as ideals,
% assuming `basis` is a Groebner basis
asserted procedure core_checkIdealInclusion1(basis: List, polys: List);
  begin scalar ans, tmp, evals, nf;
    ans := t;
    while ans and polys do <<
      p := pop(polys);
      evals := lp_eval(p);
      _ . nf := core_normalFormTopReducers(evals, basis, t);
      if not poly_iszero!?(nf) then
        ans := nil
    >>;
    return ans
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if `basis` is a Groebner basis by checking that
% all S-polynomials reduce to zero
asserted procedure core_isGroebner1(basis: List);
  begin scalar ans, tmp, evals, nf, flag;
    ans := t;
    while basis do <<
      tmp := cdr basis;
      while ans and tmp do <<
        evals := poly_spoly(lp_eval(car tmp), lp_eval(car basis));
        _ . nf := core_normalFormTopReducers(evals, basis, t);
        if not poly_iszero!?(nf) then
          ans := nil;
        tmp := cdr tmp
      >>;
      basis := cdr basis
    >>;
    return ans
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For polynomial at index k searches for a possible top-reducer from
% the previous Groebner basis.
% A possible reducer should fullfil several criteria:
% First, it should be a top-reducer of r_k in a usual sense:
%   the leading term of r_k must be divisible by the leading term of reducer
% Secondly, reducer should admit a proper signature:
%   . The signature must be not equal to the signature of r_k,
%   . Reducer should be not rewritable AND not top-reducible itself
%
% If such reducer is found, return its index into r;
% if no reducer was found, return 0
asserted procedure core_findReducer(k: Integer, Gprev: List, newGcurr: List,
                                  r: Basistracker, Rule: Vector): Integer;
  begin scalar tt, rk, rkEval, rkSgn, rj, ri, u, rjSgnMult, tj;
        integer j, reducer;
    rk     := core_getPoly(r, k);
    rkEval := lp_eval(rk);
    rkSgn  := lp_sgn(rk);
    tt := poly_leadTerm(rkEval);
    for each j in newGcurr do <<
      rj := core_getPoly(r, j);
      tj := poly_leadTerm(lp_eval(rj));
      % TODO: check order
      if (reducer #= 0) and poly_divTerm!?(tj, tt) then <<
        u := poly_divExp(tt, tj);
        rjSgnMult := lp_multSgn(lp_sgn(rj), u);
        if not lp_eqSig(rjSgnMult, rksgn) then
          if not core_isRewritable(u, j, r, Rule) then
            if not core_isTopReducibleTerm(lp_sgnTerm(rjSgnMult), Gprev, r) then
              reducer := j
      >>
    >>;
    return reducer
  end;

% Adapted from John Perry et al.
%   https://arxiv.org/abs/0906.2967
%
% The most tricky part in F5 reduction, the top-reduction function.
% If the LabeledPolynomial being examined has polynomial part 0, then there
% is no data left in that particular signed polynomial - an
% empty ordered pair is returned.
%
% Otherwise top-reduction calls upon another sub-subroutine core_findReducer.
% Essentially, if core_findReducer comes back 0, the current LabeledPolynomial
% is made monic and returned to reduction to be placed in completed.
%
% If a top-reduction is deemed possible, then there are two possible cases:
% either the reduction will increase the signature of polynomial or it won't.
% In the latter case, the signature is maintained, the polynomial part
% is top-reduced and the LabeledPolynomial is returned to reduction to be
% added back into todo.
%
% In the former case, however, the signature will change. This
% is marked by adding a new polynomial to r with appropriate
% signature based upon the reducer. A new rule is added and then
% both reduced and reducer polynomials and are sent back to
% reduction to be added back into todo. This is done because
% reduced has a different signature than reducer and reducer might
% still be reducible by another LabeledPolynomial.
asserted procedure core_topReductionF5(k: Integer, Gprev: List, newGcurr: List,
                      r: Basistracker, Rule: Vector): DottedPair;
  begin scalar rk, p, u, newpoly, flag, reduced, reducerSgn;
        integer j;
    rk := core_getPoly(r, k);
    % if reduction to zero happened in the normal form -
    % meaning the system is not regular
    if lp_iszero!?(rk) then <<
      % prin2t "Reduction to zero!";
      return nil . nil
    >>;
    j := core_findReducer(k, Gprev, newGcurr, r, Rule);
    % no top reducers found in the previous Groebner basis --
    % reduction is not possible
    if j #= 0 then <<
      core_setPoly(r, k, lp_normalize(rk));
      return k . nil
    >>;
    % reducer rj found, perform top-reduction
    rj := core_getPoly(r, j);
    _ . reduced := core_reducePolyByTop(lp_eval(rk), lp_eval(rj));
    if not poly_iszero!?(reduced) then
      reduced := poly_normalize(reduced);
    % we need to check that signature of reducer
    % is not greater than the signature of rj
    u := poly_divTerm(poly_leadTerm(lp_eval(rk)), poly_leadTerm(lp_eval(rj)));
    reducerSgn := lp_multSgn(lp_sgn(rj), u);
    return if lp_sgnCmp(reducerSgn, lp_sgn(rk)) then <<
      % signatures OK, reduction is successful.
      % We form new new polynomial at index k,
      % and return it to be reduced once again later
      lp_setEval(rk, reduced);
      core_setPoly(r, k, rk);
      nil . k . nil
    >> else <<
      % signature of reducer is greater, reduction failed.
      % BUT, we still can form a new polynomial with signature of reducer
      % for further reduction
      newpoly := lp_LabeledPolynomial2(reduced, reducerSgn);
      core_addPoly(r, newpoly);
      core_addRule(Rule, lp_sgn(newpoly), core_getBasisIdx(r));
      % TODO: do we really want to return k here?
      nil . k . core_getBasisIdx(r) . nil
    >>
  end;

% insert index j in todo, so that
% indices in todo are sorted increasingly by signatures from r
asserted procedure core_insertSorted(todo: List, j: Integer,
                                        r: Basistracker): List;
  begin scalar tmp, prev, sj, s;
    if null todo then
      return {j};
    sj := lp_sgn(core_getPoly(r, j));
    s  := lp_sgn(core_getPoly(r, car todo));
    if lp_sgnCmp(sj, s) then
      return j . todo;
    tmp := todo;
    while cdr todo and lp_sgnCmp(sj, lp_sgn(core_getPoly(r, cadr todo))) do <<
      todo := cdr todo
    >>;
    cdr todo := j . cdr todo;
    return tmp
  end;

% The main reduction function.
% Given indexes of S-polynomials in S,
% computes the F5-style reducted form for each of them.
%
% By F5-style reduced form we mean a normal form,
% with respect to the previous Groebner basis,
% combined with a top-reduction that does not change the signature of LP.
asserted procedure core_reduction(S: List, Gprev: List, Gcurr: List,
                                    r: Basistracker, Rule: Vector): List;
  begin scalar todo, S, completed, newGcurr, rk, rknfeval,
                newcompleted, redo;
        integer k, j;
    todo := S;
    newGcurr := copy(Gcurr);
    while todo do <<
      k := pop(todo);
      rk := core_getPoly(r, k);
      % Compute the normal form;
      % If full reduction is not needed, compute only top normal form.
      % Otherwise, also reduce the polynomial tail.
      _ . rknfeval := if !*f5fullreduce then
        core_normalForm(lp_eval(rk), Gprev, r, t)
      else
        core_normalForm(lp_eval(rk), Gprev, r, t);
      lp_setEval(rk, rknfeval);
      core_setPoly(r, k, rk);
      % Compute the F5-style top-reduction;
      % If the polynomial
      newcompleted . redo := core_topReductionF5(k, Gprev, newGcurr, r, Rule);
      if newcompleted then <<
        push(newcompleted, completed);
        push(newcompleted, newGcurr)
      >>;
      for each j in redo do
        todo := core_insertSorted(todo, j, r);
    >>;
    return completed
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adapted from John Perry et al.
%   https://arxiv.org/abs/0906.2967
%
% Though at first glance this subroutine may look complicated,
%  core_computeSpolys essentially does one thing: form the new
%  S-polynomials output from critical_pairs as LabeledPolynomials.
%  We note that, because core_makeCriticalPair ensured
%  that sgn(u*k) < sgn(v*l), we know that the signature of all new
%  polynomials will always be of the form u_k*sgn(r_k) in core_computeSpolys.
asserted procedure core_computeSpolys(pairs: List, r: Basistracker,
                                        Rule: Vector): List;
  begin scalar S, p, u, v, lpk, lpl, alS, seval, ssgn;
        integer l, k;
    pairs := sort(pairs, 'core_pairLcmCmp);
    while pairs do <<
      p := pop(pairs);
      k . u := core_getPairFirst(p);
      l . v := core_getPairSecond(p);
      % Rewritten criterion
      if (not core_isRewritable(u, k, r, Rule)) and (not core_isRewritable(v, l, r, Rule)) then <<
        lpk := core_getPoly(r, k);
        lpl := core_getPoly(r, l);
        seval := poly_spoly(lp_eval(lpk), lp_eval(lpl));
        ssgn  := lp_multSgn(lp_sgn(lpk), u);
        core_addPoly(r, lp_LabeledPolynomial2(evals, sgns));
        core_addRule(Rule, sgns, core_getBasisIdx(r));
        if not poly_iszero!?(seval) then
          push(core_getBasisIdx(r), S)
      >>
    >>;
    % Sort indices S w.r.t. values in r
    alS := for each i in S collect i . core_getPoly(r, i);
    alS := sort(alS, 'core_assocSgnCmp);
    S := for each pr in alS collect car pr;
    return S
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constructs a CriticalPair from polynomials at indices k and l.
% If the pair is redundant by the F5 or the rewritten criterion, returns nil.
%
asserted procedure core_makeCriticalPair(i: Integer, k: Integer,
                        l: Integer, Gprev: List, r: Basistracker): CriticalPair;
  begin scalar rk, rl, tk, tl, tt, u1, u2, sgn1, sgn2, usgn1, usgn2;
        integer i;
    rk := core_getPoly(r, k);
    rl := core_getPoly(r, l);
    tk := poly_leadTerm(lp_eval(rk));
    tl := poly_leadTerm(lp_eval(rl));
    % tt = lcm(leadTerm(rk), leadTerm(rl))
    tt := poly_lcmTerm(tk, tl);
    % so that u1*rk - u2*rl is S-polynomial
    u1 := poly_subExp(tt, tk);
    u2 := poly_subExp(tt, tl);
    % signatures of u1*rk and u2*rl
    usgn1 := lp_mulSgn(lp_sgn(rk), u1);
    usgn2 := lp_mulSgn(lp_sgn(rl), u2);
    % F5 criterion
    if (lp_indexSgn(sgn1) #= i) and core_isTopReducibleTerm(lp_termSgn(usgn1), Gprev, r) then
      return nil;
    if (lp_sgnIndex(sgn2) #= i) and core_isTopReducibleTerm(lp_termSgn(usgn2), Gprev, r) then
      return nil;
    % Rewritten criterion
    if core_isRewritable(u1, k, r, Rule) or core_isRewritable(u2, l, r, Rule) then
      return nil;
    % the pair should be ordered
    if lp_cmpSgn(usgn1, usgn2) then <<
      u1 . u2 := u2 . u1;
      k . l := l . k
    >>;
    return core_criticalPair(tt, k, u1, l, u2)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given that Gprev indexes the basis of {f1..fi-1} in `r`,
% construct the basis of {f1..fi} and return a list of new indices
asserted procedure core_incrementalBasis(i: Integer, Gprev: List,
                          r: Basistracker, Rule: Vector): List;
  begin scalar Gcurr, pairs, p, S, dpairs, reduced, k, tmp;
        integer i, j, d, currIdx;
    currIdx := core_getBasisIdx(r);
    Gcurr := copy(Gprev);
    push(currIdx, Gcurr);
    for each j in Gprev do <<
      p := core_makeCriticalPair(i, currIdx, j, Gprev, r);
      if p then
        push(p, pairs)
    >>;
    while pairs do <<
      pairs := sort(pairs, 'core_pairTotalDegreeCmp);
      p := pop(pairs);
      d := poly_totalDegTerm(core_getPairLcm(p));
      dpairs := {p};
      while pairs and (poly_totalDegTerm(core_getPairLcm(car pairs)) #= d) do <<
        push(pop(pairs), dpairs)
      >>;
      S := core_computeSpolys(dpairs, r, Rule);
      reduced := core_reduction(S, Gprev, Gcurr, r, Rule);
      reduced := reversip(reduced);
      while reduced do <<
        k := pop(reduced);
        tmp := Gcurr;
        while tmp do <<
          j := pop(tmp);
          p := core_makeCriticalPair(i, j, k, Gprev, r);
          if p then
            push(p, pairs)
        >>;
        push(k, Gcurr)
      >>
    >>;
    return Gcurr
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The heart of the package, --
% the function to take a list of LabeledPolynomials as ideal generators,
% and to return a standardized Groebner basis of this ideal as output
%
% Output contracts:
%   . the output basis is top-reduced, meaning the number of generators is minimal,
%   . the output basis is sorted by the leading term increasingly,
%   . the output basis contains only normalized polynomials,
%   . if f5fullreduce is on, then the output basis is tail-reduced
asserted procedure core_groebner1(basis: List): List;
  begin scalar f1, r, Gprev, Rule, fi;
        integer m, i;
    m  := length(basis);
    % f1 - first polynomial added to the basis
    f1 := pop(basis);
    f1 := lp_normalize(f1);
    r := core_Basistracker(core_initialBasisSize!*);
    core_addPoly(r, f1);
    % Gprev indexes generators of the current basis in the Basistracker `r`,
    % So, Gprev := {0} indexes the basis of {f1}
    Gprev := {0};
    % Vector of RewriteRules
    Rule := mkvect(m);
    i := 1;
    % Incremental construction from index 1 to index m
    while i #< m do <<
      fi := pop(basis);
      fi := lp_normalize(fi);
      core_addPoly(r, fi);
      putv(Rule, i, nil);
      % construct the basis for {f1...fi} using the basis for {f1...fi-1}
      Gprev := core_incrementalBasis(i, Gprev, r, Rule);
      i := i #+ 1
    >>;
    % filter redundant generators
    Gprev := core_filterRedundant(Gprev, r);
    % if !*f5fullreduce then
    %   Gprev := core_interreduceBasis(Gprev, r);
    basis := for each i in Gprev collect core_getPoly(r, i);
    return core_standardizeOutput(basis)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trst core_groebner1;
trst core_constructModule;

endmodule;

end; % of file
