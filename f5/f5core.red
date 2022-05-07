module f5core;
% The main module with the f5 algorithm implementation.

% fluid '(NREDUCTIONSF5 NNORMALFORMS);
% NREDUCTIONSF5 := 0;
% NNORMALFORMS := 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *** Notes on the implementation.
%
% In our implementation of the F5 algorithm all polynomials are wrapped
% in the LabeledPolynomial type (interface implemented in f5lp.red).
% LabeledPolynomial provides the `lp_eval` function to access the polynomial part itself.
% Further, we use the notions of "polynomial" and "labeled polynomial" interchangeably.
%
% Most of the functions below manipulate polynomials (represented as a LabeledPolynomial).
% Polynomials are stored in the `Basistracker` struct, which is usually passed
% as a function argument and has a name `r`. To get to these polynomials,
% integer indices are used. So, by saying "a polynomial at index k" we mean
% a polynomial stored in the Basistracker `r` at index k (our notation for this is `r_k`).
% To get the polynomial itself, `core_getPoly` is used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CriticalPair

% The struct implements the usual Critical Pair interface:
% the pair stores information about two polynomials to construct
% an S-polynomial from it later.
% Assuming there exist two polynomials `f1` and `f2` and their S-polynomial is
%   u*f1 - v*f2
% Then CriticalPair object is represented with a 6-item list
%   {'cp, tt, k, u, l, v}
% Where
%   tt is a Term for the lcm(leadTerm(f1), leadTerm(f2)),
%   k is an Integer for the index of `f1` in the current Basistracker (described below),
%   u is a Term for the multiplier of `f1`,
%   l is an Integer for the index of `f2` in the current Basistracker,
%   v is a Term for the multiplier of `f2`
%
% It is safe to assume that (k > l) or (k = l and u >= v).

asserted inline procedure core_CriticalPair(tt: Term, k: Integer, u: Term,
                                    l: Integer, v: Term): CriticalPair;
  {'cp, tt, k, u, l, v};

% Return tt
asserted inline procedure core_getPairLcm(p: CriticalPair): Term;
  cadr p;

% Return k . u1
asserted inline procedure core_getPairFirst(p: CriticalPair): DottedPair;
  (caddr p) . (cadddr p);

% Return l . u2
asserted inline procedure core_getPairSecond(p: CriticalPair): DottedPair;
  (car cddddr p) . (cadr cddddr p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RewriteRule

% The struct implements Rewrite Rule interface to keep track of signatures
% that have already been produced during the f5 execution.
%
% A RewriteRule object is a list
%   {'rr, index, term}
% Here, index is an Integer, an index of polynomial from the Basistracker structure,
% and term is a Term, a multiplier of some signature.

asserted inline procedure core_RewriteRule(idx: Integer, tt: Term): RewriteRule;
  {'rr, idx, tt};

% Return index
asserted inline procedure core_getRuleIndex(r: RewriteRule): Integer;
  cadr r;

% Return term
asserted inline procedure core_getRuleTerm(r: RewriteRule): Term;
  caddr r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basistracker

% Implements the Basistracker struct. Basistracker instance is an object
% that stores all polynomials produced during f5 execution in a single vector.
% More precicely, a Basistracker instance is represented as
%   {'bt, polys, filled, capacity}
% Where
%   polys is a Vector of `Polynomial`s,
%   filled is an Integer equal to the number of elements added to `polys`,
%   capacity is an Integer equal to the length of `polys`
%
% In a single core_groebner1 call, a single Basistracker object is created.
% A LabeledPolynomial can be added to the current Basistracker with `core_addPoly`,
% and accessed with `core_getPoly`. To get the index of the last polynomial added
% use `core_getBasisIdx`.
%
% It safe to assume that polynomials added to the current Basistracker
% are normalized (either divided by the leading coeff if f5integers if OFF,
%                 or divided by the content if f5integers is ON).

fluid '(core_initialBasisSize!*);
% We pick a big number from the start,
% so that there is no need to copy the storage vector often
% TODO: scale this with the basis growing
core_initialBasisSize!* := 10000;

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
asserted procedure core_pairTotalDegreeCmp(p1: CriticalPair,
                                            p2: CriticalPair): Boolean;
  poly_totalDegTerm(core_getPairLcm(p1)) #< poly_totalDegTerm(core_getPairLcm(p2));

% compare critical pairs by lcm term according to the current term order
asserted procedure core_pairLcmCmp(p1: CriticalPair,
                                    p2: CriticalPair): Boolean;
  poly_cmpTerm(core_getPairLcm(p1), core_getPairLcm(p2));

% compare associative list elements by their signature
asserted procedure core_assocSgnCmp(pr1: DottedPair,
                                      pr2: DottedPair): Boolean;
  lp_cmpSgn(lp_sgn(cdr pr1), lp_sgn(cdr pr2));

% compare associative list elements by their leading term
% in the current term order
asserted procedure core_assocLeadCmp(pr1: DottedPair,
                                      pr2: DottedPair): Boolean;
  poly_leadTotalDegreeCmp(lp_eval(cdr pr1), lp_eval(cdr pr2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the given list of polynomials `inputBasis` first
%   . interreduce input polynomials
%   . sort them by the order on their leading terms
%     from the smallest to the largest
% and then construct the standard basis of unit vectors corresponding to inputBasis, i.e.,
%   {f, g}  --> {(1, 0), (0, 1)}
% where each unit vector is represented as a LabeledPolynomial
asserted procedure core_constructModule(inputBasis: List): List;
  begin scalar outputModule;
        integer i;
    % If f5integers is on, then scale input polynomials, so that
    % coefficients of input polynomials become integers
    if !*f5integers then
      inputBasis := for each poly in inputBasis
                      collect poly_scaleDenominators(poly);
    % Interreducing input basis is a heuristic. The idea it to produce
    % polynomials with disjoint leading terms after interreduction if possible.
    inputBasis := core_interreduceInput(inputBasis);
    % Since F5 constructs the basis for the input polynomials {f1...fn} incrementally,
    % starting from constructing the basis for {f1,f2}, then for {f1,f2,f3}, and so on,
    % the order of input actually matters.
    % We sort polynomials w.r.t. the total degree of leading terms
    % from the smallest to the largest
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
% During reductions, the topReduce flag controls the reduction type.
% If topReduce is set, only top-reductions happen.
asserted procedure core_normalForm(f: Polynomial, Gprev: List,
                                    r: Basistracker, topReduce: Boolean): DottedPair;
  begin scalar updated, reducer, reduced, updatedToreturn;
    updated := t;
    % while polynomial f gets updated by reductions
    % scan the list Gprev in search for possible reducers
    while updated do <<
      updated := nil;
      for each g in Gprev do <<
        reducer := lp_eval(core_getPoly(r, g));
        if not poly_iszero!?(reducer) and not poly_iszero!?(f) then <<
          reduced . f := if topReduce then
            poly_tryTopReductionStep(f, reducer)
          else
            poly_tryReductionStep(f, reducer);
          updated := reduced or updated;
          updatedToreturn := reduced or updatedToreturn
        >>
      >>
    >>;
    return updatedToreturn . f
  end;

% Same as the above, but all possible reducers are
% already stored in `reducers` as Polynomial objects.
asserted procedure core_normalFormReducers(f: Polynomial, reducers: List,
                                              topReduce: Boolean): DottedPair;
  begin scalar updated, reducer, reduced, updatedToreturn;
    updated := t;
    while updated do <<
      updated := nil;
      for each reducer in reducers do <<
        if not poly_iszero!?(reducer) and not poly_iszero!?(f) then <<
          reduced . f := if topReduce then
            poly_tryTopReductionStep(f, reducer)
          else
            poly_tryReductionStep(f, reducer);
          updated := updated or reduced;
          updatedToreturn := updated or reduced
        >>
      >>
    >>;
    return updatedToreturn . f
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% returns G without the i-th element (starting from index 1), e.g.,
% for G = {8, 4, 5} and i = 2 returns {8, 5}
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

% Given the set of basis generators `input`,
% applies several (possibly, not all) passes of the Autoreduction algorithm
% until no further reductions (possibly, top-reductions) happen
asserted procedure core_interreduceInput(input: List): List;
  begin scalar reducers, updated, reduced, f, newInput;
    updated := t;
    while updated do <<
      updated := nil;
      newInput := nil;
      while input do <<
        f := pop(input);
        reducers := append(newInput, input);
        % Computing only top reductions here. The idea is that the
        % input polynomials are (usually) relatively simple,
        % and full interreduction will not be of any help;
        % Still, top-reductions can help us detect input polynomials
        % with disjoint leading terms, and discard redundant polynomials, if any
        reduced . f := core_normalFormReducers(f, reducers, t);
        updated := updated or reduced;
        if not poly_iszero!?(f) then
          push(f, newInput)
      >>;
      input := newInput
    >>;
    return newInput
  end;

asserted procedure core_interreduceBasis(Gprev: List, r: Basistracker): List;
  begin scalar reducers, updated, reduced, f, newGprev;
        integer i;
    updated := t;
    while updated do <<
      updated := nil;
      newGprev := nil;
      while Gprev do <<
        i := pop(Gprev);
        f := lp_eval(core_getPoly(r, i));
        reducers := append(newGprev, Gprev);
        reduced . f := core_normalForm(f, reducers, r, nil);
        updated := updated or reduced;
        lp_setEval(core_getPoly(r, i), f);
        if not poly_iszero!?(f) then
          push(i, newGprev)
      >>;
      Gprev := newGprev
    >>;
    return newGprev
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part of the main F5 criterion:
%   checks if Term m is top-reducible by the previously computed Groebner basis,
%   which is r_i for i ∈ Gprev;
%   If it is, returns t;
%   Otherwise, returns nil
asserted procedure core_isTopReducibleTerm(m: Term, Gprev: List,
                                            r: Basistracker): Boolean;
  begin scalar glead, isReducible;
        integer gi;
    while Gprev and (not isReducible) do <<
      gi := pop(Gprev);
      glead := poly_leadTerm(lp_eval(core_getPoly(r, gi)));
      if poly_dividesTerm!?(glead, m) then
        isReducible := t
    >>;
    return isReducible
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adds rewrite rule with index k and term term(sgn) to the set
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
  begin scalar sgn, usgnt, foundRewriter, rulesAtSgn, ruleCurrent;
        integer rewriter;
    rewriter := k;
    sgn      := lp_sgn(core_getPoly(r, k));
    usgnt    := lp_termSgn(lp_mulSgn(sgn, u));
    rulesAtSgn := getv(Rule, lp_indexSgn(sgn));
    % Being careful with order of rules in rulesAtSgn:
    % we want to traverse the rules for the signature `sgn`
    % from the last to the first added ones,
    % and terminate as soon as the first matching rule is found
    while rulesAtSgn and (not foundRewriter) do <<
      ruleCurrent := pop(rulesAtSgn);
      foundRewriter := poly_dividesTerm!?(core_getRuleTerm(ruleCurrent), usgnt);
      if foundRewriter then
        rewriter := core_getRuleIndex(ruleCurrent)
    >>;
    return rewriter
  end;

% The Main Rewrite Criterion: Is signature u*sgn(r_k) rewritable by the Rule?
asserted inline procedure core_isRewritable(u: Term, k: Integer,
                                        r: Basistracker, Rule: Vector);
  not (core_findRewriting(u, k, r, Rule) #= k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filters polynomials indexed by Gprev, leaving only those with the leading term
% not divisible by the leading term of any other polynomial
asserted procedure core_filterRedundant(Gprev: List, r: Basistracker): List;
  begin scalar Gnew, alG, Gsorted, glead;
        integer gi;
    % Sort indices Gprev w.r.t. leading terms of polynomials from in r
    alG := for each i in Gprev collect i . core_getPoly(r, i);
    alG := sort(alG, 'core_assocLeadCmp);
    Gsorted := for each pr in alG collect car pr;
    while Gsorted do <<
      gi := pop(Gsorted);
      glead := poly_leadTerm(lp_eval(core_getPoly(r, gi)));
      if not core_isTopReducibleTerm(glead, Gnew, r) then
        push(gi, Gnew)
    >>;
    return Gnew
  end;

% Normalizes each generator in the given `basis`.
% If f5integers is OFF, then each polynomial
% in the basis is divided by its leading coeff;
% Otherwise, this divides each polynomial in
% the basis by the content.
asserted procedure core_normalizeBasis(basis: List): List;
  for each x in basis collect lp_normalize(x);

% given a Groebner `basis` standardizes it,
% so that the output agrees with the following invariants:
%  . the basis is interreduced (if f5fullreduce is OFF, only top-reductions are performed)
%  . the basis is normalized
%  . the basis is sorted by the order on the leading terms of its generators
asserted procedure core_standardizeOutput(basis: List): List;
  begin scalar normalizedBasis;
    normalizedBasis := core_normalizeBasis(basis);
    % Our coefficients are integers if f5integers is ON, so
    % we transform each coefficient back to a Standard Quotient
    if !*f5integers then
      for each poly in normalizedBasis do
        lp_setEval(poly, poly_int2sqCoeffs(lp_eval(poly)));
    return sort(normalizedBasis, 'lp_cmpLPLeadReverse)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if the ideal generated by <polys> contains
% in the ideal generated by <basis>, assuming `basis` is a Groebner basis
asserted procedure core_checkIdealInclusion1(basis: List, polys: List);
  begin scalar included, tmp, poly, nf, flag_;
    % Computes the normal form of each polynomial from `polys`
    % with respect to the Groebner basis from `basis`
    included := t;
    basis    := for each f in basis collect lp_eval(f);
    while included and polys do <<
      poly := lp_eval(pop(polys));
      flag_ . nf := core_normalFormReducers(poly, basis, t);
      if not poly_iszero!?(nf) then
        included := nil
    >>;
    return included
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if `basis` is a Groebner basis by checking that
% all S-polynomials reduce to zero
asserted procedure core_isGroebner1(basis: List);
  begin scalar isBasis, tmp, spoly, nf, flag_;
    isBasis := t;
    basis   := for each f in basis collect lp_eval(f);
    while basis do <<
      tmp := cdr basis;
      while isBasis and tmp do <<
        spoly := poly_spoly(car tmp, car basis);
        flag_ . nf := core_normalFormReducers(spoly, basis, t);
        if not poly_iszero!?(nf) then
          isBasis := nil;
        pop(tmp)
      >>;
      pop(basis)
    >>;
    return isBasis
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For polynomial r_k at index k, it searches for a possible top-reducer from
% the Groebner basis computed up to the previous signature index.
% A possible top-reducer should fullfil several criteria:
% First, it should be a top-reducer of r_k in a usual sense:
%   the leading term of r_k must be divisible by the leading term of top-reducer
% Secondly, top-reducer should admit a proper signature:
%   . The signature must be not equal to the signature of r_k,
%   . Top-reducer should be not rewritable AND not top-reducible itself
%
% If such top-reducer is found, return its index into r;
% if no such top-reducer was found, return 0
asserted procedure core_findReducer(k: Integer, Gprev: List, newGcurr: List,
                                  r: Basistracker, Rule: Vector): Integer;
  begin scalar tt, rk, rkEval, rkSgn, rj, u, rjSgnMult, tj;
        integer j, reducer;
    rk     := core_getPoly(r, k);
    rkEval := lp_eval(rk);
    rkSgn  := lp_sgn(rk);
    tt := poly_leadTerm(rkEval);
    for each j in newGcurr do <<
      rj := core_getPoly(r, j);
      tj := poly_leadTerm(lp_eval(rj));
      % TODO: check the order of args
      if (reducer #= 0) and poly_dividesTerm!?(tj, tt) then <<
        u := poly_divTerm(tt, tj);
        rjSgnMult := lp_mulSgn(lp_sgn(rj), u);
        if not lp_eqSgn(rjSgnMult, rkSgn) then
          if not core_isRewritable(u, j, r, Rule) then
            if not core_isTopReducibleTerm(lp_termSgn(rjSgnMult), Gprev, r) then
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
% is no data left in that particular Labeled polynomial - an
% empty ordered pair is returned.
%
% Otherwise top-reduction calls upon another sub-routine core_findReducer.
% Essentially, if core_findReducer returns 0, the current LabeledPolynomial
% is normalized and then it is added to the list `completed`
% defined within the core_reduction function.
%
% If a top-reduction is deemed possible, then there are two cases:
% either the reduction will increase the signature of the polynomial or it will not.
% In the latter case, the signature is preserved, the polynomial part
% is top-reduced and the LabeledPolynomial is added to the list `todo`
% in the core_reduction function.
% In the former case, however, the signature will change. This
% is marked by adding a new polynomial to r with appropriate
% signature which is equal to the signature of the top-reducer.
% A new rewrite rule is added and then
% both top-reduced and top-reducer polynomials and are added to the list `todo`
% in the core_reduction function. This is done because the
% top-reduced has a different signature than top-reducer and top-reducer might
% still be reducible by another LabeledPolynomial.
asserted procedure core_topReductionF5(k: Integer, Gprev: List, newGcurr: List,
                      r: Basistracker, Rule: Vector): DottedPair;
  begin scalar rk, rj, u, newpoly, reduced, reducerSgn, flag_;
        integer j;
    rk := core_getPoly(r, k);
    % if reduction to zero happened in the normal form -
    % meaning the system is not regular
    if lp_iszero!?(rk) then <<
      %  Reduction to zero! happened;
      return nil . nil
    >>;
    j := core_findReducer(k, Gprev, newGcurr, r, Rule);
    % no top reducers found in the previous Groebner basis --
    % top-reduction is not possible
    if j #= 0 then <<
      core_setPoly(r, k, lp_normalize(rk));
      return k . nil
    >>;
    % reducer rj found, perform top-reduction
    rj := core_getPoly(r, j);
    flag_ . reduced := poly_tryTopReductionStep(lp_eval(rk), lp_eval(rj));
    if not poly_iszero!?(reduced) then
      reduced := poly_normalize(reduced);
    % we need to check that the signature of reducer
    % is not greater than the signature of rj
    u := poly_divTerm(poly_leadTerm(lp_eval(rk)), poly_leadTerm(lp_eval(rj)));
    reducerSgn := lp_mulSgn(lp_sgn(rj), u);
    return if lp_cmpSgn(reducerSgn, lp_sgn(rk)) then <<
      % signatures OK, reduction is successful.
      % We form the new polynomial at index k,
      % and return it to reduction, to be reduced once again later
      lp_setEval(rk, reduced);
      core_setPoly(r, k, rk);
      nil . k . nil
    >> else <<
      % the signature of reducer is greater, reduction cannot be performed.
      % BUT, we still can form the new polynomial with the signature of reducer
      % for further reduction
      newpoly := lp_LabeledPolynomial2(reduced, reducerSgn);
      core_addPoly(r, newpoly);
      core_addRule(Rule, lp_sgn(newpoly), core_getBasisIdx(r));
      nil . k . core_getBasisIdx(r) . nil
    >>
  end;

% inserts integer j in todo and returns the new resulting list,
% so that indices in todo are sorted increasingly by signatures
% of labeled polynomials indexed in r by todo
asserted procedure core_insertSorted(todo: List, j: Integer,
                                        r: Basistracker): List;
  begin scalar tmp, sj, s;
    if null todo then
      return {j};
    sj := lp_sgn(core_getPoly(r, j));
    s  := lp_sgn(core_getPoly(r, car todo));
    if lp_cmpSgn(sj, s) then
      return j . todo;
    tmp := todo;
    while cdr todo and lp_cmpSgn(sj, lp_sgn(core_getPoly(r, cadr todo))) do <<
      pop(todo)
    >>;
    cdr todo := j . cdr todo;
    return tmp
  end;

% The main reduction function.
% Given indices of S-polynomials in the list S, computes
% the F5-style reduction of each of the polynomial indexed by S.
%
% By F5-style reduction we mean a normal form,
% with respect to the previously computed Groebner basis,
% combined with a top-reduction that does not change the signature.
asserted procedure core_reduction(S: List, Gprev: List, Gcurr: List,
                                    r: Basistracker, Rule: Vector): List;
  begin scalar todo, S, completed, newGcurr, rk, reduced,
                newcompleted, redo, flag_;
        integer k, j;
    % For each polynomial f the reduction process is the following.
    % First, this polynomial is placed in the `todo` list from the beginning.
    % Secondly, when f is popped from the list, the normal form of f with
    % respect to the basis computed up to this module index (in Gprev) is produced.
    % Since during the normal form computation it is possible that not all feasible
    % reducers were handled, a separate top-reduction step is needed.
    % This is handled by the core_topReductionF5, which is called on the
    % result of the normal form.
    %
    % The core_topReductionF5 returns a pair of two lists:
    %   newcompleted and redo
    % Elements from redo are inserted back to the todo list,
    % and the elements of newcompleted are considered fully reduced.
    todo := S;
    newGcurr := copy(Gcurr);
    while todo do <<
      k := pop(todo);
      rk := core_getPoly(r, k);
      % Compute the normal form;
      % If full reduction is not needed, computes only the top normal form.
      % Otherwise, reduces all terms in the polynomial.
      flag_ . reduced := core_normalForm(lp_eval(rk), Gprev, r, t);
      lp_setEval(rk, reduced);
      core_setPoly(r, k, rk);
      % Compute the F5-style top-reduction;
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
%  S-polynomials output from core_makeCriticalPair as LabeledPolynomials.
%  We note that, because core_makeCriticalPair ensured
%  that sgn(u*k) > sgn(v*l), we know that the signature of all new
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
        ssgn  := lp_mulSgn(lp_sgn(lpk), u);
        core_addPoly(r, lp_LabeledPolynomial2(seval, ssgn));
        core_addRule(Rule, ssgn, core_getBasisIdx(r));
        if not poly_iszero!?(seval) then
          push(core_getBasisIdx(r), S)
      >>
    >>;
    % Sort indices S w.r.t. signatures of polynomials stored in r
    alS := for each i in S collect i . core_getPoly(r, i);
    alS := sort(alS, 'core_assocSgnCmp);
    S := for each pr in alS collect car pr;
    return S
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constructs a CriticalPair from polynomials at indices k and l.
% Returns a CriticalPair object.
% If the pair is deemed redundant by the F5 or the rewritten criterion, returns nil.
asserted procedure core_makeCriticalPair(i: Integer, k: Integer,
                        l: Integer, Gprev: List,
                        r: Basistracker, Rule: Vector): CriticalPair;
  begin scalar rk, rl, tk, tl, tt, u1, u2, usgn1, usgn2;
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
    if (lp_indexSgn(usgn1) #= i) and core_isTopReducibleTerm(lp_termSgn(usgn1), Gprev, r) then
      return nil;
    if (lp_indexSgn(usgn2) #= i) and core_isTopReducibleTerm(lp_termSgn(usgn2), Gprev, r) then
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
% constructs the basis of {f1..fi} and returns a list of new basis indices
asserted procedure core_incrementalBasis(i: Integer, Gprev: List,
                          r: Basistracker, Rule: Vector): List;
  begin scalar Gcurr, pairs, p, S, dpairs, reduced, k, tmp;
        integer i, j, d, currIdx;
    % The function is organized in the following way.
    % In the very beginning, some initial critical pairs are constructed
    % and added to the list `pairs`. Pairs are formed using
    % the r_currIdx polynomial and the generators of {f1..fi-1}.
    % Then, while the `pairs` list is not empty, all pairs of the smallest
    % degree of the lcm part are chosen and deleted from `pairs`. These pairs
    % are then transformed into S-polynomials with the `core_computeSpolys`.
    % S-polynomials are passed into the `core_reduction`, which returns
    % the indices of new basis members. These new polynomials
    % are used to produce new critical pairs, which are inserted in `pairs`.
    currIdx := core_getBasisIdx(r);
    Gcurr   := copy(Gprev);
    push(currIdx, Gcurr);
    for each j in Gprev do <<
      p := core_makeCriticalPair(i, currIdx, j, Gprev, r, Rule);
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
          p := core_makeCriticalPair(i, j, k, Gprev, r, Rule);
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
% the function takes a list of LabeledPolynomials as ideal generators,
% and returns a standardized Groebner basis of this ideal as output
%
% Output contracts:
%   . the output basis is top-reduced, meaning the number of generators is minimal,
%   . the output basis is sorted,
%   . the output basis contains only normalized polynomials,
%   . if f5fullreduce is on, then the output basis is fully reduced
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
    % Incremental construction from index 1 to index m-1
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
    if !*f5fullreduce then
      Gprev := core_interreduceBasis(Gprev, r);
    basis := for each i in Gprev collect core_getPoly(r, i);
    return core_standardizeOutput(basis)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trst poly_tryTopReductionStep;
% trst core_computeSpolys;
% trst core_topReductionF5;
% trst core_groebner1;
% trst core_constructModule;
% trst core_interreduceInput;
% trst core_normalFormReducers;
% trst core_incrementalBasis;
% trst core_makeCriticalPair;
% trst core_isTopReducibleTerm;
% trst core_reduction;
% trst core_findReducer;
% trst core_findRewriting;
% trst core_filterRedundant;
% trst core_isGroebner1;
% trst core_checkIdealInclusion1;

endmodule;

end; % of file
