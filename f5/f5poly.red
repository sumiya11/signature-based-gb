module f5poly;
% Polynomial interface module to be used in f5.
% The module provides procedures for basic operations with the `Polynomial` type.

% Polynomial `p` is stored as a list of 4 items:
%     {'p, Terms, Coeffs, Sugar}
% Where `p` is a convenience tag,
%       `Terms` is a list of `Term`s. Each `Term` is an exponent list of
%         non-negative integers of form
%              {totaldegree, pow1, pow2, ... pown}
%              `Term` supports basic arithmetic operations,
%              the interface is defined further in this file.
%       `Coeffs` is a list of `Coeff`s,
%         where `Coeff` can be either a SQ or an Integer.
%       `Sugar` is an Integer, the sugar degree.
% Some relevant functions on `Coeff` are defined further in this file:
%     . poly_addCoeff(x, y)  -- addition x + y
%     . poly_divCoeff(x, y)  -- division x / y
%     . poly_invCoeff(x)     -- inverse  x^(-1)
%     . poly_negCoeff(x)     -- negation -x
%     . poly_mulCoeff(x, y)  -- product  x * y
%     . poly_iszeroCoeff(x)  -- zero?    x
%
% `Terms` are ordered according to the current term order decreasingly,
%  and `Coeffs` are ordered respectively.
%
% For example, xy^2 + 3x is stored as
%   {'p, {{3, 1, 2}, {1, 1, 0}}, {1, 3}, 3}
% if f5integers is ON. It is stored as
%   {'p, {{3, 1, 2}, {1, 1, 0}}, {1 ./ 1, 3 ./ 1}, 3}
% otherwise (using SQ).
%
% Zero polynomial is represented as
%   {'p, nil, nil, anything}
%
% The global polynomial ring should be initialized before constructing polynomials.
% To initialize the ring in variables `vars` and term order `ord`
% `poly_initRing(vars, ord)` should be used.

off1 'allfac;

% We use the term order comparators for exponent lists from dipoly
load!-package 'dipoly;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns the current torder
asserted procedure poly_extractTorder(): List;
   begin scalar oldTorder;
      oldTorder := torder('(list));
      torder(cdr oldTorder);
      dipsortingmode(vdpsortmode!*);
      return oldTorder
   end;

% Initializes polynomial ring with the parameters from `u`,
% and returns the corresponding torder.
%
% Possible options are:
%  - u is nil, then change nothing;
%  - u is of length 1, then set variables in torder to the only element from u;
%  - u is of length 2, then set variables and order in torder to the ones from u.
asserted procedure poly_initRing(u: List): List;
   begin scalar vars, ord, oldTorder;
      if null u then <<
         oldTorder := poly_extractTorder()
      >> else if length(u) = 1 then <<
         vars := pop u;
         oldTorder := poly_extractTorder();
         global!-dipvars!* := 'list . vars
      >> else <<
         vars := pop u;
         ord  := pop u;
         oldTorder := torder({'list . vars, ord})
      >>;
      return oldTorder
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% POLYNOMIAL INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constructor of Polynomial, forms a Polynomial
% from a list of `Term`s, a list of `Coeff`s, and a sugar degree.
asserted inline procedure poly_PolynomialWithSugar(ts: Terms, cfs: Coeffs, sugar: Integer): Polynomial;
   {'p, ts, cfs, sugar};

% Same as above, but doesn't care about sugar
asserted inline procedure poly_Polynomial(ts: Terms, cfs: Coeffs): Polynomial;
   poly_PolynomialWithSugar(ts, cfs, 0);

asserted inline procedure poly_getTerms(poly: Polynomial): Terms;
   cadr poly;

asserted inline procedure poly_getCoeffs(poly: Polynomial): Coeffs;
   caddr poly;

asserted inline procedure poly_getSugar(poly: Polynomial): Integer;
   cadddr poly;

% Returns zero polynomial
asserted inline procedure poly_zero(): Polynomial;
   poly_Polynomial(nil, nil);

% Checks if `p` is zero
asserted inline procedure poly_iszero!?(p: Polynomial): Boolean;
   null poly_getTerms(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% POLYNOMIAL CONVERSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constructs a Polynomial from a SF `u`
asserted inline procedure poly_f2poly(u: SF): Polynomial;
   poly_f2poly1(u, poly_zeroExp(), poly_oneCoeff());

% Conversion to Polynomial: scan the standard form. ev and bc are the
% exponent and coefficient parts collected so far from higher parts.
asserted procedure poly_f2poly1(u: SF, ev: List, bc: Coeff): Polynomial;
  if null u then
     poly_zero()
  else if domainp u then
     poly_PolynomialWithSugar({ev}, {poly_mulCoeff(bc, poly_2Coeff(u))}, poly_totalDegExp(ev))
  else
     poly_sumPoly(poly_f2poly2(mvar u,ldeg u,lc u,ev,bc), poly_f2poly1(red u,ev,bc));

% Conversion to Polynomial: multiply leading power either into exponent vector
asserted procedure poly_f2poly2(var,dg,c,ev,bc): Polynomial;
   poly_f2poly1(c,poly_insertExp(ev, var, dg, cdr global!-dipvars!*), bc);

% Returns prefix equivalent to the sum of elements of u
asserted procedure poly_replus(u: List): List;
   if atom u then u else if null cdr u then car u else 'plus . u;

% Returns prefix equivalent to the product of elements of u.
% u is a list of prefix expressions the first of which is a number.
procedure poly_retimes(u: List): List;
   if car u = 1 then
      if cdr u then poly_retimes cdr u else 1
   else if null cdr u then
      car u
   else
      'times . u;

% Returns prefix equivalent to the Polynomial u.
asserted procedure poly_2a1(u: Polynomial): List;
   begin scalar x,y;
      if poly_iszero!?(u) then
         return nil;
      x := poly_leadCoeff u;
      y := poly_2aExp poly_leadTerm u;
      if poly_isNegCoeff!?(x) then <<
         return {'minus,poly_retimes(poly_2aCoeff(poly_negCoeff(x)) . y)} . poly_2a1 poly_tail(u)
      >>;
      return poly_retimes(poly_2aCoeff x . y) . poly_2a1 poly_tail(u)
   end;

% Returns prefix equivalent to the Polynomial f.
asserted procedure poly_2a(f: Polynomial): List;
   if poly_iszero!?(f) then 0 else poly_replus poly_2a1(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EXPONENT LISTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invariant: the first entry in the exponent list is the sum of subsequent entries

% Returns the first entry in the exponent list - the total degree
asserted inline procedure poly_totalDegExp(e1: List): Integer;
   car e1;

% Returns exponent list of zeros of the appropriate length
asserted inline procedure poly_zeroExp(): List;
   for x := 1:length(global!-dipvars!*) collect 0;

asserted procedure poly_sumExp(e1: List, e2: List): List;
   <<
      ASSERT(not null e1 and not null e2);
      if car e1 #= 0 then
         e2
      else if car e2 #= 0 then
         e1
      else
         for each x in e1 collect x #+ pop e2
   >>;

% Return the elementwise subtraction of exponent lists e1, e2
asserted procedure poly_subExp(e1: List, e2: List): List;
   <<
      ASSERT(not null e1 and not null e2);
      for each x in e1 collect x #- pop e2
   >>;

% Returns the elementwise maximum of exponent lists e1, e2
asserted procedure poly_elmaxExp(e1: List, e2: List): List;
   begin scalar w;
      ASSERT(not null e1 and not null e2);
      e1 := cdr e1;
      e2 := cdr e2;
      w := for each x in e1 collect max(x, pop e2);
      return (for each x in w sum x) . w
   end;

% Checks if e1 is elementwise less or equal than e2
asserted procedure poly_divExp!?(e1: List, e2: List): Boolean;
   begin scalar brk;
      % We deliberately also check the total degree.
      while not brk and e1 do
         brk := pop e1 #> pop e2;
      return not brk
   end;

% Checks that at least one of e1[i] or e2[i] is zero for each i
asserted procedure poly_disjExp!?(e1: List, e2: List): Boolean;
   begin scalar ok;
      e1 := cdr e1;
      e2 := cdr e2;
      ok = t;
      while ok and e1 do
         ok := pop e1 #= 0 or pop e2 #= 0;
      return ok
   end;

% Insert the "dg" into the ev in the place of variable v
asserted procedure poly_insertExp(ev: List, v: Any, dg: Integer, vars: List): List;
   (car ev #+ dg) . poly_insertExp1(cdr ev, v, dg, vars);

asserted procedure poly_insertExp1(ev: List, v: Any, dg: Integer, vars: List): List;
   if null ev or null vars then
      nil
   else if car vars eq v then
      dg . cdr ev
   else
      car ev . poly_insertExp1(cdr ev, v, dg, cdr vars);

asserted inline procedure poly_2aExp(e);
   % Returns list of prefix equivalents of exponent vector e.
   ev_2aExp1(cdr e,  cdr global!-dipvars!*);

procedure ev_2aExp1(u,v);
   if null u then
      nil
   else if car u #= 0 then
      ev_2aExp1(cdr u,cdr v)
   else if car u #= 1 then
      car v . ev_2aExp1(cdr u,cdr v)
   else
      {'expt,car v,car u} . ev_2aExp1(cdr u,cdr v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparators for exponent lists.
% We mainly fall back to functions from dipoly package here

% Compares exponent lists e1, e2 w.r.t. lex term order,
% and returns e1 < e2
asserted inline procedure poly_cmpExpLex(e1: List, e2: List): Boolean;
   evlexcomp(cdr e1, cdr e2) #= -1;

% Compares exponent lists e1, e2 w.r.t. graded lex term order,
% and returns e1 < e2
asserted inline procedure poly_cmpExpGradLex(e1: List, e2: List): Boolean;
  evgradlexcomp(e1, e2) #= -1;

% Compares exponent lists e1, e2 w.r.t. graded reversed lex term order,
% and returns e1 < e2
asserted inline procedure poly_cmpExpRevGradLex(e1: List, e2: List): Boolean;
   if car e1 #= car e2 then
      evinvlexcomp(cdr e1, cdr e2) #= -1
   else
      car e1 #< car e2;

% Compares exponent lists e1, e2 w.r.t. the current order in torder
asserted inline procedure poly_cmpExpGeneric(e1: List, e2: List): Boolean;
   evcomp(cdr e1, cdr e2) = -1;

% Compares exponent lists e1, e2 w.r.t. the current term order
asserted procedure poly_cmpExp(e1: List, e2: List): Boolean;
   <<
   if vdpsortmode!* eq 'lex then
      poly_cmpExpLex(e1, e2)
   else if vdpsortmode!* eq 'gradlex then
      poly_cmpExpGradLex(e1, e2)
   else if vdpsortmode!* eq 'revgradlex then
      poly_cmpExpRevGradLex(e1, e2)
   else  % Trouble with vdpmatrix!*
      poly_cmpExpGeneric(e1, e2)
   >>;

% Compares exponent lists w.r.t. the total degree
asserted inline procedure poly_tdegCmpExp(e1: List, e2: List): Boolean;
   car e1 #< car e2;

% Checks that e1 = e2 elementwise
asserted inline procedure poly_eqExp!?(e1: List, e2: List): Boolean;
   e1 = e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Term is a polynomial term, which is implemented as an exponent list.

% Returns the identity Term (just one).
asserted inline procedure poly_identityTerm(): Term;
   poly_zeroExp();

asserted inline procedure poly_isIdentityTerm!?(tm: Term): Boolean;
   poly_totalDegTerm(tm) #= 0;

% Returns the total degree of the Term
asserted inline procedure poly_totalDegTerm(a: Term): Integer;
   poly_totalDegExp(a);

% Returns a * b
asserted inline procedure poly_mulTerm(a: Term, b: Term): Term;
   poly_sumExp(a, b);

% Returns a / b,
% Assuming a >= b.
asserted inline procedure poly_divTerm(a: Term, b: Term): Term;
   poly_subExp(a, b);

% Checks that a | b
asserted inline procedure poly_dividesTerm!?(a: Term, b: Term): Boolean;
   poly_divExp!?(a, b);

% Returns lcm(a, b)
asserted inline procedure poly_lcmTerm(a: Term, b: Term): Term;
   poly_elmaxExp(a, b);

% Returns a < b in the current term order term order
asserted inline procedure poly_cmpTerm(a: Term, b: Term): Boolean;
   poly_cmpExp(a, b);

% Checks if gcd(a, b) is one
asserted inline procedure poly_disjTerm!?(a: Term, b: Term): Boolean;
   poly_disjExp!?(a, b);

% Checks if a = b
asserted inline procedure poly_eqTerm!?(a: Term, b: Term): Boolean;
   poly_eqExp!?(a, b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% POLYNOMIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns the tail of the polynomial `poly`.
% That is, `poly - lead(poly)`,
% and preserves the sugar degree.
asserted inline procedure poly_tail(poly: Polynomial): Polynomial;
   poly_PolynomialWithSugar(poly_tailTerms(poly), poly_tailCoeffs(poly), poly_getSugar(poly));

% Returns the leading term of `poly`
asserted inline procedure poly_leadTerm(poly: Polynomial): Term;
   car poly_getTerms(poly);

% Returns the leading coefficient of `poly`
asserted inline procedure poly_leadCoeff(poly: Polynomial): Coeff;
   car poly_getCoeffs(poly);

% Returns the tail terms of `poly`
asserted inline procedure poly_tailTerms(poly: Polynomial): Terms;
   cdr poly_getTerms(poly);

% Returns the tail coefficients of `poly`
asserted inline procedure poly_tailCoeffs(poly: Polynomial): Coeffs;
   cdr poly_getCoeffs(poly);

% Returns the length of `poly`, i.e., the number of terms
asserted inline procedure poly_length(poly: Polynomial): Integer;
   length(poly_getTerms(poly));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ADAPTIVE COEFFICIENT ARITHMETIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted inline procedure poly_iszeroCoeff!?(a: Coeff): Boolean;
   if !*f5modular then
      a #= 0
   else if !*f5integers then
      a = 0
   else
      not numr(a);

asserted inline procedure poly_isoneCoeff!?(a: Coeff): Boolean;
   if !*f5modular then
      a #= 1
   else if !*f5integers then
      a = 1
   else
      eqn(numr(a), 1) and eqn(denr(a), 1);

asserted inline procedure poly_oneCoeff(): Coeff;
   poly_2Coeff(1);

asserted inline procedure poly_2Coeff(c: Any): Coeff;
   if !*f5modular then
      modular!-number(1)
   else if !*f5integers then
      c
   else
      c ./ 1;

asserted inline procedure poly_2aCoeff(c: Coeff);
   prepsq c;

asserted inline procedure poly_addCoeff(a: Coeff, b: Coeff): Coeff;
   if !*f5modular then
      modular!-plus(a, b)
   else if !*f5integers then
     a + b
   else
      addsq(a, b);

asserted inline procedure poly_subCoeff(a: Coeff, b: Coeff): Coeff;
   if !*f5modular then
      modular!-sub(a, b)
   else if !*f5integers then
      a - b
   else
      subsq(a, b);

asserted inline procedure poly_mulCoeff(a: Coeff, b: Coeff): Coeff;
   if !*f5modular then
      modular!-times(a, b)
   else if !*f5integers then
      a * b
   else
      multsq(a, b);

asserted inline procedure poly_negCoeff(a: Coeff): Coeff;
   if !*f5modular then
      modular!-minus(a)
   else if !*f5integers then
      - a
   else
      negsq(a);

asserted inline procedure poly_isNegCoeff!?(a: Coeff): Boolean;
   if !*f5integers then
      a < 0
   else
      minusf numr a;

asserted inline procedure poly_divCoeff(a: Coeff, b: Coeff): Coeff;
   if !*f5modular then
      modular!-quotient(a, b)
   else if !*f5integers then
      a / b
   else
      quotsq(a, b);

asserted inline procedure poly_invCoeff(a: Coeff): Coeff;
   if !*f5modular then
      modular!-reciprocal(a)
   else if !*f5integers then
      rederr "*****  Trying to take an inverse with f5integers ON."
   else <<
      % ASSERT(not null numr a);
      denr(a) ./ numr(a)
   >>;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% POLYNOMIAL LOW LEVEL OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the only section where polynomial coefficient arithmetic happens

% Returns s = fmult*f - fcoeff*gmult*g
%
% !! Assuming the leading terms of `fmult*f` and `fcoeff*gmult*g`
%    mutually cancel each other
asserted inline procedure poly_paircombTail(f: Polynomial, fmult: Term,
                                                fcoeff: Coeff, g: Polynomial,
                                                gmult: Term): Polynomial;
   % Assuming the leading monomials of `fmult*f` and `fcoeff*gmult*g`
   % are mutually canceled, we can use the general reduction applied to
   % the tails of input polynomials
   poly_paircomb(poly_tail(f), fmult, fcoeff, poly_tail(g), gmult);

% Returns s = fmult*f - fcoeff*gmult*g
asserted procedure poly_paircomb(f: Polynomial,  fmult: Term,
                                 fcoeff: Coeff, g: Polynomial,
                                 gmult: Term): Polynomial;
   begin scalar fterms, fcoeffs, gterms, gcoeffs, gmultcoeff,
                sterms, scoeffs, ft, gt, fc, gc, newc, sugar;
      % We return s, a new polynomial,
      % constructed as s = gcoeff*fmult*f - fcoeff*gmult*g.
      % We form two lists, sterms and scoeffs, which would be the list of
      % terms and the list of coefficients of s.
      % The sterms list is formed by merging two sorted lists:
      % the list of terms of f each multiplied by fmult,
      % with the list of terms of g each multiplied by gmult.
      % In parallel (in the same loop), scoeffs list is formed
      % by merging the list of coefficients of f each multiplied by fmultcoeff,
      % with the list of coefficients of g each multiplied by gmultcoeff,
      % in the same merge order as for the terms.
      %
      fterms  := poly_getTerms(f);
      fcoeffs := poly_getCoeffs(f);
      gterms  := poly_getTerms(g);
      gcoeffs := poly_getCoeffs(g);
      gmultcoeff := poly_negCoeff(fcoeff);

      % Merge two sorted lists: fterms and gterms, multiplied by fmult and gmult, respectively.
      % Merge in the same order two other lists:
      % fcoeffs and gcoeffs, multiplied by gmultcoeff and fmultcoeff, respectively.
      %
      % if n = length(fterms), m = length(gterms), then in the worst case
      %   2(n+m)*E  +    m*E     +     2m*C    +   (n+m)
      % comparison    term mult.    coeff op.     reverse
      %
      % where E is the cost of iterating exponent vector,
      % and C is the cost of one arithmetic operation on coefficients
      while fterms and gterms do <<
         if null ft then <<
            ft := car fterms;
            % identity check is 1 car and 1 comparison
            if not poly_isIdentityTerm!?(fmult) then
               ft := poly_mulTerm(ft, fmult)
         >>;
         if null gt then
            gt := poly_mulTerm(car gterms, gmult);
         % Optimization: return -1,0,1 just as C comparator;
         if poly_cmpTerm(gt, ft) then <<   % if term gt < term ft
            push(ft, sterms);
            push(car fcoeffs, scoeffs);
            pop(fterms); pop(fcoeffs);
            ft := nil
         >> else if poly_eqTerm!?(gt, ft) then <<  % if term gt = term ft
            fc := car fcoeffs;
            gc := poly_mulCoeff(car gcoeffs, gmultcoeff);
            newc := poly_addCoeff(fc, gc);
            if not poly_iszeroCoeff!?(newc) then <<
               push(gt, sterms);
               push(newc, scoeffs)
            >>;
            pop(fterms); pop(fcoeffs);
            pop(gterms); pop(gcoeffs);
            gt := nil;
            ft := nil
         >> else <<   % if term gt > term ft
            push(gt, sterms);
            push(poly_mulCoeff(car gcoeffs, gmultcoeff), scoeffs);
            pop(gterms); pop(gcoeffs);
            gt := nil
         >>
      >>;
      if null gterms and null fterms then <<
         scoeffs := reversip(scoeffs);
         sterms  := reversip(sterms)
      >>;
      % Merge what is left from gterms and gcoeffs
      if not null gterms then <<
         if poly_isIdentityTerm!?(gmult) then
            sterms := nconc(reversip sterms, gterms)
         else <<
            while gterms do
               push(poly_mulTerm(pop(gterms), gmult), sterms);
            sterms := reversip(sterms)
         >>;
         while gcoeffs do
            push(poly_mulCoeff(pop(gcoeffs), gmultcoeff), scoeffs);
         scoeffs := reversip(scoeffs);
      >>;
      % Merge what is left from fterms and fcoeffs
      if not null fterms then <<
         if poly_isIdentityTerm!?(fmult) then
            sterms := nconc(reversip sterms, fterms)
         else <<
            while fterms do
               push(poly_mulTerm(pop(fterms), fmult), sterms);
            sterms := reversip(sterms)
         >>;
         scoeffs := nconc(reversip scoeffs, fcoeffs)
      >>;
      sugar := max(poly_getSugar(f) + poly_totalDegTerm(fmult),
                        poly_getSugar(g) + poly_totalDegTerm(gmult));
      return poly_PolynomialWithSugar(sterms, scoeffs, sugar)
   end;

% for history
procedure copyList(l);
  begin scalar queue, newPair;
    if null l then
      return nil;
    queue := car l . nil;
    queue := queue . queue;
  % car queue points to the end of the queue
  % cdr queue points to the beginning of the queue
    l := cdr l;
    while not null l do <<
      newPair := car l . nil;
      cdar queue := newPair;
      car queue := newPair;
      l := cdr l
    >>;
    return cdr queue
  end;

% Constructs a new polynomial as a copy of `poly` with normalized coefficients.
% If !*f5integers is ON, this will divide all coefficients by the content of `poly`
% (resulting in a polynomial with unit content).
% Otherwise, this will divide all coefficient by the leading coefficient of `poly`.
% (resulting in a monic polynomial).
asserted procedure poly_normalize(poly: Polynomial): Polynomial;
   if !*f5integers then
      poly_normalizeByContent(poly)
   else
      poly_normalizeByLead(poly);

% Constructs a new polynomial as a copy of `poly`
% with coefficients divided by the content of `poly`
asserted procedure poly_normalizeByContent(poly: Polynomial): Polynomial;
   begin scalar newcoeffs, cnt;
      cnt := poly_content(poly);
      % The leading coefficient will be positive after division
      if poly_leadCoeff(poly) < 0 then
         cnt := -cnt;
      newcoeffs := for each cf in poly_getCoeffs(poly)
         collect poly_divCoeff(cf, cnt);
      return poly_PolynomialWithSugar(poly_getTerms(poly), newcoeffs, poly_getSugar(poly))
   end;

% Constructs a new polynomial as a copy of `poly`
% with coefficients divided by the leading coeff of `poly`
asserted procedure poly_normalizeByLead(poly: Polynomial): Polynomial;
   begin scalar newcoeffs, mult1;
      mult1 := poly_invCoeff(poly_leadCoeff(poly));
      newcoeffs := for each cf in poly_getCoeffs(poly)
         collect poly_mulCoeff(cf, mult1);
      return poly_PolynomialWithSugar(poly_getTerms(poly), newcoeffs, poly_getSugar(poly))
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High level functions: Spolynomial computation and reduction

asserted inline procedure poly_sumPoly(f: Polynomial, g: Polynomial): Polynomial;
   poly_paircomb(f, poly_identityTerm(), poly_negCoeff(poly_oneCoeff()), g, poly_identityTerm());

% Computes S-polynomial of f and g
asserted procedure poly_spoly(f: Polynomial, g: Polynomial): Polynomial;
   begin scalar e1, e2, elcm, mult1, mult2;
      e1 := poly_leadTerm(f);
      e2 := poly_leadTerm(g);
      elcm := poly_lcmTerm(e1, e2);
      mult1 := poly_divTerm(elcm, e2);
      mult2 := poly_divTerm(elcm, e1);
      % using poly_paircombTail since leading terms vanish
      return poly_paircombTail(f, mult2, poly_leadCoeff(f), g, mult1)
   end;

% Tries to top-reduce f using g (reduces only the leading term of f, only once).
% Two cases are possible:
% 1. the leading term of g divides the leading term of f. Then we perform
%    the top-reduction f - u*g for some u and return the t flag
%    and the result of reduction.
% 2. top-reduction is not possible, we return the nil flag and f
asserted procedure poly_tryTopReductionStep(f: Polynomial,
                                              g: Polynomial): DottedPair;
   begin scalar glead, flead, fmult, gmult, updated;
      glead := poly_leadTerm(g);
      flead := poly_leadTerm(f);
      if poly_dividesTerm!?(glead, flead) then <<
         fmult := poly_identityTerm();
         gmult := poly_divTerm(flead, glead);
         % using poly_paircombTail since leading terms vanish.
         % Optimization: do not pass fmult;
         % Optimization: poly_leadCoeff(g) is always one;
         f := poly_paircombTail(f, fmult, poly_leadCoeff(f), g, gmult);
         updated := t
      >>;
      return updated . f
   end;

% Tries to reduce f with g (only once). Two cases are possible:
% 1. the leading term of g divides some term of f. Then we perform
%    the reduction of that term, and return t flag and the result of reduction.
% 2. reduction is not possible, we return nil flag and f
asserted procedure poly_tryReductionStep(f: Polynomial,
                                          g: Polynomial): DottedPair;
   begin scalar updated, fterms, fcoeffs, glead,
                fcoef, fterm, fmult, gmult;
      fterms  := poly_getTerms(f);
      fcoeffs := poly_getCoeffs(f);
      glead := poly_leadTerm(g);
      while (not updated) and fterms do <<
         fterm := pop(fterms);
         fcoef := pop(fcoeffs);
         if poly_dividesTerm!?(glead, fterm) then <<
            fmult := poly_identityTerm();
            gmult := poly_divTerm(fterm, glead);
            f := poly_paircomb(f, fmult, fcoef, g, gmult);
            updated := t
         >>
      >>;
      return updated . f
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional polynomial coefficient manipulations

% Assuming coefficients of `poly` are Integers,
% transforms them into Standard Quotients and returns the resulting polynomial
asserted procedure poly_int2sqCoeffs(poly: Polynomial): Polynomial;
   poly_PolynomialWithSugar(
      poly_getTerms(poly),
      for each c in poly_getCoeffs(poly) collect c . 1,
      poly_getSugar(poly)
   );

% Returns the lcm of denominators of coefficients of `f`
asserted procedure poly_commonDenominator(f: Polynomial): Integer;
   begin scalar fcoeffs, den;
      den := 1;
      fcoeffs := poly_getCoeffs(f);
      while fcoeffs do <<
         den := lcm(den, denr(pop(fcoeffs)))
      >>;
      return den
   end;

% Returns the content of `f`
asserted procedure poly_content(f: Polynomial): Integer;
  begin scalar fcoeffs, cnt;
    cnt := poly_leadCoeff(f);
    fcoeffs := poly_tailCoeffs(f);
    while fcoeffs do <<
      cnt := mod_euclid(cnt, pop(fcoeffs))
    >>;
    return abs(cnt)
  end;

% Constructs a new polynomial in the following way:
%   f * inv(poly_commonDenominator(f))
asserted procedure poly_scaleDenominators(f: Polynomial): Polynomial;
  begin scalar fcoeffs, newcoeffs, c, den;
    den := poly_commonDenominator(f);
    fcoeffs := poly_getCoeffs(f);
    while fcoeffs do <<
      c := pop(fcoeffs);
      push(numr(c) * (den / denr(c)), newcoeffs)
    >>;
    return poly_PolynomialWithSugar(poly_getTerms(f), reversip(newcoeffs), poly_getSugar(f))
  end;

% Reduces coefficients of `f` by the given `prime`
asserted procedure poly_reduceCoeffs(f: Polynomial, prime: Integer): Polynomial;
  begin scalar fcoeffs, newcoeffs, c;
    % note that prime is not used here, since `modular!-number` works globally,
    % to elide warning
    prime := prime;
    fcoeffs := poly_getCoeffs(f);
    while fcoeffs do <<
         c  := pop(fcoeffs);
         % ASSERT(denr(c) = 1);
         c := modular!-number(c);
         push(c, newcoeffs)
      >>;
      return poly_PolynomialWithSugar(poly_getTerms(f), reversip(newcoeffs), poly_getSugar(f))
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional modular polynomial coefficient manipulation

% Reconstructs each coefficient of `poly` modulo the given prime
asserted procedure poly_reconstructCoeffs(poly: Polynomial,
                                          prime: Integer): Polynomial;
  begin scalar newcoeffs;
    newcoeffs := for each cf in poly_getCoeffs(poly)
      collect mod_reconstruction(cf, prime);
      return poly_PolynomialWithSugar(poly_getTerms(poly), newcoeffs, poly_getSugar(poly))
  end;

% Applis CRT to coefficients of (polyaccum mod modulo) and (polycomp mod prime)
% to obtain new polynomial over modulo*prime
asserted procedure poly_crtCoeffs(polyaccum: Polynomial, modulo: Integer,
                          polycomp: Polynomial, prime: Integer): Polynomial;
  begin scalar coeffsaccum, coeffscomp, newcoeffs, c;
    coeffsaccum := poly_getCoeffs(polyaccum);
    coeffscomp  := poly_getCoeffs(polycomp);
    while coeffsaccum do <<
      c := mod_crt(pop(coeffsaccum), modulo, pop(coeffscomp), prime);
      push(c, newcoeffs)
    >>;
    return poly_PolynomialWithSugar(poly_getTerms(polyaccum), reversip(newcoeffs), poly_getSugar(polyaccum))
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial sorting ad-hoc

% Returns leadTerm(poly1) < leadTerm(poly2)
asserted inline procedure poly_cmpPolyLead(poly1: Polynomial,
                                            poly2: Polynomial): Boolean;
   poly_cmpTerm(poly_leadTerm(poly1), poly_leadTerm(poly2));

% Returns t if
% the total degree of lead term of poly1 < the total degree of lead term poly2,
% Otherwise, if degrees are equal compares terms with the current order,
% Otherwise, nil.
asserted procedure poly_leadTotalDegreeCmp(poly1: Polynomial,
                                            poly2: Polynomial): Boolean;
   begin integer t1, t2;
      t1 := poly_totalDegTerm(poly_leadTerm(poly1));
      t2 := poly_totalDegTerm(poly_leadTerm(poly2));
      return if t1 #= t2 then
         poly_cmpPolyLead(poly1, poly2)
      else
         t1 #< t2
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endmodule;  % end of module f5poly

end;  % of file
