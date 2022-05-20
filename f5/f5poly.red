module f5poly;
% Polynomial interface module to be used in f5.
% The module provides procedures for basic operations with the `Polynomial` type.

% Polynomial `p` is stored as a list of 3 items:
%     {'p, Terms, Coeffs}
% Where `p` is a convenience tag,
%       `Terms` is a list of `Term`s. Each `Term` is an exponent list of
%         non-negative integers of form
%              {totaldegree, pow1, pow2, ... pown}
%              `Term` supports basic arithmetic operations,
%              the interface is defined further in this file.
%       `Coeffs` is a list of `Coeff`s,
%         where `Coeff` can be either a SQ or an Integer.
%   Some relevant functions on `Coeff` are defined further in this file:
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
%   {'p, {{3, 1, 2}, {1, 1, 0}}, {1, 3}}
% if f5integers is ON. It is stored as
%   {'p, {{3, 1, 2}, {1, 1, 0}}, {1 ./ 1, 3 ./ 1}}
% otherwise (using SQ).

% The global polynomial ring should be initialized before constructing polynomials.
% To initialize the ring in variables `vars` and term order `ord`
% `poly_initRing(vars, ord)` should be used.

off1 'allfac;

% We use parsing StandardFrom -> DIP routine from dp
% and then convert DIP to our own polynomial type.
load!-package 'dp;

% We use the term order comparators for exponent lists from dipoly
load!-package 'dipoly;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize polynomial ring with variables `vars` and term order `ord`
asserted procedure poly_initRing(vars: List, ord: Any): List;
   <<
      % initialize polynomials from cgb/dp for parsing input/output
      dip_init(vars, ord, nil);
      % initialize term order from dipoly/torder
      dipsortingmode(ord);
      torder({'list . vars, ord})
   >>;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% POLYNOMIAL INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard Polynomial constructor, forms a polynomial
% from a list of `Term`s (`Terms`) and list of `Coeff`s (`Coeffs`)
asserted inline procedure poly_Polynomial(ts: Terms, cfs: Coeffs): Polynomial;
   {'p, ts, cfs};

asserted inline procedure poly_getTerms(poly: Polynomial): Terms;
   cadr poly;

asserted inline procedure poly_getCoeffs(poly: Polynomial): Coeffs;
   caddr poly;

% Constructs a Polynomial from a SF `f`
asserted procedure poly_f2poly(f: SF): Polynomial;
   begin scalar fterms, fcoeffs, dpoly, ev, cf;
      integer deg;
    % construct dpoly..
      dpoly := dip_f2dip(f);
    % and parse it into our polynomial structure
      while dpoly do <<
         ev := pop(dpoly);
         cf := pop(dpoly);
      % also store the total degree
         % as a first element of each exponent list
         deg := for each x in ev sum x;
         push(deg, ev);
         push(ev, fterms);
         push(cf, fcoeffs)
      >>;
      return poly_Polynomial(reversip(fterms), reversip(fcoeffs))
   end;

% Constructs a Standard form for a Polynomial `f`
asserted procedure poly_poly2a(f: Polynomial): SF;
   begin scalar fterms, fcoeffs, dpoly, ev, cf;
      fterms  := poly_getTerms(f);
      fcoeffs := poly_getCoeffs(f);
      % First, we construct a `dip` polynomial from the input polynomial.
      % Then, we use `dip_2a` for the final conversion
      while fterms do <<
         % ev - a term in our Polynomial, represented as a list of integers
         ev := pop(fterms);
         % pop the first element, which stands for the total degree
         pop(ev);
         cf := pop(fcoeffs);
         push(ev, dpoly);
      push(cf, dpoly)
      >>;
      return dip_2a(reversip(dpoly))
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EXPONENT LISTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invariant: the first entry in the exponent list is the sum of subsequent entries

% Returns the first entry in the exponent list - the total degree
asserted inline procedure poly_totalDegExp(e1: List): Integer;
   car e1;

% Returns exponent list of zeros of the appropriate length
asserted inline procedure poly_zeroExp(): List;
   for x := 1:length(global!-dipvars!*) collect 0;

% Returns the elementwise sum of exponent lists e1, e2
asserted procedure poly_sumExp(e1: List, e2: List): List;
   if null e1 then
      nil
   else
      (car e1 #+ car e2) . poly_sumExp(cdr e1, cdr e2);

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
   if null e1 then
      nil
   else
      (car e1 #- car e2) . poly_subExp(cdr e1, cdr e2);

% Returns the elementwise maximum of exponent lists e1, e2
asserted procedure poly_elmaxExp(e1: List, e2: List): List;
   begin scalar ans;
      ans := poly_elmaxExp1(e1, e2);
      car ans := for each x in (cdr ans) sum x;
      return ans
   end;

asserted procedure poly_elmaxExp1(e1: List, e2: List): List;
   if null e1 then
      nil
   else
      max(car e1, car e2) . poly_elmaxExp1(cdr e1, cdr e2);

% Returns the elementwise maximum of exponent lists e1, e2
asserted procedure poly_elmaxExp(e1: List, e2: List): List;
   begin scalar w;
      ASSERT(not null e1 and not null e2);
      e1 := cdr e1;
      e2 := cdr e2;
      w := for each x in e1 collect max(x, pop e2);
      return (for each x in w sum x) . w
   end;

% Checks if e1 is elementwise less ore equal than e2
asserted procedure poly_divExp!?(e1: List, e2: List): Boolean;
   if null e1 then
      t
   else if car e1 #> car e2 then
      nil
   else
      poly_divExp!?(cdr e1, cdr e2);

% Checks if e1 is elementwise less ore equal than e2
asserted procedure poly_divExp!?(e1: List, e2: List): Boolean;
   begin scalar brk;
      % We deliberately also check the total degree.
      while not brk and e1 do
         brk := pop e1 #> pop e2;
      return brk
   end;

% Checks if exponent e1 is disjoint with e2
asserted procedure poly_disjExp!?(e1: List, e2: List): Boolean;
   poly_disjExp1(cdr e1, cdr e2);

asserted procedure poly_disjExp1(e1: List, e2: List): Boolean;
   if null e1 then
      t
   else if (car e1 #* car e2) #> 0 then
      nil
   else
      poly_disjExp1(cdr e1, cdr e2);

asserted procedure poly_disjExp!?(e1: List, e2: List): Boolean;
   begin scalar ok;
      e1 := cdr e1;
      e2 := cdr e2;
      ok = t;
      while ok and e1 do
         ok := pop e1 #= 0 or pop e2 #= 0;
      return ok
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparators for exponent lists.
% We mainly fall back to functions from dipoly package here

% Compares exponent lists e1, e2 w.r.t. lex term order,
% and returns e1 < e2
asserted procedure poly_cmpExpLex(e1: List, e2: List): Boolean;
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

% Compares exponent lists e1, e2 w.r.t. the current term order
asserted procedure poly_cmpExp(e1: List, e2: List): Boolean;
   if vdpsortmode!* eq 'lex then
      poly_cmpExpLex(e1, e2)
   else if vdpsortmode!* eq 'gradlex then
      poly_cmpExpGradLex(e1, e2)
   else if vdpsortmode!* eq 'revgradlex then
      poly_cmpExpRevGradLex(e1, e2)
   else
      evcomp(e1, e2);

% Compares exponent lists w.r.t. the total degree
asserted inline procedure poly_tdegCmpExp(e1: List, e2: List): Boolean;
   car e1 #< car e2;

% Checks that e1 = e2 elementwise
asserted procedure poly_eqExp!?(e1: List, e2: List): Boolean;
   evlexcomp(e1, e2) #= 0;

% Checks that e1 = e2 elementwise
asserted procedure poly_eqExp!?(e1: List, e2: List): Boolean;
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

% Returns zero polynomial, represented as
%   {'p, nil, nil}
% Ideally this should NEVER be called
asserted inline procedure poly_zero(): Polynomial;
   poly_Polynomial(nil, nil);

% Checks if `p` is zero
asserted inline procedure poly_iszero!?(p: Polynomial): Boolean;
   null poly_getTerms(p);

% Returns the tail of the polynomial `poly`.
% That is, `poly - lead(poly)`
asserted inline procedure poly_tail(poly: Polynomial): Polynomial;
   poly_Polynomial(poly_tailTerms(poly), poly_tailCoeffs(poly));

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
      denr(a) ./ numr(a);
   >>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% POLYNOMIAL LOW LEVEL OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This is the only section where polynomial coefficient arithmetic happens

% Returns s = gcoeff*fmult*f - fcoeff*gmult*g
%
% !! Assuming the leading terms of `gcoeff*fmult*f` and `fcoeff*gmult*g`
% mutually cancel each other
asserted inline procedure poly_paircombTail(f: Polynomial, fmult: Term,
                                                fcoeff: Coeff, g: Polynomial,
                                                gmult: Term,   gcoeff: Coeff): Polynomial;
   % Assuming the leading monomials of `gcoeff*fmult*f` and `fcoeff*gmult*g`
   % are mutually canceled, we can use the general reduction applied to
   % the tails of input polynomials
   <<
   % prin2t {"lc(f) =", poly_leadCoeff(f)};
   % prin2t {"lc(g) =", poly_leadCoeff(g)};
   poly_paircomb(poly_tail(f), fmult, fcoeff, poly_tail(g), gmult, gcoeff)
   >>;

% Returns s = gcoeff*fmult*f - fcoeff*gmult*g
asserted procedure poly_paircomb(f: Polynomial,  fmult: Term,
                                 fcoeff: Coeff, g: Polynomial,
                                 gmult: Term,   gcoeff: Coeff): Polynomial;
   begin scalar fterms, fcoeffs, gterms, gcoeffs, gmultcoeff, fmultcoeff,
                sterms, scoeffs, ft, gt, fc, gc, newc;
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
      % *in the same merge order as for the terms*.
      %
      fterms  := poly_getTerms(f);
      fcoeffs := poly_getCoeffs(f);
      gterms  := poly_getTerms(g);
      gcoeffs := poly_getCoeffs(g);
      gmultcoeff := poly_negCoeff(fcoeff);
      % fmultcoeff := gcoeff;

      if not poly_isoneCoeff!?(gcoeff) then
         rederr "Beda!";

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
         % prin2t {fterms, gterms};
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
      return poly_Polynomial(sterms, scoeffs)
   end;

% trst poly_paircomb;

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
      return poly_Polynomial(poly_getTerms(poly), newcoeffs)
   end;

% Constructs a new polynomial as a copy of `poly`
% with coefficients divided by the leading coeff of `poly`
asserted procedure poly_normalizeByLead(poly: Polynomial): Polynomial;
   begin scalar newcoeffs, mult1;
      mult1 := poly_invCoeff(poly_leadCoeff(poly));
      newcoeffs := for each cf in poly_getCoeffs(poly)
         collect poly_mulCoeff(cf, mult1);
      return poly_Polynomial(poly_getTerms(poly), newcoeffs)
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High level functions: Spolynomial computation and reduction

% Computes S-polynomial of f and g
asserted procedure poly_spoly(f: Polynomial, g: Polynomial): Polynomial;
   begin scalar e1, e2, elcm, mult1, mult2;
      e1 := poly_leadTerm(f);
      e2 := poly_leadTerm(g);
      elcm := poly_lcmTerm(e1, e2);
      mult1 := poly_divTerm(elcm, e2);
      mult2 := poly_divTerm(elcm, e1);
      % using poly_paircombTail since leading terms vanish
      return poly_paircombTail(f, mult2, poly_leadCoeff(f), g, mult1, poly_leadCoeff(g))
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
         f := poly_paircombTail(f, fmult, poly_leadCoeff(f), g, gmult, poly_leadCoeff(g));
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
   begin scalar updated, fterms, fcoeffs, glead, gcoef,
                fcoef, fterm, fmult, gmult;
      fterms  := poly_getTerms(f);
      fcoeffs := poly_getCoeffs(f);
      glead := poly_leadTerm(g);
      gcoef := poly_leadCoeff(g);
      while (not updated) and fterms do <<
         fterm := pop(fterms);
         fcoef := pop(fcoeffs);
         if poly_dividesTerm!?(glead, fterm) then <<
            fmult := poly_identityTerm();
            gmult := poly_divTerm(fterm, glead);
            f := poly_paircomb(f, fmult, fcoef, g, gmult, gcoef);
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
   poly_Polynomial(
      poly_getTerms(poly),
      for each c in poly_getCoeffs(poly) collect c . 1
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
    return poly_Polynomial(poly_getTerms(f), reversip(newcoeffs))
  end;

% Reduces coefficients of `f` by the given `prime`
asserted procedure poly_reduceCoeffs(f: Polynomial, prime: Integer): Polynomial;
  begin scalar fcoeffs, newcoeffs, c;
    % note that prime is not used here, since `modular!-number` works globally
    fcoeffs := poly_getCoeffs(f);
    while fcoeffs do <<
         c  := pop(fcoeffs);
         % ASSERT(denr(c) = 1);
         c := modular!-number(c);
         push(c, newcoeffs)
      >>;
      return poly_Polynomial(poly_getTerms(f), reversip(newcoeffs))
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional modular polynomial coefficient manipulation

% Reconstructs each coefficient of `poly` modulo the given prime
asserted procedure poly_reconstructCoeffs(poly: Polynomial,
                                          prime: Integer): Polynomial;
  begin scalar newcoeffs;
    newcoeffs := for each cf in poly_getCoeffs(poly)
      collect mod_reconstruction(cf, prime);
      return poly_Polynomial(poly_getTerms(poly), newcoeffs)
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
    return poly_Polynomial(poly_getTerms(polyaccum), reversip(newcoeffs))
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
