module f5poly;
% Polynomial interface module to be used in f5.
% The module provides procedures for basic operations with the `Polynomial` type.

% Polynomial `p` is stored as a list of 3 items:
%     {'p, Terms, Coeffs}
% Where `p` is a convenience tag,
%       `Terms` is a list of `Term`. Each `Term` is an exponent list of
%         non-negative integers of form
%              {totaldegree, pow1, pow2, ... pown}
%       `Coeffs` is a list of `Coeff`,
%         where `Coeff` can be either an SQ or an Integer
%   f5mod.red defines the following relevant functions on `Coeff`:
%     . mod_add(x, y)  -- addition x + y
%     . mod_div(x, y)  -- division x / y
%     . mod_inv(x)     -- inverse  x^(-1)
%     . mod_neg(x)     -- negation -x
%     . mod_mul(x, y)  -- product  x * y
%     . mod_iszero(x)  -- zero?    x
%
% `Terms` are ordered according to the current term order decreasingly,
%  and `Coeffs` are ordered respectively.
%
% For example, xy^2 + 3x is stored as
%   {'p, {{3, 1, 2}, {1, 1, 0}}, {1, 3}}
% if f5modular switch is on. It is stored as
%   {'p, {{3, 1, 2}, {1, 1, 0}}, {1 ./ 1, 3 ./ 1}}
% if f5modular is off (using SQ).
%
% Possible term orders are
%     lex, revgradlex

% The global polynomial ring should be initialized before constructing polynomials.
% To initialize the ring in variables `vars` and term order `ord` one
% should use `poly_initRing(vars, ord)`
% Initialization will set the following globals accordingly
fluid '(poly_ord!* poly_nvars!* poly_vars!*);

off1 'allfac;

% We use parsing StandardFrom -> DIP routine from dp
% and convert DIP to our own polynomial then.
% TODO: write our own?
load!-package 'dp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The current term order,
% `revgradlex` by default.
poly_ord!* := 'revgradlex;

% The current number of variables
poly_nvars!* := 0;

% The list of current variable identifiers
poly_vars!* := '(list);

% Initialize polynomial ring with variables `vars` and term order `ord`
asserted procedure poly_initRing(vars: List, ord: Id);
  <<
    poly_nvars!* := length(vars);
    poly_ord!*   := ord;
    poly_vars!*  := vars;
    dip_init(vars, ord, nil)
  >>;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% POLYNOMIAL INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard Polynomial ctor, construct a Polynomial
% from a list of Terms (exponent lists) and list of coefficients
asserted inline procedure poly_Polynomial(ts: Terms, cfs: Coeffs): Polynomial;
  {'p, ts, cfs};

asserted inline procedure poly_getTerms(poly: Polynomial): Terms;
  cadr poly;

asserted inline procedure poly_getCoeffs(poly: Polynomial): Coeffs;
  caddr poly;

% Construct a Polynomial from a SF
asserted procedure poly_f2poly(f: SF): Polynomial;
  begin scalar exps, coeffs, ans, dpoly, ev, cf;
        integer deg;
    % construct dpoly..
    dpoly := dip_f2dip(f);
    % and parse it into our polynomial
    while dpoly do <<
      ev := pop(dpoly);
      cf := pop(dpoly);
      % also store total degree in each exponent list
      deg := for each x in ev sum x;
      push(deg, ev);
      push(ev, exps);
      push(cf, coeffs)
    >>;
    exps := reversip(exps);
    coeffs := reversip(coeffs);
    return poly_Polynomial(exps, coeffs)
  end;

% Construct a Standard form a Polynomial
asserted procedure poly_poly2a(poly: Polynomial): SF;
  begin scalar ans, ts, coeffs, dpoly, ev, cf;
    ts   := poly_getTerms(poly);
    coeffs := poly_getCoeffs(poly);
    while ts do <<
      ev := pop(ts);
      pop(ev);
      cf := pop(coeffs);
      push(ev, dpoly);
      push(cf, dpoly)
    >>;
    dpoly := reversip(dpoly);
    return dip_2a(dpoly)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EXPONENT LISTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invariant: the first entry in the exponent list is the sum of subsequent entries

asserted inline procedure poly_totalDegExp(e1: List): Integer;
  car e1;

% return exponent list of zeros,
% Ideally this should rarely be called
asserted inline procedure poly_zeroExp(): List;
  for x := 0:poly_nvars!* collect 0;

% return the sum of exponent lists e1, e2
asserted procedure poly_sumExp(e1: List, e2: List): List;
  if null e1 then
    nil
  else
    (car e1 #+ car e2) . poly_sumExp(cdr e1, cdr e2);

% return the subtraction of exponent lists e1, e2
asserted procedure poly_subExp(e1: List, e2: List): List;
  if null e1 then
    nil
  else
    (car e1 #- car e2) . poly_subExp(cdr e1, cdr e2);

% return the elementwise maximum of exponent lists e1, e2
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

% check if exponent list e1 is elementwise greater than e2
asserted procedure poly_divExp!?(e1: List, e2: List);
  if null e1 then
    t
  else if car e1 #> car e2 then
    nil
  else
    poly_divExp!?(cdr e1, cdr e2);

% check if exponent e1 is disjoint with e2
asserted procedure poly_disjExp!?(e1: List, e2: List);
  poly_disjExp1(cdr e1, cdr e2);

asserted procedure poly_disjExp1(e1: List, e2: List);
  if null e1 then
    t
  else if (car e1 #* car e2) #> 0 then
    nil
  else
    poly_disjExp1(cdr e1, cdr e2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparators for exponent lists

% Compare exponent lists e1, e2 w.r.t. lex term order,
% returns e1 <ₗₑₓ e2
asserted procedure poly_cmpExpLex(e1: List, e2: List);
  begin integer ep1, ep2;
        scalar flag;
    flag := t;
    e1 := cdr e1;
    e2 := cdr e2;
    while e1 and flag do <<
      ep1 := pop(e1);
      ep2 := pop(e2);
      flag := (ep1 #= ep2)
    >>;
    return if flag then nil else ep1 #< ep2
  end;

% Compare exponent lists e1, e2 w.r.t. revgradlex term order
asserted inline procedure poly_cmpExpRevgradlex(e1: List, e2: List);
  if (car e1) #< (car e2) then
    t
  else if (car e1) #= (car e2) then
    poly_cmpExpRevLexHelper(cdr e1, cdr e2) #= 1
  else
    nil;

% This looks not good, to be changed at the first opportunity
asserted procedure poly_cmpExpRevLexHelper(e1: List, e2: List): Integer;
  begin integer ep1, ep2, cmp, rec;
        scalar last;
    ep1 := car e1;
    ep2 := car e2;
    cmp := if ep1 #> ep2 then
      1
    else if ep1 #= ep2 then
      2
    else
      3;
    last := null (cdr e1);
    if not last then
      rec := poly_cmpExpRevLexHelper(cdr e1, cdr e2)
    else
      rec := 2;
    return if ((not last) and rec #= 1) or (last and cmp #= 1) then
      1
    else if rec #= 2 then
      cmp
    else
      3
  end;

% Compare exponent lists e1, e2 w.r.t. the current order poly_ord!*
asserted procedure poly_cmpExp(e1: List, e2: List);
  if poly_ord!* eq 'lex then
    poly_cmpExpLex(e1, e2)
  else if poly_ord!* eq 'gradlex then
    poly_cmpExpGradlex(e1, e2)
  else
    poly_cmpExpRevgradlex(e1, e2);

% Compare exponent lists w.r.t. total degree
asserted inline procedure poly_tdegCmpExp(e1: List, e2: List);
  poly_totalDegExp(e1) #< poly_totalDegExp(e2);

% Check that e1 = e2 elementwise
asserted procedure poly_eqExp!?(e1: List, e2: List);
  if null e1 then
    t
  else if car e2 equal car e1 then
    poly_eqExp!?(cdr e1, cdr e2)
  else
    nil;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A term is essentially the same as exponent lists.
% This additional abstraction was made to handle exponent lists
% safely in other modules.

asserted inline procedure poly_identityTerm(): Term;
  poly_zeroExp();

asserted inline procedure poly_totalDegTerm(a: Term): Integer;
  poly_totalDegExp(a);

asserted inline procedure poly_mulTerm(a: Term, b: Term): Term;
  poly_sumExp(a, b);

asserted inline procedure poly_divTerm(a: Term, b: Term): Term;
  poly_subExp(a, b);

asserted inline procedure poly_lcmTerm(a: Term, b: Term): Term;
  poly_elmaxExp(a, b);

asserted inline procedure poly_cmpTerm(a: Term, b: Term): Term;
  poly_cmpExp(a, b);

asserted inline procedure poly_disjTerm!?(a: Term, b: Term);
  poly_disjExp!?(a, b);

asserted inline procedure poly_eqTerm!?(a: Term, b: Term);
  poly_eqExp!?(a, b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% POLYNOMIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns zero polynomial, internally represented as
%   {'p, nil, nil}
% Ideally this should NEVER be called
asserted inline procedure poly_zero(): Polynomial;
  poly_Polynomial(nil, nil);

% Checks if `p` is zero
asserted inline procedure poly_iszero!?(p: Polynomial);
  null poly_getTerms(p);

% Returns the tail of the polynomial `poly`.
% Essentially, `poly - lead(poly)`
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
%%%%%%%% POLYNOMIAL LOW LEVEL OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the only section where polynomial coefficient arithmetic happens

% returns s = fmult*f - C*gmult*g,
% where C is fcoeff/gcoeff,
%
% !! Assuming the leading terms of `f` and `g` do cancel each other
asserted inline procedure poly_paircombTail(f: Polynomial, fmult: Term,
                                    fcoeff: Coeff, g: Polynomial,
                                    gmult: Term, gcoeff: Coeff): Polynomial;
  poly_paircomb(poly_tail(f), fmult, fcoeff, poly_tail(g), gmult, gcoef);

% returns s = fmult*f - C*gmult*g,
% where C is fcoeff/gcoeff,
asserted procedure poly_paircomb(f: Polynomial, fmult: Term, fcoeff: Coeff,
                        g: Polynomial, gmult: Term, gcoeff: Coeff): Polynomial;
  begin scalar fterms, fcoeffs, gterms, gcoeffs, gmultcoeff, sterms, scoeffs,
                t1, t2, c1, c2, cf;
    fterms  := poly_getTerms(f);
    fcoeffs := poly_getCoeffs(f);
    gterms  := poly_getTerms(g);
    gcoeffs := poly_getCoeffs(g);
    gmultcoeff := mod_neg(mod_div(fcoeff, gcoeff));
    % Merge two sorted lists: fterms and gterms,
    % together with fcoeffs and gcoeffs
    while fterms and gterms do <<
      t1 := poly_mulTerm(car fterms, fmult);
      t2 := poly_mulTerm(car gterms, gmult);
      c1 := car fcoeffs;
      c2 := mod_mul(car gcoeffs, gmultcoeff);
      if poly_cmpTerm(t2, t1) then <<    % if t2 < t1
        push(t1, sterms);
        push(c1, scoeffs);
        pop(fterms);
        pop(fcoeffs)
      >> else if poly_eqTerm!?(t1, t2) then <<  % if t1 = t2
        cf := mod_add(c1, c2);
        if not mod_iszero!?(cf) then <<
          push(t1, sterms);
          push(cf, scoeffs)
        >>;
        pop(fterms);
        pop(fcoeffs);
        pop(gterms);
        pop(gcoeffs)
      >> else <<   % if t1 < t2
        push(t2, sterms);
        push(c2, scoeffs);
        pop(gterms);
        pop(gcoeffs)
      >>
    >>;
    % Merge what is left from fterms
    while fterms do <<
      push(poly_mulTerm(pop(fterms), fmult), sterms);
      push(pop(fcoeffs), scoeffs)
    >>;
    % Merge what is left from gterms
    while gterms do <<
      push(poly_mulTerm(pop(gterms), gmult), sterms);
      push(pop(gcoeffs), scoeffs)
    >>;
    return poly_Polynomial(reversip(sterms), reversip(scoeffs))
  end;

% Construct a new poly with all coefficients divided by the leading one
asserted procedure poly_normalize(poly: Polynomial): Polynomial;
  begin scalar newcoeffs, mult1, cf;
    mult1 := mod_inv(poly_leadCoeff(poly));
    newcoeffs := for each cf in poly
      collect mod_mul(cf, mult1);
    return poly_Polynomial(poly_getTerms(poly), newcoeffs)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some high level functions: Spolynomial computation and reduction

% Compute S-polynomial of f and g
asserted procedure core_spoly(f: Polynomial, g: Polynomial): Polynomial;
  begin scalar e1, e2, elcm, mult1, mult2;
    e1 := poly_leadTerm(f);
    e2 := poly_leadTerm(g);
    elcm := poly_lcmTerm(e1, e2);
    mult1 := poly_divTerm(elcm, e2);
    mult2 := poly_divTerm(elcm, e1);
    return poly_paircombTail(f, mult2, poly_leadCoeff(f), g, mult1, poly_leadCoeff(g))
  end;

% Tries to top-reduce f with g (only the leading term of f, only once).
% Two cases are possible:
% 1. the leading term of g divides the leading term of f. Then we perform
%    top-reduction and return true flag and the result of reduction.
% 2. top-reduction is not possible, we return false flag and f
asserted procedure poly_tryTopReductionStep(f: Polynomial,
                                              g: Polynomial): Polynomial;
  begin scalar glead, flead, fmult, gmult, updated;
    glead := poly_leadTerm(g);
    flead := poly_leadTerm(f);
    if poly_divTerm!?(ge, fe) then <<
      fmult := poly_identityTerm();
      gmult := poly_divTerm(fe, ge);
      f := poly_paircombTail(f, fmult, poly_leadCoeff(f), g, gmult, poly_leadCoeff(g));
      updated := t;
    >>;
    return updated . f
  end;

% Tries to reduce f with g. Two cases are possible:
% 1. the leading term of g divides some term of f. Then we perform
%    reduction and return true flag and the result of reduction.
% 2. reduction is not possible, we return false flag and f
asserted procedure poly_tryReductionStep(f: Polynomial,
                                          g: Polynomial): Polynomial;
  begin scalar updated, fterms, fcoeffs, glead, gcoef,
                fcoef, fterm, fmult, gmult, updatedToreturn;
    fterms  := poly_getTerms(f);
    fcoeffs := poly_getCoeffs(f);
    glead := poly_leadTerm(g);
    gcoef := poly_leadCoeff(g);
    while (not updated) and fterms do <<
      fterm := pop(fexps);
      fcoef := pop(fcoeffs);
      if poly_divTerm!?(glead, fterm) then <<
        fmult := poly_identityTerm();
        gmult := poly_divTerm(fterm, glead);
        f := poly_paircomb(f, fmult, fc, g, gmult, gc);
        updated := t
      >>
    >>;
    return updated . f
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional polynomial coefficient manipulations

% Return the lcm of denominators of coefficients of `f`
asserted procedure poly_commonDenominator(f: Polynomial): Integer;
  begin scalar fcoeffs;
        integer den;
    den := 1;
    fcoeffs := poly_getCoeffs(f);
    while fcoeffs do <<
      den := lcm(den, denr(pop(fcoeffs))
    >>;
    return den
  end;

% Construct a new polynomial `f` in the following way:
%   f * inv(poly_commonDenominator(f))
asserted procedure poly_scaleDenominators(f: Polynomial): Polynomial;
  begin scalar fcoeffs, newcoeffs;
        integer den;
    den := poly_commonDenominator(f);
    fcoeffs := poly_getCoeffs(f);
    while fcoeffs do <<
      c := pop(coeffs);
      push(numr(c) * (den / denr(c)), newcoeffs)
    >>;
    return poly_Polynomial(poly_getTerms(f), reversip(newcoeffs))
  end;

% Reduce coefficients of `f` by the given `prime`
asserted procedure poly_reduceCoeffs(f: Polynomial, prime: Integer): Polynomial;
  begin scalar fcoeffs, newcoeffs, c;
    % note that prime is not used here as `modular!-number` works globally
    fcoeffs := poly_getCoeffs(poly);
    while fcoeffs do <<
         c  := pop(coeffs);
         % ASSERT(denr(c) = 1);
         c := modular!-number(c);
         push(c, newcoeffs)
      >>;
      return poly_Polynomial(poly_getTerms(poly), reversip(newcoeffs))
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional modular polynomial coefficient manipulation

% Reconstruct each coefficient of `poly` modulo the given prime
asserted procedure poly_reconstructCoeffs(poly: Polynomial,
                                          prime: Integer): Polynomial;
  begin scalar newcoeffs;
    newcoeffs := for each cf in poly_getCoeffs(poly)
      collect mod_reconstruction(cf, prime);
    return poly_Polynomial(poly_getTerms(poly), newcoeffs)
  end;

% Apply CRT to coefficients of (polyaccum mod modulo) and (polycomp mod prime)
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
asserted inline procedure poly_cmpPolyLead(poly1: Polynomial, poly2: Polynomial);
  poly_cmpTerm(poly_leadTerm(poly1), poly_leadTerm(poly2));

% Returns t if
% the total degree of lead term of poly1 < the total degree of lead term poly2,
% Otherwise if these degrees are equal compare term with the current order,
% Otherwise nil.
asserted procedure poly_leadTotalDegreeCmp(poly1: Polynomial, poly2: Polynomial);
  begin integer t1, t2;
    t1 := poly_totalDegTerm(poly_leadTerm(poly1));
    t2 := poly_totalDegTerm(poly_leadTerm(poly2));
    return if t1 #= t2 then
      poly_cmpPolyLead(poly1, poly2)
    else
      t1 #< t2
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endmodule;

end;  % of file
