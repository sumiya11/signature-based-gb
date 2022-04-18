module f5poly;

% Polynomial interface module to be used in f5 computation.
% The module provides functions for basic operations with `Polynomial` type

% do we need this hmm?
off1 'allfac;

% we use parsing StandardFrom -> Polynomial routine from dp
% TODO: write our own?
load_package 'dp;

fluid '(poly_ord!* poly_nvars!* poly_vars!*);

% Polynomial `p` is stored internally as a list of 3 items:
%     {'p, monomials, coeffs}
% Where `p` is a convenience tag,
%         `monomials` is an array of monomials stored as exponent vectors,
%         `coeffs` is an array of coefficients
% with `monomials` and `coeffs` ordered from lead to tail. Exponent vectors
% are stored internally as integer lists:
%		{totaldegree, pow1, pow2, ... pown}
%
% For example, xy^2 + 3x is stored as
%   {'p, {{3, 1, 2}, (1, 1, 0)}, {1, 3}}
%
%	Possible monomial orderings are
%     lex, revgradlex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the current monomial ordering
poly_ord!* := 'revgradlex;

% the current number of variables
poly_nvars!* := 0;

% the current variables
poly_vars!* := '(list);

% initialize polynomial ring with variables `vars` and monomial ordering `ord`
asserted procedure poly_initRing(vars, ord);
  begin;
    poly_nvars!* := length(vars);
    poly_ord!*   := ord;
    poly_vars!*  := vars;

		dip_init(vars, ord, nil)
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inline procedure poly_getExps(poly: Polynomial): List;
	cadr poly;

inline procedure poly_getCoeffs(poly: Polynomial): List;
	caddr poly;

% Standard ctor from list of exponent vectors and list of coefficients
inline procedure poly_init(exps: List, coeffs: List): Polynomial;
	'p . exps . coeffs . nil;

% Standard form -> Polynomial
asserted procedure poly_f2poly(f);
	begin scalar exps, coeffs, ans, dpoly, ev, cf, deg;
		dpoly := dip_f2dip(f);
		while dpoly do <<
			ev . dpoly := dpoly;
			cf . dpoly := dpoly;

			deg := for each x in ev sum x;
			ev := deg . ev;
			exps   := ev . exps;
			coeffs := cf . coeffs;
		>>;
		exps := reversip(exps);
		coeffs := reversip(coeffs);
		return poly_init(exps, coeffs)
	end;

% Polynomial -> Standard form
asserted procedure poly_poly2a(poly);
	begin scalar ans, exps, coeffs, dpoly,
                ev, cf;
		exps := poly_getExps(poly);
		coeffs := poly_getCoeffs(poly);
		dpoly := nil;
		while exps do <<
			ev . exps   := exps;
			ev := cdr ev;
			cf . coeffs := coeffs;
			dpoly := cf . ev . dpoly;
		>>;
		dpoly := reversip(dpoly);
		return dip_2a(dpoly)
	end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EXPONENT VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invariant: the first entry in the vector is the sum of subsequent entries

% return exponent vector of zeros,
% Ideally this should rarely be called
inline procedure poly_zeroExp(): List;
	for x := 0:poly_nvars!* collect 0;

% return sum of exponent vectors e1, e2
asserted procedure poly_sumExp(e1: List, e2: List): List;
	if null e1 then
		nil
	else
		(car e1 + car e2) . poly_sumExp(cdr e1, cdr e2);

% return difference of exponent vectors e1, e2
asserted procedure poly_subExp(e1: List, e2: List): List;
	if null e1 then
		nil
	else
		(car e1 - car e2) . poly_subExp(cdr e1, cdr e2);

% return elementwise maximum of exponent vectors e1, e2
asserted procedure poly_lcmExp(e1: List, e2: List): List;
	begin scalar ans;
		ans := poly_lcmExp1(e1, e2);
		car ans := 0;
		car ans := for each x in ans sum x;
		return ans
	end;

asserted procedure poly_lcmExp1(e1: List, e2: List): List;
	if null e1 then
		nil
	else
		max(car e1, car e2) . poly_lcmExp1(cdr e1, cdr e2);

% check if exponent e1 divides e2
asserted procedure poly_divExp!?(e1: List, e2: List);
	if null e1 then
		t
	else if car e1 > car e2 then
		nil
	else
		poly_divExp!?(cdr e1, cdr e2);

% check if exponent e1 is disjoint with e2
asserted procedure poly_disjExp!?(e1: List, e2: List);
	poly_disjExp1(cdr e1, cdr e2);

asserted procedure poly_disjExp1(e1: List, e2: List);
	if null e1 then
		t
	else if (car e1 * car e2) > 0 then
		nil
	else
		poly_disjExp1(cdr e1, cdr e2);

% comparator for exponent vectors e1, e2 w.r.t. lex monomial ordering
% returns e1 <ₗₑₓ e2
asserted procedure poly_cmpExpLex(e1: List, e2: List);
	begin integer ep1, ep2;
        scalar flag;
		flag := t;
		e1 := cdr e1;
		e2 := cdr e2;
		while e1 and flag do <<
			ep1 . e1 := e1;
			ep2 . e2 := e2;
			flag := (ep1 = ep2);
		>>;
		return if flag then nil else ep1 < ep2
	end;

% comparator for exponent vectors e1, e2 w.r.t. gradlex monomial ordering
asserted procedure poly_cmpExpGradlex(e1: List, e2: List);
	begin integer ep1, ep2;
        scalar flag;
    flag := t;
    e1 := cdr e1;
    e2 := cdr e2;
    while e1 and flag do <<
      ep1 . e1 := e1;
      ep2 . e2 := e2;
      flag := (ep1 = ep2);
    >>;
    return if flag then nil else ep1 < ep2
  end;

% comparator for exponent vectors e1, e2 w.r.t. revgradlex monomial ordering
asserted procedure poly_cmpExpRevgradlex(e1: List, e2: List);
  begin integer ep1, ep2;
			ep1 := car e1;
			ep2 := car e2;
			return if ep1 #< ep2 then
				t
			else if ep1 #= ep2 then
				poly_cmpExpRevLex(cdr e1, cdr e2)
			else
				nil
  end;

% comparator for exponent vectors e1, e2 w.r.t. revlex monomial ordering
asserted procedure poly_cmpExpRevLex(e1: List, e2: List);
  poly_cmpExpRevLexHelper(e1, e2) #= 1;

asserted procedure poly_cmpExpRevLexHelper(e1: List, e2: List);
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


% comparator for exponent vectors e1, e2
asserted procedure poly_cmpExp(e1: List, e2: List);
	if poly_ord!* eq 'lex then
		poly_cmpExpLex(e1, e2)
	else if poly_ord!* eq 'gradlex then
		poly_cmpExpGradlex(e1, e2)
	else
		poly_cmpExpRevgradlex(e1, e2);

% checks that e1 = e2 elementwise
asserted procedure poly_eqExp!?(e1: List, e2: List);
	if null e1 then
		t
	else if car e2 equal car e1 then
		poly_eqExp!?(cdr e1, cdr e2)
	else
		nil;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% returns zero polynomial, internally represented as
%   ('p, nil, nil)
% Ideally this should NEVER be called
inline procedure poly_zero(): Polynomial;
	poly_init({}, {});

% if p is zero
inline procedure poly_iszero!?(p);
	null poly_getExps(p);

% returns the tail of poly
% essentially, returns poly-lead(poly)
asserted procedure poly_tail(poly: Polynomial): Polynomial;
	poly_init(poly_tailExps(poly), poly_tailCoeffs(poly));

% returns the leading exponent of poly
inline procedure poly_leadExp(poly: Polynomial): List;
  car poly_getExps(poly);

% returns the leading coefficient of poly
asserted procedure poly_leadCoeff(poly: Polynomial);
	car poly_getCoeffs(poly);

% returns the tail exponents of poly
inline procedure poly_tailExps(poly: Polynomial): List;
	cdr poly_getExps(poly);

% returns the tail coefficients of poly
asserted procedure poly_tailCoeffs(poly: Polynomial);
	cdr poly_getCoeffs(poly);

% checks if polynomials leading terms are disjoint
asserted procedure poly_disjLead!?(p1: Polynomial, p2: Polynomial);
	poly_disjExp!?(poly_leadExp(p1), poly_leadExp(p2));

% returns length of poly, i.e., the number of terms
asserted procedure poly_length(poly: Polynomial): Integer;
  length(poly_getExps(poly));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% POLYNOMIAL ARITHMETIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Contract: this is the only section that uses polynomial coefficient arithmetic

% returns s = multf*f - C*multg*g
% where C is -fcoeff/gcoeff
asserted procedure poly_paircomb(f, fmult, fcoeff, g, gmult, gcoeff): Polynomial;
	begin scalar fexps, fcoeffs, gexps, gcoeffs, gmultcoeff, sexps, scoeffs,
                e1, e2, c1, c2, c;
    sexps   := {};
    scoeffs := {};

    fexps   := poly_getExps(f);
    fcoeffs := poly_getCoeffs(f);
    gexps   := poly_getExps(g);
    gcoeffs := poly_getCoeffs(g);

    gmultcoeff := mod_neg(mod_div(fcoeff, gcoeff));

    while fexps and gexps do <<
      e1 := poly_sumExp(car fexps, fmult);
      e2 := poly_sumExp(car gexps, gmult);

      c1 := car fcoeffs;
      c2 := car gcoeffs;
      c2 := mod_mul(c2, gmultcoeff);

      if poly_cmpExp(e2, e1) then <<
        sexps   := e1 . sexps;
        scoeffs := c1 . scoeffs;
        fexps   := cdr fexps;
        fcoeffs := cdr fcoeffs
			>> else if poly_eqExp!?(e1, e2) then <<
				c := mod_add(c1, c2);
				if not mod_iszero!?(c) then <<
          sexps   := e1 . sexps;
          scoeffs := c . scoeffs
				>>;
        fexps   := cdr fexps;
        fcoeffs := cdr fcoeffs;
        gexps   := cdr gexps;
        gcoeffs := cdr gcoeffs
			>> else <<
        sexps   := e2 . sexps;
        scoeffs := c2 . scoeffs;
        gexps   := cdr gexps;
        gcoeffs := cdr gcoeffs
      >>
    >>;

		while fexps do <<
		  e1 := poly_sumExp(car fexps, fmult);
			c1 := car fcoeffs;
      sexps   := e1 . sexps;
      scoeffs := c1 . scoeffs;
      fexps   := cdr fexps;
      fcoeffs := cdr fcoeffs
		>>;

		while gexps do <<
			e2 := poly_sumExp(car gexps, gmult);
			c2 := car gcoeffs;
			c2 := mod_mul(c2, gmultcoeff);
      sexps   := e2 . sexps;
      scoeffs := c2 . scoeffs;
      gexps   := cdr gexps;
      gcoeffs := cdr gcoeffs
		>>;

    sexps := reversip(sexps);
    scoeffs := reversip(scoeffs);

    return poly_init(sexps, scoeffs)
	end;

% divide all coefficients by the leading one
asserted procedure poly_normalize(poly: Polynomial): Polynomial;
  begin scalar mult1, sexps, scoeffs, exps, coeffs, ex, cf;
    mult1 := mod_inv(poly_leadCoeff(poly));
    sexps   := {};
    scoeffs := {};
    exps    := poly_getExps(poly);
    coeffs  := poly_getCoeffs(poly);
    while exps do <<
      ex . exps   := exps;
      cf . coeffs := coeffs;
      cf := mod_mul(cf, mult1);

      sexps   := ex . sexps;
      scoeffs := cf . scoeffs;
    >>;

    sexps := reversip(sexps);
    scoeffs := reversip(scoeffs);

    return poly_init(sexps, scoeffs)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure poly_commonDenominator(f);
  begin scalar den, coeffs, c;
    den := 1;
    coeffs := poly_getCoeffs(f);
    while coeffs do <<
      c . coeffs := coeffs;
      den := lcm(den, denr(c))
    >>;
    return den
  end;

asserted procedure poly_scaleDenominatorsInplace(f);
  begin scalar den, coeffs, c;
    den := poly_commonDenominator(f);
    coeffs := poly_getCoeffs(f);
    while coeffs do <<
      c := car coeffs;
      c := numr(c) * (den / denr(c));
      car coeffs := c;
      coeffs := cdr coeffs
    >>;
    return f
  end;

asserted procedure poly_scaleDenominators(f);
  begin scalar copyf;
    copyf := copy(f);
    return poly_scaleDenominatorsInplace(copyf)
  end;

% reduce coefficients of poly by the given prime and return new polynomial
asserted procedure poly_reduceCoeffs(poly: Polynomial, prime): Polynomial;
   begin scalar coeffs, ansCoeffs, c;
		coeffs := poly_getCoeffs(poly);
		ansCoeffs := nil;
      while coeffs do <<
         c  . coeffs := coeffs;

         % ASSERT(denr(c) = 1);

         c := modular!-number(c);
         ansCoeffs :=  c . ansCoeffs
      >>;

      return poly_init(poly_getExps(poly), reversip(ansCoeffs))
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reconstruct coefficients of poly by the given prime and return new polynomial
asserted procedure poly_reconstructCoeffs(poly: Polynomial, prime): Polynomial;
  begin scalar coeffs, c, ansCoeffs;
    coeffs := poly_getCoeffs(poly);
		ansCoeffs := nil;
    while coeffs do <<
      c  . coeffs := coeffs;
      c := mod_reconstruction(c, prime);
      ansCoeffs := c . ansCoeffs
    >>;

		return poly_init(poly_getExps(poly), reversip(ansCoeffs))
   end;

% Apply CRT to (polyaccum mod modulo) and (polycomp mod prime)
% to obtain new polynomial over modulo*prime
asserted procedure poly_crtCoeffs(polyaccum, modulo, polycomp, prime): Polynomial;
  begin scalar coeffsaccum, coeffscomp, ansCoeffs, ca, cc, c;
    coeffsaccum := poly_getCoeffs(polyaccum);
    coeffscomp  := poly_getCoeffs(polycomp);
		ansCoeffs := nil;
    while coeffsaccum do <<
      ca . coeffsaccum := coeffsaccum;
      cc . coeffscomp  := coeffscomp;
      c := mod_crt(ca, modulo, cc, prime);
      ansCoeffs := c . ansCoeffs
    >>;
    return poly_init(poly_getExps(polyaccum), reversip(ansCoeffs))
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial sorting ad-hoc

% return true if lead(poly1) < lead(poly2)
asserted procedure poly_cmpPolyLead(poly1, poly2);
 	poly_cmpExp(poly_leadExp(poly1), poly_leadExp(poly2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trst poly_f2poly;
% trst poly_poly2a;

% trst poly_leadCmp;
% trst poly_lcmExp;
% trst poly_normalize;

% trst poly_paircomb;
% trst poly_unsafePaircombInplace;
% trst poly_cmpExpRevgradlex;
% trst poly_cmpExpRevlex;
% trst poly_cmpExpLex;
% trst poly_cmpExpRevLexHelper;

% trst poly_reduceCoeffs;
% trst poly_crtCoeffs;
% trst poly_reconstructCoeffs;
% trst poly_scaleDenominatorsInplace;

endmodule;

end;  % of file
