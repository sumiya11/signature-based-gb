
% Main module of the f5 package
module f5;

% The f5 module provides Groebner basis computation algorithm `f5` implementation
% based on the Fougere's F5 algorithm

% The interface consists of the `f5` function for computing Groebner bases.
%
% - The signature of `f5` is `f5(polynomials, variables, ordering)` where
%   . `polynomials` is the list of ideal generators;
%   . `variables` is the list of kernels to compute basis w.r.t.,
%     kernels not present in `variables` are treated as constants;
%   . `ordering` is the monomial ordering to compute the basis in,
%      possible options are `lex`, `revgradlex`;
%
% - The algorithm uses modular computations and, thus, is randomized.
%   Obtained result will be correct with high probability.
%   If *guaranteed correctness* is needed, use `on f5certify`;
%
% - If using modular computation is not desirable, set `off f5modular`.

create!-package('(f5 f5lp f5poly f5mod f5primes), nil);

load_package 'assert;
on1 'assert;

% If f5 should certify the correctness of result.
% False by default, meaning that the algorithm is randomized
% and may output incorrect answer with a small probability (~1/2^22)
switch f5certify;
off1 'f5certify;

% Maybe place other switches also here?.. 

put('f5, 'psopfn, 'f5_groebner);

% interface implemented in f5poly.red
struct Polynomial;

% interface implemented in f5lp.red
struct LabeledPolynomial;

% interface implemented in f5primes.red
struct Primetracker;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The main function to parse input arguments and call the f5 routine
asserted procedure f5_groebner(u: List): List;
   begin scalar inputBasis, variables, sortMode,
                inputModule, outputModule;
      if null u then
         f5_argumentError();
      inputBasis := reval pop u;
      if not (pop inputBasis eq 'list) then
         f5_argumentError();
      variables := reval pop u;
      if not (pop variables eq 'list) then
         f5_argumentError();
      sortMode := pop u;
      if not null u then
         f5_argumentError();

      % initialize ground polynomial ring
      poly_initRing(variables, sortMode);

      inputBasis := for each f in inputBasis collect
         poly_f2poly numr simp f;
      % construct module basis from input polynomials
      inputModule := f5_constructModule(inputBasis);

      if !*f5modular then
         outputModule := f5_groebnerModular1(inputModule)
      else
         outputModule := f5_groebner1(inputModule);

      outputModule := 'list . for each f in outputModule collect
                        poly_poly2a lp_eval f;
      return outputModule
   end;

% Argument error
asserted procedure f5_argumentError();
   rederr "usage: f5(polynomials: List, variables: List, sortmode: Id)";

% All of the functions bellow work with module elements
% instead of accessing polynomials directly
%
% This allows us to stick with LabeledPolynomials
% and never perform conversions LabeledPolynomials <-> Polynomial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR CORRECTNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verify that the given `recBasis` is indeed the Groebner basis
% of the original input set `intBasis`
% By default, probabilistic correctness check is used
asserted procedure f5_correctnessCheck(pt, intBasis, recBasis);
   if not f5_heuristicCorrectnessCheck(recBasis) then
    nil
  else if not f5_randomizedCorrectnessCheck(pt, intBasis, recBasis) then
    nil
  else if !*f5certify then
    f5_guaranteedCorrectnessCheck(reconstructedBasis)
  else
    t;

% Heuristic correctness check based on the bitsize of coefficients
% NOT TODO
asserted procedure f5_heuristicCorrectnessCheck(reconstructedBasis: List);
  t;

% Randomized correctness check, checks that
%   . <intBasis> in <recBasis> as ideals
%   . recBasis is a Groebner basis
% modulo chosen `reliablePrime`
asserted procedure f5_randomizedCorrectnessCheck(pt, intBasis, recBasis);
  begin integer reliablePrime;
        scalar intBasisReduced, recBasisScaled, recBasisReduced,
                isgroebner;
    reliablePrime := primes_nextReliablePrime(pt, intBasis);

    intBasisReduced := f5_modularReduction(intBasis, reliablePrime);

    recBasisScaled := f5_scaleDenominators(recBasis);
    recBasisReduced := f5_modularReduction(recBasisScaled, reliablePrime);

    isgroebner := f5_isGroebner1(recBasisReduced);

    return if isgroebner then
      f5_checkIdealInclusion1(recBasisReduced, intBasisReduced)
    else
      nil
  end;

% Guaranteed correctness check. Same as f5_randomizedCorrectnessCheck,
% but all computations are carried in the original ground ring
% NOT TODO
asserted procedure f5_guaranteedCorrectnessCheck(reconstructedBasis: List);
  t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All functions below that take a basis and transform its coefficients
% preserve the order of generators in the original basis

% Given the set of generators `inputBasis` with integer coefficients construct
% a new basis with coefficients reduced modulo `prime`
asserted procedure f5_modularReduction(inputBasis: List, prime: Integer): List;
  begin scalar prime, ans, poly;
    while inputBasis do <<
      poly . inputBasis := inputBasis;
      ans := lp_reduceCoeffs(poly, prime) . ans
    >>;
    return reversip(ans)
  end;

% Given the basis `inputBasis` computed modulo some prime,
% reconstructs each coefficient into the rational field.
% Rationals are represented as Standard Quotients (pairs) of integers
% This function is assumed to be called only
% when coeffs of `inputBasis` basis are integers
asserted procedure f5_rationalReconstruction(
                                      inputBasis: List,
                                      pt: Primetracker): List;
   begin scalar poly, ans;
          integer prime;
      prime := primes_getAccumModulo(pt);
      while inputBasis do <<
         poly . inputBasis := inputBasis;
         ans := lp_reconstructCoeffs(poly, prime) . ans
      >>;
      return reversip(ans)
   end;

% Given the basis `accumBasis` computed modulo some `modulo`,
% and the basis `computedBasis` computed modulo some `prime`,
% constructs a new basis modulo `modulo × prime`.
% This function is assumed to be called only
% when coeffs of input bases are integers
asserted procedure f5_crtReconstruction(
                                  accumBasis: List,
                                  computedBasis: List,
                                  pt: Primetracker);
   begin scalar ans, polyaccum, polycomp;
          integer modulo, prime;
      modulo := primes_getAccumModulo(pt);
      prime  := primes_getLuckyPrime(pt);
      if modulo #= 1 then
        % essentially, accumBasis is garbage here
        ans := computedBasis
      else <<
        while accumBasis do <<
          polyaccum . accumBasis := accumBasis;
          polycomp  . computedBasis := computedBasis;
          ans := lp_crtCoeffs(polyaccum, modulo, polycomp, prime) . ans
        >>;
        ans := reversip(ans)
      >>;
      primes_setAccumModulo(primetracker, modulo * prime);
      return ans
   end;

% Given the set of generators `inputBasis` with rational numbers
% as polynomial coefficients, scale INPLACE
% each generator by the common denominator to obtain integer coefficients.
% Assumed to be called only for numeric rational number input.
asserted procedure f5_scaleDenominatorsInplace(inputBasis: List): List;
   begin scalar tmp, poly;
      tmp := inputBasis;
      while tmp do <<
         poly := car tmp;
         car tmp := lp_scaleDenominatorsInplace(poly);
         tmp := cdr tmp
      >>;
      return inputBasis
   end;

% Given the set of generators `inputBasis` with rational numbers
% as polynomial coefficients, scale
% each generator by the common denominator to obtain integer coefficients.
% Assumed to be called only for numeric rational number input.
asserted procedure f5_scaleDenominators(inputBasis: List): List;
  begin scalar ans, poly;
    while inputBasis do <<
      poly . inputBasis := inputBasis;
      ans := lp_scaleDenominators(poly) . ans
    >>;
    return reversip(ans)
  end;

% Main function for modular F5 Groebner basis computation
%
% !! Takes module elements in input and returns module elements in output
asserted procedure f5_groebnerModular1(inputBasis: List): List;
  begin scalar integerBasis, reducedBasis, computedBasis,
                reconstructedBasis, accumBasis,
                correctness, primetracker;
        integer iter, prime;

      % contains current lucky prime for reduction and accumulates
      % the product of all previous ones
      primetracker := primes_Primetracker();

      % scale coefficients to integers
      integerBasis := f5_scaleDenominators(inputBasis);
      % now all denominators are 1 and coefficients are actually "big" integers

      % this basis will store accumulated coefficients
      % from all previous CRT calls
      accumBasis := nil;

      % so, while the correctness check failes
      while not correctness do <<
         % select next lucky prime
         prime := primes_nextLuckyPrime(primetracker, integerBasis);
         ASSERT(primep prime);

         % reduce the basis w.r.t the prime
         reducedBasis := f5_modularReduction(integerBasis, prime);
         % now all coefficients in reducedBasis should be small integers

         % compute the basis in the finite field
         computedBasis := f5_groebner1(reducedBasis);

         % CRT (prevprime, prime) --> (prevprime*prime)
         accumBasis := f5_crtReconstruction(accumBasis, computedBasis, primetracker);

         % reconstruct modulo modulo
         reconstructedBasis := f5_rationalReconstruction(accumBasis, primetracker);

         correctness := f5_correctnessCheck(primetracker,
                                            integerBasis,
                                            reconstructedBasis);
         iter := iter + 1
      >>;

      return reconstructedBasis
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERIC F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the given set of polynomials `inputBasis` construct
% the standard basis of unit vectors corresponding to F
% , i.e.,
% {f, g}  --> {(1, 0), (0, 1)}
asserted procedure f5_constructModule(inputBasis: List): List;
  begin scalar poly, elem, outputModule;
        integer i;
    i := 1;
    while inputBasis do <<
      poly . inputBasis := inputBasis;
      elem := lp_LabeledPolynomial1(poly, i);
      outputModule := elem . outputModule;
      i := i + 1;
    >>;
    return outputModule
  end;

% for the given set of LabeledPolynomials `inputBasis` construct
% the set of all possible principal syzygies
%
% TODO: idea for optimization - sort inputBasis first
asserted procedure f5_constructSyzygies(inputBasis: List): List;
  begin scalar syzygies, inputSlice, p1, p2;
        integer i, n;
    n := length(inputBasis);
    for i := 1:n do <<
      inputSlice := cdr inputBasis;
      for j := i+1:n do <<
        p1 := car inputBasis;
        p2 := car inputSlice;
        syzygies := lp_principalSyzygy(p1, p2) . syzygies
      >>;
      inputBasis := cdr inputBasis
    >>;
    return syzygies
  end;

% for the given set of LabeledPolynomials `inputBasis` construct
% the set of all possible S-polynomials, while also applying 1st criterion
asserted procedure f5_constructSpolys(inputBasis: List): List;
  begin scalar spolys, inputSlice, p1, p2;
        integer i, n;
    n := length(inputBasis);
    for i := 1:n do <<
      inputSlice := cdr inputBasis;
      for j := i+1:n do <<
        p1 := car inputBasis;
        p2 := car inputSlice;
        if not poly_disjLead!?(lp_eval(p1), lp_eval(p2)) then
            spolys := lp_spoly(p1, p2) . spolys;
        inputSlice := cdr inputSlice
      >>;
      inputBasis := cdr inputBasis
    >>;
    return spolys
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if computing normal form of `p` is redundant
% w.r.t. already known `syzygies`
%
% TODO: idea for optimization - sort syzygies to perform less division checks
asserted procedure f5_syzygyCriterion(p: LabeledPolynomial, syzygies: List);
   begin scalar flag, syz, sgn;
      if null syzygies then
         flag := nil
      else <<
         syz := lp_sgn(car syzygies);
         sgn := lp_sgn(p);

         if lp_sgnDivides(syz, sgn) then
            flag := t
         else
            flag := f5_syzygyCriterion(p, cdr syzygies)
      >>;
      return flag
   end;

% Select the polynomial with the lowest signature in `spolys` and return
% it together with its index in `spolys`
%
% TODO: idea for optimization - store spolys in a balanced tree
asserted procedure f5_selectNext(spolys: List): List;
  begin integer i, idx;
        scalar sgn, elem, s, p;
    i := 1;
    s := lp_sgn(car spolys); % nonempty
    while spolys do <<
      p . spolys := spolys;
      if (i equal 1) or lp_potCmpSignature(s, sgn) then <<
        idx  := i;
        elem := p;
        sgn  := s
      >>;
      i := i + 1
    >>;
    return idx . elem
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given the set of generators `basis` apply the autoreduction algorithm
% to interreduce generators w.r.t. each other
asserted procedure f5_interreduceBasis(basis: List): List;
  begin scalar isregular, updated, poly, newbasis, flag,
                reducers;
    isregular := nil;
    updated := t;
    while updated do <<
      updated := nil;
      newbasis := nil;
      while basis do <<
        poly . basis := basis;
        % TODO: do this inplace
        reducers := append(newbasis, basis);
        flag . poly := lp_tryReduce(poly, reducers, isregular);
        updated := flag or updated;
        if not lp_iszero!?(poly) then <<
          newbasis := poly . newbasis
        >>
      >>;
      basis := newbasis
    >>;
    return basis
  end;

% normalize each generator in the given `basis` by the leading coefficient
asserted procedure f5_normalizeBasis(basis: List): List;
  for each x in basis collect lp_normalize(x);

% given a Groebner `basis`
% returns the unique Groebner basis of the corresponding ideal <basis>
%
% Output invariants:
%  . the basis is interreduced
%  . the basis is normalized by leading coefficients
%  . the basis is sorted by leading terms increasingly
asserted procedure f5_standardizeOutput(basis: List): List;
  begin;
    basis := f5_interreduceBasis(basis);
    basis := f5_normalizeBasis(basis);
    return sort(basis, 'lp_cmpLPLead)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if <polys> contains in <basis>,
% assuming `basis` is a Groebner basis
asserted procedure f5_checkIdealInclusion1(basis: List, polys: List);
  begin scalar ans, isRegular, p, pNf;
    isRegular := nil;
    ans := t;

    while polys do <<
      p . polys := polys;
      pNf := lp_normalForm(p, basis, isRegular);
      if not lp_iszero!?(pNf) then
        ans := nil
    >>;

    return ans
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if `basis` is a Groebner basis by checking all S-polynomials
asserted procedure f5_isGroebner1(basis: List);
  begin scalar spolys, ans, isRegular, p, pNf;
    % construct S-polynomials
    spolys := f5_constructSpolys(basis);

    isRegular := nil;
    ans := t;

    while spolys do <<
      p . spolys := spolys;
      pNf := lp_normalForm(p, basis, isRegular);
      if not lp_iszero!?(pNf) then
        ans := nil
    >>;

    return ans
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The heart of the package, --
% the function to take the list of ideal generators `basis`
% and to return a Groebner basis of this idel
asserted procedure f5_groebner1(basis: List): List;
  begin scalar basis, spolys, syzygies, isRegular, p, pNf,
                msi, msj, gg;
        integer i, iter;

    % form principal syzygies
    syzygies := f5_constructSyzygies(basis);

    % construct S-polynomials
    spolys   := f5_constructSpolys(basis);

    isRegular := t;
    while spolys do <<
      % select one S-poly
      i . p := f5_selectNext(spolys);

      % inplace ?
      spolys := remove(spolys, i);

      % if syzygy criterion then discard
      if f5_syzygyCriterion(p, syzygies) then
        % pass
      else <<
        % signature-safe reduction
        pNf := lp_normalForm(p, basis, isRegular);

        if lp_isSyzygy(pNf) then
          syzygies := pNf . syzygies
        else if (not lp_isSingularlyTopReducible(pNf, basis)) then <<
          for each gg in basis do <<
            msi . msj := lp_spolyMultSignatures(pNf, gg);
            if not (msi equal msj) then
              spolys := lp_spoly(pNf, gg) . spolys
          >>;
          basis := pNf . basis
        >>
      >>;
      iter := iter #+ 1
    >>;

    basis := f5_standardizeOutput(basis);
    return basis
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trst f5_groebner;

% trst f5_selectNext;
% trst f5_constructSyzygies;
% trst f5_constructSpolys;
% trst f5_constructModule;
% trst f5_syzygyCriterion;
% trst f5_standardizeOutput;
% trst f5_interreduceBasis;

% trst f5_groebner1;
% trst f5_checkIdealInclusion1;
% trst f5_isGroebner1;

% trst f5_groebnerModular1;
% trst f5_randomizedCorrectnessCheck;
% trst f5_crtReconstruction;
% trst f5_modularReduction;
% trst f5_rationalReconstruction;

endmodule;

% f5({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

end;  % of file
