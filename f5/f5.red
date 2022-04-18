% F5 algorithm implementation module
%
% The module provides Groebner basis computation algorithm `f5` based on
% Fougeres F5 algorithm
module f5;

% The interface consists of the `f5` function for computing Groebner bases.
%
% The signature of `f5` is `f5(polynomials, variables, ordering)` where
%   . `polynomials` is the list of generators of ideal;
%   . `variables` is the list of variables to w.r.t. with,
%     kernels not present in variables are treated as constants.
%   . `ordering` is the monomial ordering to compute the basis in.
%     Possible options are `lex`, `revgradlex`,
%
% The algorithm is randomized. By default, the obtained result
% will be correct with high probability.
create!-package('(f5 f5lp f5poly f5mod f5primes), nil);

load_package 'assert;
on1 'assert;

% If f5 should certify the answer.
% False by default, meaning that the algorithm is randomized
% and may output incorrect answer with a very small probability (~1/2^22)
switch f5certify;
off1 'f5certify;

put('f5, 'psopfn, 'f5_groebner);

struct Polynomial;
struct LabeledPolynomial;
struct Primetracker;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure f5_groebner(u: List): List;
   begin scalar inputBasis, variables, sortMode,
                inputModule, outputModule;
      if null u then
         f5_error();
      inputBasis := reval pop u;
      if not (pop inputBasis eq 'list) then
         f5_error();
      variables := reval pop u;
      if not (pop variables eq 'list) then
         f5_error();
      sortMode := pop u;
      if not null u then
         f5_error();
      poly_initRing(variables, sortMode);
      inputBasis := for each f in inputBasis collect
         poly_f2poly numr simp f;

      inputModule := f5_constructModule(inputBasis);

      if !*f5modular then
         outputModule := f5_groebnerModular1(inputModule)
      else
         outputModule := f5_groebner1(inputModule);

      outputModule := 'list . for each f in outputModule collect
                        poly_poly2a lp_evaluation f;
      return outputModule
   end;

% Void return type fails assert check
asserted procedure f5_error();
   rederr "usage: f5(polynomials: List, variables: List, sortmode: Id)";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR CORRECTNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verify that the given reconstructed basis is indeed the Groebner basis
% of the original input set
% By default, probabilistic and heuristic correctness checks are used
asserted procedure f5_correctnessCheck(pt, intBasis, recBasis);
   if not f5_heuristicCorrectnessCheck(recBasis) then
    nil
  else if not f5_randomizedCorrectnessCheck(pt, intBasis, recBasis) then
    nil
  else if !*f5certify then
    f5_guaranteedCorrectnessCheck(reconstructedBasis)
  else
    t;

asserted procedure f5_heuristicCorrectnessCheck(reconstructedBasis: List);
  t;

asserted procedure f5_randomizedCorrectnessCheck(pt, intBasis, recBasis);
  begin integer reliable_prime;
        scalar intBasisReduced, recBasisScaled, recBasisReduced,
                isgroebner;
    reliable_prime := primes_nextReliablePrime(pt, intBasis);

    intBasisReduced := f5_modularReduction(intBasis, reliable_prime);

    recBasisScaled := f5_scaleDenominators(recBasis);
    recBasisReduced := f5_modularReduction(recBasisScaled, reliable_prime);

    isgroebner := f5_isGroebner(recBasisReduced);

    return if isgroebner then
      f5_checkIdealInclusion(recBasisReduced, intBasisReduced)
    else
      nil
  end;

asserted procedure f5_guaranteedCorrectnessCheck(reconstructedBasis: List);
  t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given a set of generators with integer coefficients construct
% a new basis with coefficients reduced modulo current prime
asserted procedure f5_modularReduction(inputBasis: List, prime): List;
  begin scalar prime, ans, poly;
    while inputBasis do <<
      poly . inputBasis := inputBasis;
      ans := lp_reduceCoeffs(poly, prime) . ans
    >>;
    return reversip(ans)
    end;

% Given a basis computed modulo some prime,
% reconstructs each coefficient into rational field.
% Rationals are represented as Standard Quotients of integers
% (since this function is assumed to called only for number case)
asserted procedure f5_rationalReconstruction(inputBasis: List, primetracker): List;
   begin scalar prime, poly, ans;
      prime := primes_getAccumModulo(primetracker);
      while inputBasis do <<
         poly . inputBasis := inputBasis;
         ans := lp_reconstructCoeffs(poly, prime) . ans
      >>;
      return reversip(ans)
   end;

asserted procedure f5_crtReconstruction(accumBasis, computedBasis, primetracker);
   begin scalar modulo, prime, ans, polyaccum, polycomp;
      modulo := primes_getAccumModulo(primetracker);
      prime  := primes_getLuckyPrime(primetracker);
      if modulo #= 1 then
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

% Given a set of generators with rational numbers as polynomial coefficients,
% scale each generator w.r.t. a common denominator to obtain integer coefficients
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

asserted procedure f5_scaleDenominators(inputBasis: List): List;
  begin scalar ans, poly;
    while inputBasis do <<
      poly . inputBasis := inputBasis;
      ans := lp_scaleDenominators(poly) . ans
    >>;
    return reversip(ans)
  end;

% Main function for modular F5 Groebner basis computation
asserted procedure f5_groebnerModular1(inputBasis: List): List;
  begin scalar integerBasis, reducedBasis, computedBasis,
                reconstructedBasis, accumBasis,
                correctness, prime, primetracker;
        integer iter;

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
        % if iter #> 100 then correctness := t
      >>;

      return reconstructedBasis
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERIC F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for a given set of polynomials inputBasis construct
% a standard basis of unit vectors corresponding to F
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

% for a given set of LabeledPolynomials inputBasis construct
% the set of all principal syzygies
asserted procedure f5_constructSyzygies(inputBasis: List): List;
  begin scalar syzygies, inputSlice, p1, p2;
        integer i, n;
    n := length(inputBasis);
    for i := 1:n do <<
      inputSlice := cdr inputBasis;
      for j := i+1:n do <<
        p1 := car inputBasis;
        p2 := car inputSlice;
        syzygies := lp_principalSyzygy(p1, p2) . syzygies;
        syzygies := cdr syzygies
      >>;
      inputBasis := cdr inputBasis
    >>;
    return syzygies
  end;

% for a given set of LabeledPolynomials inputBasis construct
% the set of all S-polynomials
asserted procedure f5_constructSpolys(inputBasis: List): List;
  begin scalar spolys, inputSlice, p1, p2;
        integer i, n;
    n := length(inputBasis);
    for i := 1:n do <<
      inputSlice := cdr inputBasis;
      for j := i+1:n do <<
        p1 := car inputBasis;
        p2 := car inputSlice;
        if not poly_disjLead!?(lp_evaluation(p1), lp_evaluation(p2)) then
            spolys := lp_spoly(p1, p2) . spolys;
        inputSlice := cdr inputSlice
      >>;
      inputBasis := cdr inputBasis
    >>;
    return spolys
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% checks if computing normal form of `p` is redundant
% with already known `syzygies`
asserted procedure f5_syzygyCriterion(p: LabeledPolynomial, syzygies: List);
   begin scalar flag, syz, sgn;
      if null syzygies then
         flag := nil
      else <<
         syz := lp_signature(car syzygies);
         sgn := lp_signature(p);

         if lp_signatureDivides(syz, sgn) then
            flag := t
         else
            flag := f5_syzygyCriterion(p, cdr syzygies)
      >>;
      return flag
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the polynomial with the lowest signature in spolys and return
% it together with its index in spolys
asserted procedure f5_selectNext(spolys: List): List;
  begin integer i, idx;
        scalar sgn, elem, s, p;
    i := 1;
    s := lp_signature(car spolys); % nonempty
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

asserted procedure f5_normalizeBasis(basis: List): List;
  for each x in basis collect lp_normalize(x);

% given a Groebner basis
% returns the unique Groebner basis of the corresponding ideal
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

asserted procedure f5_checkIdealInclusion(basis: List, polys: List);
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

asserted procedure f5_isGroebner(basis: List);
  begin scalar spolys, ans, isRegular, p, pNf;
    % construct S-polynomials
    spolys   := f5_constructSpolys(basis);

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
      % if iter #> 1 then spolys := nil
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
% trst f5_checkIdealInclusion;
% trst f5_isGroebner;

% trst f5_groebnerModular1;
% trst f5_randomizedCorrectnessCheck;
% trst f5_crtReconstruction;
% trst f5_modularReduction;
% trst f5_rationalReconstruction;

endmodule;

% f5({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

end;  % of file
