
module f5mod;

% Modular computationы for Groebner bases.
% The module also provides rational reconstruction
% and Chinese reminder theorem routines
% together with generic number arithmetic interface

% Using modular arithmetic backend from smallmod
load!-package 'smallmod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ADAPTIVE ARITHMETIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inline procedure mod_iszero!?(a);
  if !*f5modular then
    a = 0
  else
    numr(a) = nil;

inline procedure mod_add(a, b);
    if !*f5modular then
        modular!-plus(a, b)
    else
        addsq(a, b);

inline procedure mod_mul(a, b);
    if !*f5modular then
        modular!-times(a, b)
    else
        multsq(a, b);

inline procedure mod_neg(a);
    if !*f5modular then
        modular!-minus(a)
    else
        negsq(a);

inline procedure mod_div(a, b);
    if !*f5modular then
        modular!-quotient(a, b)
    else
        quotsq(a, b);

inline procedure mod_inv(a);
    if !*f5modular then
        modular!-reciprocal(a)
    else
        denr(a) . numr(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR CORRECTNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verify that the given `recBasis` is indeed the Groebner basis
% of the original input set `intBasis`
% By default, probabilistic correctness check is used
asserted procedure mod_correctnessCheck(pt, intBasis, recBasis);
   if not mod_heuristicCorrectnessCheck(recBasis) then
    nil
  else if not mod_randomizedCorrectnessCheck(pt, intBasis, recBasis) then
    nil
  else if !*f5certify then
    mod_guaranteedCorrectnessCheck(reconstructedBasis)
  else
    t;

% Heuristic correctness check based on the bitsize of coefficients
% NOT TODO
asserted procedure mod_heuristicCorrectnessCheck(reconstructedBasis: List);
  t;

% Randomized correctness check, checks that
%   . <intBasis> in <recBasis> as ideals
%   . recBasis is a Groebner basis
% modulo chosen `reliablePrime`
asserted procedure mod_randomizedCorrectnessCheck(pt, intBasis, recBasis);
  begin integer reliablePrime;
        scalar intBasisReduced, recBasisScaled, recBasisReduced,
                isgroebner;
    reliablePrime := primes_nextReliablePrime(pt, intBasis);

    intBasisReduced := mod_modularReduction(intBasis, reliablePrime);

    recBasisScaled := mod_scaleDenominators(recBasis);
    recBasisReduced := mod_modularReduction(recBasisScaled, reliablePrime);

    isgroebner := core_isGroebner1(recBasisReduced);

    return if isgroebner then
      core_checkIdealInclusion1(recBasisReduced, intBasisReduced)
    else
      nil
  end;

% Guaranteed correctness check. Same as mod_randomizedCorrectnessCheck,
% but all computations are carried in the original ground ring
% NOT TODO
asserted procedure mod_guaranteedCorrectnessCheck(reconstructedBasis: List);
  t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All functions below that take a basis and transform its coefficients
% preserve the order of generators in the original basis

% Given the set of generators `inputBasis` with integer coefficients construct
% a new basis with coefficients reduced modulo `prime`
asserted procedure mod_modularReduction(inputBasis: List, prime: Integer): List;
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
asserted procedure mod_rationalReconstruction(
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
asserted procedure mod_crtReconstruction(
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
asserted procedure mod_scaleDenominatorsInplace(inputBasis: List): List;
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
asserted procedure mod_scaleDenominators(inputBasis: List): List;
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
asserted procedure mod_groebnerModular1(inputBasis: List): List;
  begin scalar integerBasis, reducedBasis, computedBasis,
                reconstructedBasis, accumBasis,
                correctness, primetracker;
        integer iter, prime;

      % contains current lucky prime for reduction and accumulates
      % the product of all previous ones
      primetracker := primes_Primetracker();

      % scale coefficients to integers
      integerBasis := mod_scaleDenominators(inputBasis);
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
         reducedBasis := mod_modularReduction(integerBasis, prime);
         % now all coefficients in reducedBasis should be small integers

         % compute the basis in the finite field
         computedBasis := core_groebner1(reducedBasis);

         % CRT (prevprime, prime) --> (prevprime*prime)
         accumBasis := mod_crtReconstruction(accumBasis, computedBasis, primetracker);

         % reconstruct modulo modulo
         reconstructedBasis := mod_rationalReconstruction(accumBasis, primetracker);

         correctness := mod_correctnessCheck(primetracker,
                                             integerBasis,
                                             reconstructedBasis);
         % correctness := t;
         iter := iter + 1
      >>;

      return reconstructedBasis
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rational number reconstruction implementation from CLUE
%   https://arxiv.org/abs/2004.11961
% modified a bit to suit the 'Modern Computer Algebra' definitions.
% Returns rational r // h in canonical form such that
%   r // h ≡ a (mod m)
%
% Here, a and m can be *arbitrary large* integers
asserted procedure mod_reconstruction(a: Integer, m: Integer);
  begin integer ans, x, com, u1, u2, u3, q, bnd,
                v1, v2, v3, t1, t2, t3, tt, r;
    if a #= 0 then
        ans := 0 . 1
    else if a #= 1 then
        ans := 1 . 1
    else <<
        bnd := isqrt(m / 2);

        u1 := 1;
        u2 := 0;
        u3 := m;

        v1 := 0;
        v2 := 1;
        v3 := a;
        while abs(v3) >= bnd do <<
            q := u3 / v3;
            t1 := u1 - q * v1;
            t2 := u2 - q * v2;
            t3 := u3 - q * v3;

            u1 := v1;
            u2 := v2;
            u3 := v3;

            v1 := t1;
            v2 := t2;
            v3 := t3
        >>;

        tt := abs(v2);
        r := v3 * sgn(v2);

        if tt < 0 then <<
          tt := - tt;
          r  := - r
        >>;

        com := mod_euclid(abs(r), tt);
        ans := (r / com) . (tt / com)
    >>;

    return ans
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard Euclidean algorithm
%
% Here, a and b can be *arbitrary large* integers
asserted procedure mod_euclid(a: Integer, b: Integer);
  if b #= 0 then
      a
  else
    mod_euclid(b, a mod b);

% Standard Extended Euclidean algorithm,
% Returns x and y such that a*x + b*y = 1.
% Assumes that a and b are coprime and *arbitrary large* integers
asserted procedure mod_extendedEuclid(a, b);
  begin integer x, y, k;
    return if b #= 0 then
      1 . 0
    else <<
      x . y := mod_extendedEuclid(b, a mod b);
      k := x - (a / b) * y;
      y . k
    >>
  end;

% Chinese reminder theorem reconstruction implementation
% Returns integer x such that
%   x ≡ a1 (mod m1)
%   x ≡ a2 (mod m2)
asserted procedure mod_crt(a1, m1, a2, m2);
  begin integer x, y, n, m;
    m := m1 * m2;
    x . y := mod_extendedEuclid(m1, m2);
    n := a2 * x * m1 + a1 * y * m2;
    return ((n mod m) + m) mod m
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trst mod_reconstruction;
% trst mod_crt;
% trst mod_extendedEuclead;
trst mod_groebnerModular1;
trst mod_randomizedCorrectnessCheck;

endmodule;

end;  % of file
