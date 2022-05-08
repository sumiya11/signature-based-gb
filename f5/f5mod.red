module f5mod;
% Modular computations for Groebner bases.
%
% The module provides rational reconstruction
% and Chinese reminder theorem implementations
% together with main modular GB procedure `mod_groebnerModular1`.

% Using small modular arithmetic backend from smallmod
load!-package 'smallmod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR CORRECTNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verifies that the given `recBasis` is indeed the Groebner basis
% of the original input list `intBasis`.
% By default, randomized correctness check is used
asserted procedure mod_correctnessCheck(pt: Primetracker, intBasis: List,
                                          recBasis: List): Boolean;
  if not mod_heuristicCorrectnessCheck(recBasis) then
    nil
  else if not mod_randomizedCorrectnessCheck(pt, intBasis, recBasis) then
    nil
  else if !*f5certify then
    mod_guaranteedCorrectnessCheck(reconstructedBasis)
  else
    t;

% Heuristic correctness check based on the bitsize of coefficients,
% Not needed currently
asserted procedure mod_heuristicCorrectnessCheck(reconstructedBasis: List);
  t;

% Randomized correctness check, checks that
%   . ideal <intBasis> contains in ideal <recBasis>
%   . recBasis is a Groebner basis
% modulo chosen `reliablePrime`
asserted procedure mod_randomizedCorrectnessCheck(pt: Primetracker, intBasis: List,
                                                    recBasis: List): Boolean;
  begin integer reliablePrime;
        scalar intBasisReduced, recBasisScaled, recBasisReduced;
    reliablePrime := primes_nextReliablePrime(pt, intBasis);
    % reduce coefficients in the integer basis intBasis modulo reliablePrime..
    intBasisReduced := mod_modularReduction(intBasis, reliablePrime);
    % ..and do the same for recBasis
    recBasisScaled := mod_scaleDenominators(recBasis);
    recBasisReduced := mod_modularReduction(recBasisScaled, reliablePrime);
    % all S-polys reduce to zero ?
    return if core_isGroebner1(recBasisReduced) then
      % ideal inclusion ?
      core_checkIdealInclusion1(recBasisReduced, intBasisReduced)
    else
      nil
  end;

% Guaranteed correctness check. Same as mod_randomizedCorrectnessCheck,
% but all computations are carried in the original ground ring,
% Not needed currently
asserted procedure mod_guaranteedCorrectnessCheck(reconstructedBasis: List);
  t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All functions below that take a list of basis generators as an input
% should preserve the original order of generators in the output

% Given a list of generators `inputBasis` with integer coefficients constructs
% a new basis with generators coefficients reduced modulo `prime`
asserted procedure mod_modularReduction(inputBasis: List, prime: Integer): List;
  for each poly in inputBasis collect lp_reduceCoeffs(poly, prime);

% Given the basis `inputBasis` with generators coefficients modulo some prime,
% reconstructs each coefficient into the rational field.
% Rationals are represented as Standard Quotients (pairs) of integers.
% This function is assumed to be called only when the
% coefficients of generators of `inputBasis` are integers
asserted procedure mod_rationalReconstruction(
                                      inputBasis: List,
                                      pt: Primetracker): List;
   begin scalar ans;
          integer prime;
      prime := primes_getAccumModulo(pt);
      while inputBasis do <<
         push(lp_reconstructCoeffs(pop(inputBasis), prime), ans)
      >>;
      return reversip(ans)
   end;

% Given the basis `accumBasis` computed modulo some `modulo`,
% and the basis `computedBasis` computed modulo some `prime`,
% constructs a new basis modulo `modulo × prime`.
% This function is assumed to be called only
% when coefficients of input bases are integers
asserted procedure mod_crtReconstruction(
                                  accumBasis: List,
                                  computedBasis: List,
                                  pt: Primetracker);
  begin scalar ans, polyaccum, polycomp;
        integer modulo, prime;
    modulo := primes_getAccumModulo(pt);
    prime  := primes_getLuckyPrime(pt);
    if modulo = 1 then
      % essentially, accumBasis is garbage here,
      % and we return computedBasis
      ans := computedBasis
    else <<
      while accumBasis do <<
        polyaccum := pop(accumBasis);
        polycomp  := pop(computedBasis);
        push(lp_crtCoeffs(polyaccum, modulo, polycomp, prime), ans)
      >>;
      ans := reversip(ans)
    >>;
    primes_setAccumModulo(primetracker, modulo * prime);
    return ans
  end;

% Given the set of generators `inputBasis` with SQ
% as coefficients, multiplies each generator by the common denominator
% of its coefficients to obtain integer coefficients.
% Assumed to be called only for rational number coefficients in input.
asserted procedure mod_scaleDenominators(inputBasis: List): List;
  for each poly in inputBasis collect lp_scaleDenominators(poly);

% Main procedure for modular F5 Groebner basis computation,
% !! Takes list of module elements as input and returns a list of module elements as output
asserted procedure mod_groebnerModular1(inputBasis: List): List;
  begin scalar integerBasis, reducedBasis, computedBasis,
                reconstructedBasis, accumBasis,
                correctness, primetracker;
        integer prime;
    % contains current lucky prime for reduction and accumulates
    % the product of all previous ones
    primetracker := primes_Primetracker();
    % scale coefficients to integers
    % if not !*f5integers then
    %  integerBasis := mod_scaleDenominators(inputBasis)
    % else
      integerBasis := inputBasis;
    % now all denominators are 1 and coefficients are actually "big" integers
    % in a sense that it's not safe to use machine arithmetic.

    % this basis will store accumulated coefficients
    % from all previous CRT calls
    accumBasis := nil;
    % so, while the correctness check failes
    while not correctness do <<
      % select the next lucky prime
      prime := primes_nextLuckyPrime(primetracker, integerBasis);
      ASSERT(primep prime);
      % reduce the basis w.r.t the prime
      reducedBasis := mod_modularReduction(integerBasis, prime);
      % now all coefficients in reducedBasis should be small integers.
      % Compute the basis in the finite field
      computedBasis := core_groebner1(reducedBasis);
      % CRT (prevprime, prime) --> (prevprime*prime)
      accumBasis := mod_crtReconstruction(accumBasis, computedBasis, primetracker);
      % reconstruct modulo modulo
      reconstructedBasis := mod_rationalReconstruction(accumBasis, primetracker);
      correctness := mod_correctnessCheck(primetracker,
                                            integerBasis,
                                            reconstructedBasis)
    >>;
    return reconstructedBasis
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rational number reconstruction implementation from CLUE
%   https://arxiv.org/abs/2004.11961
% modified a bit to suit the 'Modern Computer Algebra' definitions.
% Returns a rational r // h in canonical form such that
%   r // h ≡ a (mod m)
%
% Here, a and m can be *arbitrary large* integers
asserted procedure mod_reconstruction(a: Integer, m: Integer): SQ;
  begin integer ans, x, com, u1, u2, u3, q, bnd,
                v1, v2, v3, t1, t2, t3, tt, r;
    if a = 0 then
        ans := 0 . 1
    else if a = 1 then
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

% Standard Euclidean algorithm, returns the gcd of a and b
%
% Here, a and b can be *arbitrary large* integers
asserted procedure mod_euclid(a: Integer, b: Integer): Integer;
  if b = 0 then
      a
  else
    mod_euclid(b, a mod b);

% Standard Extended Euclidean algorithm,
% Returns x and y such that a*x + b*y = 1.
% Assumes that a and b are coprime and *arbitrary large* integers
asserted procedure mod_extendedEuclid(a: Integer, b: Integer): DottedPair;
  begin integer x, y, k;
    return if b = 0 then
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
asserted procedure mod_crt(a1: Integer, m1: Integer,
                            a2: Integer, m2: Integer): Integer;
  begin integer x, y, n, m;
    m := m1 * m2;
    x . y := mod_extendedEuclid(m1, m2);
    n := a2 * x * m1 + a1 * y * m2;
    return ((n mod m) + m) mod m
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endmodule;  % end of module f5mod

end;  % of file
