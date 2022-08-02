module f5mod;
% Modular computations for Groebner bases.
%
% The module provides rational reconstruction
% and Chinese reminder theorem implementations
% together with main modular GB procedure `mod_groebnerModular1`.

% added
% The f5mod and f5primes files provide rational reconstruction
% and lucky prime numbers manipulations, which extends existing
% F5 to a modular setting.

% The files f5mod.red and f5primes.red extend the existing f5
% with the modular Groebner basis computation.
% ("extend" in a sense that they can be removed and main f5 still works)
%
% The structure of the Groebner bases modular computation is the following.
%
%   A list of ideal generators is passed as an input to `mod_groebnerModular1`.
%   A "lucky" prime number is then selected. The "lucky" properties of this number
% and the interface for manipulating primes is implemented in f5primes.
% Each coefficients of input generators is reduced to the finite field modulo
% selected prime number.
%   Then, a Groebner basis is computed over the finite field from the previously
% reduced generators. Coefficients in this basis are reconstructed to integers
% using the Chinese Reminder Theorem. Then, integer coefficients are reconstructed
% to rationals using the Rational Reconstruction.
%   From the previous step, we have obtained a Groebner basis over rationals.
% We check that the obtained basis is indeed a Groebner basis of the input ideal,
% and, if it is, return it. Otherwise, we select the next lucky prime,
% and repeat all steps.

% Using small modular arithmetic backend from smallmod
load!-package 'smallmod;

% Not exported, not documented;
% Currently, f5modular is not available as an option mainly for two reasons:
%  1. f5modular is not compatible with f5fractionfree;
%  2. computation with f5modular can be slower for some examples.
%
% f5modular - If f5 should use modular algorithms during computation.
%             Is set OFF by default, so all arithmetic operations
%             take place in the original coefficient domain.
%
%             Currently, f5modular ON assumes there are no parameters
%             in input coefficients.
% switch f5modular;
% off1 'f5modular;

% Not exported, not documented;
% f5certify - If f5 should certify the correctness of result
%             during modular computation (when f5modular is ON).
%             Is OFF dy default.
% switch f5certify;
% off1 'f5certify;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR CORRECTNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verifies that the given `recBasis` is indeed a Groebner basis
% of the original input list `intBasis`.
% By default, randomized correctness check is used
asserted procedure mod_correctnessCheck(pt: Primetracker, intBasis: List,
                                          recBasis: List): Boolean;
   % if not mod_heuristicCorrectnessCheck(recBasis) then
   %   nil
   if not mod_randomizedCorrectnessCheck(pt, intBasis, recBasis) then
      nil
   else
      t;
   % if !*f5certify then
   %   mod_guaranteedCorrectnessCheck(recBasis)
   % else
   %   t;

% Heuristic correctness check based on the bitsize of coefficients,
% Not needed currently
% asserted procedure mod_heuristicCorrectnessCheck(reconstructedBasis: List);
%  t;

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
    on1 'f5modular_internal;
    return if core_isGroebner1(recBasisReduced) then
      % ideal inclusion ?
      core_checkIdealInclusion1(recBasisReduced, intBasisReduced)
    else
      nil
  end;

% Guaranteed correctness check. Same as mod_randomizedCorrectnessCheck,
% but all computations are carried in the original ground ring,
% Not needed currently
% asserted procedure mod_guaranteedCorrectnessCheck(reconstructedBasis: List);
%   t;

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
% constructs a new basis modulo modulo * prime.
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
    primes_setAccumModulo(pt, modulo * prime);
    return ans
  end;

% asserted procedure mod_standardizeInput(inputBasis: List): List;
%    begin scalar f, newcoeffs;
%       if not !*f5fractionfree then
%          for each f in inputBasis do <<
%             newcoeffs := for each cf in lp_setCoeffs
%          >>;
%    end;

% Given the set of generators `inputBasis` with SQ
% as coefficients, multiplies each generator by the common denominator
% of its coefficients to obtain integer coefficients.
% Assumed to be called only for rational number coefficients in input.
asserted procedure mod_scaleDenominators(inputBasis: List): List;
  for each poly in inputBasis collect lp_scaleDenominators(poly);

% Main procedure for modular F5 Groebner basis computation
asserted procedure mod_groebnerModular1(inputBasis: List): List;
  begin scalar integerBasis, reducedBasis, computedBasis,
                reconstructedBasis, accumBasis,
                correctness, primetracker;
        integer prime;
    % contains current lucky prime for reduction and accumulates
    % the product of all previous ones
    primetracker := primes_Primetracker();
    % scale coefficients to integers
    if not !*f5fractionfree then
      integerBasis := mod_scaleDenominators(inputBasis)
    else
      integerBasis := inputBasis;
    % now all denominators are 1 and coefficients are actually "big" integers
    % in a sense that it's not safe to use machine arithmetic.

    on1 'f5modular_internal;

    prin2t {"Int basis", integerBasis};

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
      off1 'f5modular_internal;
      correctness := mod_correctnessCheck(primetracker,
                                            integerBasis,
                                            reconstructedBasis);
      on1 'f5modular_internal
      % correctness := t
    >>;
    off1 'f5modular_internal;
    if !*f5fractionfree then <<
      reconstructedBasis := mod_scaleDenominators(reconstructedBasis);
      reconstructedBasis := for each p in reconstructedBasis collect lp_normalize(p)
   >>;
    return reconstructedBasis
  end;

trst mod_crtReconstruction;
trst mod_rationalReconstruction;
trst mod_correctnessCheck;
trst mod_groebnerModular1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rational number reconstruction implementation from CLUE
%   https://arxiv.org/abs/2004.11961
% modified a bit to suit the 'Modern Computer Algebra, Edition 3' definitions
%   https://doi.org/10.1017/CBO9781139856065
% Returns a rational r // h in canonical form such that
%   r // h = a (mod m)
%
% Here, a and m can be *arbitrary large* integers
asserted procedure mod_reconstruction(a: Integer, m: Integer): SQ;
  begin integer ans, com, u1, u2, u3, q, bnd,
                v1, v2, v3, t1, t2, t3, tt, r;
    if a = 0 then
        ans := 0 ./ 1
    else if a = 1 then
        ans := 1 ./ 1
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
        ans := (r / com) ./ (tt / com)
    >>;
    return ans
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard Euclidean algorithm, returns the gcd of a and b
%
% Here, a and b can be *arbitrary large* integers
asserted procedure mod_euclid(a: Integer, b: Integer): Integer;
  if b = 0 or (null b) then
      a
  else
    mod_euclid(b, a mod b);

% Standard Extended Euclidean algorithm,
% Returns x and y such that a*x + b*y = 1.
% Assumes that a and b are coprime.
%
% a and b can be *arbitrary large* integers (and they are usually large)
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

% Chinese reminder theorem reconstruction implementation.
% Returns integer x such that
%   x = a1 (mod m1)
%   x = a2 (mod m2)
%
% Here, a1, a2, m1, m2 can be *arbitrary large* integers
asserted procedure mod_crt(a1: Integer, m1: Integer,
                            a2: Integer, m2: Integer): Integer;
  begin integer x, y, n, m;
    m := m1 * m2;
    x . y := mod_extendedEuclid(m1, m2);
    n := a2 * x * m1 + a1 * y * m2;
    return ((n mod m) + m) mod m
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reduces coefficients of `f` by the given `prime`
asserted procedure poly_reduceCoeffs(f: Polynomial, prime: Integer): Polynomial;
   begin scalar fcoeffs, newcoeffs, c;
      % note that prime is not used here, since `modular!-number` works globally,
      % to omit warning
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

endmodule;  % end of module f5mod

end;  % of file
