module f5primes;
% The module provides lucky prime numbers manipulation

% The helper module to keep track of prime numbers used in F5 modular computations.
% Provides the structure `Primetracker` to store intermediate "lucky" and "reliable" primes.
% Primetracker is represented as a list of integers with a convenience tag
%   {'pt, lucky prime, reliable prime, accum modulo}
% where the first two are guaranteed to be small integers < largest!-small!-modulus.

% We call the prime "lucky" if it is used directly for modular computation.
% We call the ptime "reliable" if it is used to verify the correctness of
% the modular computation.
% This module stores current primes globally,
% directly accessing these from other modules is forbidden
fluid '(primes_initialLuckyPrime!* primes_initialReliablePrime!*);

% The standard largest!-small!-modulus value is 2^23, which is 8388608
% Thus, lets take primes_initialLuckyPrime!* to be
%   largest!-small!-modulus / 2
%
% We will use lucky primes in range
%   [largest!-small!-modulus / 2, largest!-small!-modulus]
% for modular computations as candidates for lucky primes.
% There are 268216 primes between 2^22 and 2^23,
% so it should be enough in most cases
primes_initialLuckyPrime!* := largest!-small!-modulus / 2;

% All reliable primes would be at any moment bigger than lucky ones,
% Namely, the range would be
%   [largest!-small!-modulus * 3/2, largest!-small!-modulus]
primes_initialReliablePrime!* := (largest!-small!-modulus / 3) * 2;

% We expect this to hold.
% If the initial prime is too small,
% modular computations are not very efficient
if largest!-small!-modulus < 2^10 then
  rederr {"***** Strange largest!-small!-modulus value:",
            "expected at least 2^10, found",
            largest!-small!-modulus};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Primetracker
asserted procedure primes_Primetracker(): Primetracker;
  begin integer luckyPrime, reliablePrime, accumModulo;
    luckyPrime    := primes_initialLuckyPrime!*;
    reliablePrime := primes_initialReliablePrime!*;
    accumModulo   := 1;
    return {'pt, luckyPrime, reliablePrime, accumModulo}
  end;

% Return the current lucky prime.
% This should not be called before the first lucky prime was produced
% (see primes_nextLuckyPrime below)
asserted inline procedure primes_getLuckyPrime(p: Primetracker): Integer;
  cadr p;

% Return the current reliable prime.
% This should not be called before the first reliable prime was produced
% (see primes_nextReliablePrime below)
asserted inline procedure primes_getReliablePrime(p: Primetracker): Integer;
  caddr p;

% Return the accumulated modulo (may be arbitrary large)
asserted inline procedure primes_getAccumModulo(p: Primetracker): Integer;
  cadddr p;

asserted inline procedure primes_setLuckyPrime(p: Primetracker, i: Integer);
  cadr p := i;

asserted inline procedure primes_setReliablePrime(p: Primetracker, i: Integer);
  caddr p := i;

asserted inline procedure primes_setAccumModulo(p: Primetracker, i: Integer);
  cadddr p := i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% LUCKY PRIMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if the given `prime` is lucky w.r.t. the input `basis` coefficients,
% It is sufficient to check that generators coefficients
% do not vanish under reduction modulo `prime`
asserted procedure primes_isLuckyPrime(basis: List, prime: Integer): Boolean;
  begin scalar islucky, poly, cf;
    islucky := t;
    for each poly in basis do <<
      if islucky then <<
        for each cf in poly_getCoeffs(lp_eval(poly)) do <<
          if modular!-number(cf) = 0 then
            islucky := nil
        >>
      >>
    >>;
    return islucky
  end;

% Returns (and sets for 'smallmod) the next lucky prime number for pt
asserted procedure primes_nextLuckyPrime(pt: Primetracker, basis: List): Integer;
  begin integer nextprime;
    nextprime := primes_nextPrime(primes_getLuckyPrime(pt));
    set!-small!-modulus nextprime;
    while not primes_isLuckyPrime(basis, nextprime) do <<
      nextprime := primes_nextPrime(nextprime);
      set!-small!-modulus nextprime
    >>;
    primes_setLuckyPrime(pt, nextprime);
    return nextprime
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% RELIABLE PRIMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if the given `prime` is reliable w.r.t. the input `basis` coefficients,
asserted procedure primes_isReliablePrime(basis: List, prime: Integer): Boolean;
  primes_isLuckyPrime(basis, prime);

% Returns (and sets in 'smallmod) the next reliable prime number for pt
asserted procedure primes_nextReliablePrime(pt: Primetracker, basis: List): Integer;
  begin integer nextprime;
    nextprime := primes_nextPrime(primes_getReliablePrime(pt));
    set!-small!-modulus nextprime;
    while not primes_isReliablePrime(basis, nextprime) do <<
      nextprime := primes_nextPrime(nextprime);
      set!-small!-modulus nextprime
    >>;
    primes_setReliablePrime(pt, nextprime);
    return nextprime
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns the next prime number greater than `prime`,
% recursion depth should rarely exceed 15 for `prime` ~ 2^22
asserted procedure primes_nextPrime(prime: Integer): Integer;
  <<
    if evenp prime then
      prime := prime + 1
    else
      prime := prime + 2;
    if primep prime then
      prime
    else
      primes_nextPrime(prime)
  >>;

endmodule;

end; % eof
