module f5primes;

% The helper module to keep track of prime numbers used in F5 modular computations.
% Provides the structure `Primetracker` to store intermediate lucky and reliable primes.
% Primetracker is represented internally as a list of integers
%   {lucky prime, reliable prime, accum modulo}
% where the first two guaranteed to be small integers < largest!-small!-modulus.

% We call the prime "lucky" if it is used directly for modular computation.
% We call the ptime "reliable" if it is used to verify the correctness of
% the modular computation.

fluid '(initial_lucky_prime!* initial_reliable_prime!*);

% The standard largest!-small!-modulus value is 2^23, which is 8388608
% Thus, lets take initial_lucky_prime!* to be
%   largest!-small!-modulus / 2
%
% We will use lucky primes in range
%   [largest!-small!-modulus / 2, largest!-small!-modulus / 2]
% for modular computations as candidates for lucky primes.
% There are 268216 primes between 2^22 and 2^23,
% so it should be enough in most cases
initial_lucky_prime!* := largest!-small!-modulus / 2;

% All reliable primes would be at any moment bigger than lucky ones,
% Namely, the range would be
%   [largest!-small!-modulus * 3/2, largest!-small!-modulus]
initial_reliable_prime!* := (largest!-small!-modulus / 3) * 2;

% We expect this to hold
% TODO: make this more verbose
if not (largest!-small!-modulus = 2^23) then
  prin2t {"*** Strange largest!-small!-modulus value: expected 2^23, found", largest!-small!-modulus};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize primetracker
asserted procedure primes_Primetracker();
  begin integer lucky_prime, reliable_prime;
        scalar accum_modulo;
    lucky_prime := initial_lucky_prime!*;
    reliable_prime := initial_reliable_prime!*;
    accum_modulo := 1;
    return {lucky_prime, reliable_prime, accum_modulo}
  end;

asserted inline procedure primes_getLuckyPrime(p: Primetracker): Integer;
  car p;

asserted inline procedure primes_getReliablePrime(p: Primetracker): Integer;
  cadr p;

asserted inline procedure primes_getAccumModulo(p: Primetracker): Integer;
  caddr p;

asserted inline procedure primes_setLuckyPrime(p: Primetracker, i);
  car p := i;

asserted inline procedure primes_setReliablePrime(p: Primetracker, i);
  cadr p := i;

asserted inline procedure primes_setAccumModulo(p: Primetracker, i);
  caddr p := i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Lucky primes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if the given `prime` is lucky w.r.t. the input `basis` coefficients,
% It is sufficient to check that generator's terms do not vanish under reduction
asserted procedure primes_isLuckyPrime(basis, prime);
  begin scalar flag, poly, coeffs, cf;
    flag := t;
    for each poly in basis do <<
      if flag then <<
        coeffs := poly_getCoeffs(lp_getPoly(poly));
        for each cf in coeffs do <<
          if modular!-number(cf) = 0 then
            flag := nil
        >>
      >>
    >>;
    return flag
  end;

% Returns the next lucky prime number for primetracker
asserted procedure primes_nextLuckyPrime(primetracker, basis);
  begin integer nextprime;
    nextprime := primes_nextPrime(primes_getLuckyPrime(primetracker));
    set!-small!-modulus nextprime;
    while not primes_isLuckyPrime(basis, nextprime) do <<
      nextprime := primes_nextPrime(nextprime);
      set!-small!-modulus nextprime
    >>;
    primes_setLuckyPrime(primetracker, nextprime);
    return nextprime
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Reliable primes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if the given `prime` is reliable w.r.t. the input `basis` coefficients,
asserted procedure primes_isReliablePrime(basis, prime: Integer);
  primes_isLuckyPrime(basis, prime);

% Returns the next reliable prime number for primetracker
asserted procedure primes_nextReliablePrime(primetracker, basis);
  begin integer nextprime;
    nextprime := primes_nextPrime(primes_getReliablePrime(primetracker));
    set!-small!-modulus nextprime;
    while not primes_isReliablePrime(basis, nextprime) do <<
      nextprime := primes_nextPrime(nextprime);
      set!-small!-modulus nextprime
    >>;
    primes_setReliablePrime(primetracker, nextprime);
    return nextprime
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns the next prime number greater than `prime`
asserted procedure primes_nextPrime(prime);
  begin;
    % recursion depth should rarely exceed 15 for prime ~ 2^22
    ASSERT(prime > 1);
    if evenp prime then
      prime := prime + 1
    else
      prime := prime + 2;
    return if primep prime then
      prime
    else
      primes_nextPrime(prime);
  end;

% trst primes_nextPrime;
% trst primes_nextReliablePrime;
% trst primes_nextLuckyPrime;
% trst primes_isLuckyPrime;

endmodule;

end; % eof
