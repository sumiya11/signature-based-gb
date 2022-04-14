module f5primes;

% The helper module to keep track of prime numbers used in F5 modular computations

% The standard largest!-small!-modulus value is 2^23, which is 8388608
% Thus, lets take initial_prime!* to be
%   largest!-small!-modulus / 2
%
% We will use primes in range [initial_prime!*, 2*initial_prime!*]
% for modular computations as candidates for lucky primes.
% There are 268216 primes between 2^22 and 2^23,
% so it should be enough in most cases
initial_prime!* := largest!-small!-modulus / 2;

% We expect this to hold
if not (largest!-small!-modulus = 2^23) then
  prin2t {"*** Strange largest!-small!-modulus value: expected 2^23, found", largest!-small!-modulus};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if the given `prime` is lucky w.r.t. the input `basis` coefficients,
% It is sufficient to check that generator's terms do not vanish under reduction
asserted procedure primes_isLuckyPrime(basis, prime);
  begin;
    flag := t;
    for each poly in basis do <<
      if flag then <<
        coeffs := poly_getCoeffs(poly);
        for each cf in coeffs do <<
          if modular!-number(cf) = 0 then
            flag := nil
        >>
      >>
    >>;
    return flag
  end;

% Returns the next lucky prime number greater than `prime`
asserted procedure primes_nextLuckyPrime(basis, prime);
  begin integer nextprime;
    nextprime := primes_nextPrime(prime);
    set!-small!-modulus nextprime;
    if not primes_isLuckyPrime(basis, nextprime) then
      nextprime := primes_nextLuckyPrime(basis, nextprime);
    return nextprime
  end;

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
% trst primes_nextLuckyPrime;
% trst primes_isLuckyPrime;

endmodule;

end; % eof
