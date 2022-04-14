
module f5mod;

% The module provides rational reconstruction routine
% together with generic arithmetic interface

% Using modular arithmetic backend from smallmod
load_package 'smallmod;

% If f5 should use modular arithmetic.
% True by default
switch f5modular;
off1 'f5modular;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ADAPTIVE ARITHMETIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inline procedure mod_iszero!?(a);
    if !*f5modular then
        a = 0
    else
        numr(a) equal nil;

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
%%%%%%%% RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rational number reconstruction implementation modified a bit
% to suit the 'Modern Computer Algebra' definitions
% Returns a rational r // h of Rational field in a canonical form such that
%   r // h ≡ a (mod m)
asserted procedure mod_reconstruction(a, m);
  begin integer ans, x;
    if a equal 0 then
        ans := 0 . 1
    else if a equal 1 then
        ans := 1 . 1
    else <<
        bnd := sqrt(float(m) / 2);

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

        ans := r . tt
    >>;

    return ans
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns x and y such that a*x + b*y = 1.
% Assumes that a and b are coprime
asserted procedure mod_extendedPrimeGcd(a, b, x, y);
  % The function is adapted from the `reciprocal!-by!-gcd` function at
  %   reduce-algebra/packages/rtools/smallmod.red
  if b > a then
    mod_extendedPrimeGcd(b, a, x, y)
  else if b = 1 then
    % TODO: is this actually correct?
    x . y
  else begin scalar w;
    w := a #/ b; % Truncated integer division;
    return mod_extendedPrimeGcd(b, a #- b #* w, y, x #- y #* w)
  end;

% Chinese reminder theorem reconstruction implementation
% Returns integer x such that
%   x ≡ a1 (mod m1)
%   x ≡ a2 (mod m2)
asserted procedure mod_crt(a1, m1, a2, m2);
  begin scalar x, n1, n2, m;
    m := m1 #* m2;
    n1 . n2 := mod_extendedPrimeGcd(m1, m2, 0, 1);
    x := n1 #* m1 #* a2 #+ n2 #* m2 #* a1;
    return x #/ m
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trst mod_reconstruction;
trst mod_crt;
trst mod_extendedPrimeGcd;

endmodule;

end;  % of file
