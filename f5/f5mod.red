module f5mod;

% The module provides rational reconstruction routine
% together with generic arithmetic interface

% Using modular arithmetic backend from smallmod
load_package 'smallmod;

% If f5 should use modular arithmetic
switch f5modular;
on1 'f5modular;

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

% trst mod_reconstruction;

endmodule;

end;  % of file
