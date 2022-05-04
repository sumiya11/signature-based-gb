module f5;
% The F5 Algorithm for computing Groebner bases.
%
% The f5 module provides the implementation of the Faugère's F5 algorithm,
% which is in detail documented in the work by Christian Eder and John Perry
%   https://arxiv.org/abs/0906.2967
%
% The interface contains the `f5` operator of signature
%     `f5(polynomials: List, vars: List, order: Id)`
% Where
%   . `polynomials` is a list of expressions of Ideal generators;
%   . `vars` is a list of Reduce kernels to compute the basis w.r.t. with;
%     kernels not present in `vars` are treated as elements from the ground domain;
%   . `order` is a term order to compute the basis in,
%      possible options are `lex`, `revgradlex`;
%
% For example, one can use f5 in the following way:
%  > load_package f5;
%  > f5({x*y + 1, y*z + 1}, {x, y, z}, lex);

create!-package('(f5 f5core f5lp f5poly f5core f5primes f5mod), nil);

% f5modular - If f5 should use modular algorithms during computation.
%             Is set OFF by default, so all arithmetic operations
%             take place in the original coefficient domain.
switch f5modular;
off1 'f5modular;

% f5certify - If f5 should certify the correctness of result
%             during modular computation. Is set ON by default,
%             meaning that the output is guaranteed to be correct.
%             Otherwise, the algorithm is randomized and may
%             produce incorrect answer with a small probability (~1/2^22)
%                                                          TODO?: Proof
switch f5certify;
off1 'f5certify;

% f5fullreduce - If the resulting basis should be fully interreduced.
%                If this is ON, generators in the output basis are
%                both head and tail reduced and the output is unique.
%                Otherwise, the basis is head reduced, but not tail reduced.
%                Is OFF by default.
switch f5fullreduce;
off1 'f5fullreduce;

% TODO?: do this
% f5integer - If to perform computations assuming ring arithmetic
%             for coefficients. If set ON, polynomials are not divided
%             by the leading coefficient during reductions.
% switch f5integer;
% on1 'f5integer;

load!-package 'assert;
on1 'assert;

% Not used for now
% load!-package 'cali;

% The only function in the interface
put('f5, 'psopfn, 'f5_groebner);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STRUCTS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interface implemented in f5poly.red
% asserted inline procedure f5_isPolynomial
struct Polynomial; % checked by function(lambda x; eqcar(x, 'p));
struct Terms checked by 'listp;
struct Term checked by 'listp;
struct Coeffs checked by 'listp;
struct Coeff;

% interface implemented in f5lp.red
struct LabeledPolynomial; % checked by function(lambda x; eqcar(x, 'lp));
struct Signature; % checked by function(lambda x; eqcar(x, 'sgn));

% interface implemented in f5primes.red
struct Primetracker; % checked by function(lambda x; eqcar(x, 'pt));

% interface implemented in f5core.red
struct Basistracker; % checked by function(lambda x; eqcar(x, 'bt));
struct CriticalPair; % checked by function(lambda x; eqcar(x, 'cp));
struct RewriteRule; % checked by function(lambda x; eqcar(x, 'rr));

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
      % initialize base polynomial ring
      poly_initRing(variables, sortMode);
      % convert input expressions to polynomials of our type
      inputBasis := for each f in inputBasis collect
         poly_f2poly numr simp f;
      % construct the basis of the module as LabeledPolynomials
      inputModule := core_constructModule(inputBasis);
      % call the main routine
      if !*f5modular then
         outputModule := mod_groebnerModular1(inputModule)
      else
         outputModule := core_groebner1(inputModule);
      % convert polynomials back to SFs
      outputModule := 'list . for each f in outputModule collect
                        poly_poly2a lp_eval f;
      return outputModule
   end;

% Argument error
asserted procedure f5_argumentError();
   rederr "usage: f5(polynomials: List, variables: List, sortmode: Id)";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ADAPTIVE ARITHMETIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These arithmetic functions are defined here to
% because they are inline and are needed in polynomial arithmetic.
% Suggestion: move it to f5poly.red
inline procedure mod_iszero!?(a);
  if !*f5modular then
    a #= 0
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

endmodule;

end;  % of file
