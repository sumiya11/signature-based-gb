module f5;
% The F5 Algorithm for computing Groebner bases.
%
% The f5 module provides the implementation of the Faugère's F5 algorithm
%     https://www-polsys.lip6.fr/~jcf/Papers/F02a.pdf
%
% The interface contains the operator `f5` with the following signature
%     `f5(polynomials: List, vars: List, order: Id)`
% Where
%   . `polynomials` is a list of expressions, the Ideal generators,
%   . `vars` is a list of Reduce kernels with respect to which the basis is computed;
%     kernels not present in `vars` are treated as elements from the ground domain,
%   . `order` is the identifier of the term order to compute the basis in,
%      possible options are `lex`, `revgradlex`;
%
% For example, one can use f5 to compute the Groebner basis in the following way:
%  > load_package f5;
%  > f5({x*y + 1, y*z + 1}, {x, y, z}, lex);

% The f5core file is the heart of the package, it contains
% the implementation of the F5 algorithm with the Rewritten Criterion.
% The f5poly implements `Polynomial` interface, that is used by the
% f5lp to implement a `LabeledPolynomial`. f5core operates mainly
% on `LabeledPolynomial` objects.
% The f5mod and f5primes files provide rational reconstruction
% and lucky prime numbers manipulations, which extends existing
% F5 to a modular setting.
% create!-package('(f5 f5core f5lp f5poly f5primes f5mod), nil);

% Currently, there are two switches available:
% . f5fullreduce (default is OFF)
% . f5integers   (default is OFF)
% They are described below

% f5fullreduce - If the output basis should be fully interreduced.
%                If this is ON, each generator in the output basis is
%                in the normal form with respect to other generators.
%                Otherwise, only head terms of polynomials in the basis
%                are reduced (and size of the basis is minimal).
%                Generally, f5 with f5fullreduce ON is considerably slower.
%                Is OFF by default.
switch f5fullreduce;
off1 'f5fullreduce;

% f5integers - If this is ON, then coefficients of polynomials
%              in the output basis have denominator 1, and the numerator
%              parts have unit content within one polynomial.
%              Otherwise, each polynomial in the output is monic.
%              Is OFF by default.
% For example,
%  > off f5integers;
%  > f5({5x + y, 2x + 1}, {x, y}, lex);
%          1       5
%    {x + ---,y - ---}
%          2       2
%
%  > on f5integers;
%  > f5({5x + y, 2x + 1}, {x, y}, lex);
%
%    {2*x + 1,2*y - 5}
%
switch f5integers;
off1 'f5integers;

% Not exported and should not be used directly.
% Currently, f5modular is not available as an option mainly for two reasons:
%  1. f5modular is not compatible with f5integers;
%  2. computation with f5modular can be slower for some examples.
%
% f5modular - If f5 should use modular algorithms during computation.
%             Is set OFF by default, so all arithmetic operations
%             take place in the original coefficient domain.
switch f5modular;
off1 'f5modular;

% Not exported and should not be used directly.
% f5certify - If f5 should certify the correctness of result
%             during modular computation (when f5modular is ON).
switch f5certify;
off1 'f5certify;

load!-package 'assert;
off1 'assert;
off1 'assert_procedures;

% The only function in the interface
put('f5, 'psopfn, 'f5_groebner);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STRUCTS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interface implemented in f5poly.red
inline procedure f5_isPolynomial(x); eqcar(x, 'p);
struct Polynomial checked by 'f5_isPolynomial;
struct Terms checked by 'listp;
struct Term checked by 'listp;
struct Coeffs checked by 'listp;
% Coeff can be either an Integer or a Standard Quotient.
% If there are parameters in the input coefficients,
% f5integers should be set OFF.
struct Coeff;

% interface implemented in f5lp.red
inline procedure f5_isLabeledPolynomial(x); eqcar(x, 'lp);
inline procedure f5_isSignature(x); eqcar(x, 'sgn);
struct LabeledPolynomial checked by 'f5_isLabeledPolynomial;
struct Signature checked by 'f5_isSignature;

% interface implemented in f5primes.red
inline procedure f5_isPrimetracker(x); eqcar(x, 'pt);
struct Primetracker checked by 'f5_isPrimetracker;

% interface implemented in f5core.red
inline procedure f5_isBasistracker(x); eqcar(x, 'bt);
inline procedure f5_isCriticalPair(x); eqcar(x, 'cp);
inline procedure f5_isRewriteRule(x); eqcar(x, 'rr);
struct Basistracker checked by 'f5_isBasistracker;
struct CriticalPair checked by 'f5_isCriticalPair;
struct RewriteRule checked by 'f5_isRewriteRule;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The main function to parse input arguments and call the f5 routine
asserted procedure f5_groebner(u: List): List;
   begin scalar inputBasis, vars, ord,
                inputModule, outputModule;
      if null u or not (listp u) or not (length u = 3) then
         f5_argumentError();
      inputBasis := reval pop u;
      if not (listp inputBasis)then
         f5_argumentError();
      vars := reval pop u;
      if not (listp vars) then
         f5_argumentError();
      ord := pop u;
      if not (idp ord) then
        f5_argumentError();
      % initialize base polynomial ring
      poly_initRing(vars, ord);
      % convert input expressions to polynomials of our type
      inputBasis := for each f in inputBasis collect
         poly_f2poly numr simp f;
      % construct the basis of the module,
      % elements of inputModule are `LabeledPolynomial`s
      inputModule := core_constructModule(inputBasis);
      % call the main groebner routine
      if !*f5modular then
         outputModule := mod_groebnerModular1(inputModule)
      else
         outputModule := core_groebner1(inputModule);
      % convert `LabeledPolynomial`s back to expressions
      outputModule := 'list . for each f in outputModule collect
                        poly_poly2a lp_eval f;
      return outputModule
   end;

% Argument error
asserted procedure f5_argumentError();
   rederr "usage: f5(polynomials: List, vars: List, order: Id). For example,

          f5({x*y + 1, y*z + 1}, {x, y, z}, lex);";

endmodule;  % end of module f5

end;  % of file
