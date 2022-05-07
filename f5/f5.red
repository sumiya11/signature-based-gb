module f5;
% The F5 Algorithm for computing Groebner bases.
%
% The f5 module provides the implementation of the Faugère's F5 algorithm
%     https://www-polsys.lip6.fr/~jcf/Papers/F02a.pdf
%
% The interface contains the operator `f5` with the following signature
%     `f5(polynomials: List, vars: List, order: Id)`
% Where
%   . `polynomials` is a list of expressions of Ideal generators;
%   . `vars` is a list of Reduce kernels with respect to which the basis is computed;
%     kernels not present in `vars` are treated as elements from the ground domain;
%   . `order` is a term order to compute the basis in,
%      possible options are `lex`, `revgradlex`;
%
% For example, one can use f5 in the following way:
%  > load_package f5;
%  > f5({x*y + 1, y*z + 1}, {x, y, z}, lex);

create!-package('(f5 f5core f5lp f5poly f5core f5primes f5mod), nil);

% f5fullreduce - If the resulting basis should be fully interreduced.
%                If this is ON, each generator in the output basis is
%                in the normal form with respect to other generators.
%                Otherwise, only head terms in the basis are reduced.
%                Is OFF by default.
switch f5fullreduce;
off1 'f5fullreduce;

% f5integers - If this is ON, then coefficients of polynomials
%              in the output basis have denominator 1, and the numerator
%              parts have unit content.
%              Otherwise, each polynomial in the output is divided
%              by the leading coefficient.
switch f5integers;
off1 'f5integers;

% Not exported and should not be used directly.
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

endmodule;

end;  % of file
