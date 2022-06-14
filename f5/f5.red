module f5;
% The F5 Algorithm for computing Groebner bases.

% The f5 module provides the implementation of the Faugère's F5 algorithm
%     https://www-polsys.lip6.fr/~jcf/Papers/F02a.pdf
%
% The interface contains the operator `f5` with the following signature
%     `f5(polynomials: List, vars: List, order: Id)`
% Where
%   . `polynomials` is a list of expressions, the Ideal generators,
%   . `vars` is a list of Reduce identifiers with respect to which the basis is computed;
%     identifiers not present in `vars` are treated as elements from the ground domain,
%   . `order` is the identifier of the term order to compute the basis in,
%      possible options are `lex`, `revgradlex`;
%
% For example, one can use f5 to compute the Groebner basis of
% x*y + 1, y*z + 1 in lex term order with x > y > z in the following way:
%  > load_package f5;
%  > f5({x*y + 1, y*z + 1}, {x, y, z}, lex);
%
%  or
%
%  > load_package f5;
%  > torder({x, y, z}, lex);
%  > f5({x*y + 1, y*z + 1});
%
% For different combinations of torder and input arguments f5 works
% much the same way as groebner:
%
% If torder was called before, then
%     f5(system) uses the order set by torder;
%     f5(system, vars, ord) temporarily shadows torder and uses input arguments.
% If torder was not called before, then
%     f5(system) extracts variables as identifiers present in input system,
%                 and uses default order from torder (lex);
%     f5(system, vars, ord) uses input arguments.

% The f5core file is the heart of the package, it contains
% the implementation of the F5 algorithm with the Rewritten Criterion.
% The f5poly implements `Polynomial` interface, that is used by the
% f5lp to implement a `LabeledPolynomial`. f5core operates mainly
% on `LabeledPolynomial` objects.
% The f5mod and f5primes files provide rational reconstruction
% and lucky prime numbers manipulations, which extends existing
% F5 to a modular setting.
% The f5stat allows recording and printing useful statistics
% for each f5 call.
%
% f5mod and f5primes are not used at the moment; 
create!-package('(f5 f5core f5lp f5poly f5primes f5mod f5stat), nil);
% create!-package('(f5 f5radical f5core f5lp f5poly f5primes f5mod f5stat), nil);


fluid '(!*backtrace);

% Needed for compatibility with torder
fluid '(global!-dipvars!*);
fluid '(vdpsortmode!*);

% Currently, there are three switches available, these are described below
% . f5interreduce (default is OFF)
% . f5integers   (default is OFF)
% . f5statistics (default is OFF)

% f5interreduce - If the output basis should be fully interreduced.
%                If this is ON, each generator in the output basis is
%                in the normal form with respect to other generators.
%                Otherwise, only head terms of polynomials in the basis
%                are reduced (and the size of the basis is minimal).
%                Is OFF by default.
%                Generally, f5 with f5interreduce ON is considerably slower.
%
%  For example,
%  > off f5interreduce;
%  > f5({x^2 + x + y, x*y + y, x^3 + x}, {x, y}, lex);
%
%    {x + y,y}
%
%  > on f5interreduce;
%  > f5({x^2 + x + y, x*y + y, x^3 + x}, {x, y}, lex);
%
%    {x,y}
%
switch f5interreduce;
off1 'f5interreduce;

% f5integers - If this is ON, then coefficients of polynomials
%              in the output basis have denominator 1, and the numerator
%              parts have unit gcd within one polynomial.
%              Otherwise, each polynomial in the output is monic.
%              Is OFF by default.
%
%              Currently, f5integers ON assumes there are no parameters
%              in input coefficients.
%
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

switch f5sugar;
on1 'f5sugar;

% f5statistics - If this is ON, collects and prints the following statistics
%                after each f5 call:
%                the number of polynomials reduced,
%                the number of polynomials reduced to zero,
%                the number of calls to normal form,
%                the range of produced critical pairs degrees.
%                All statistics above are differentiated by the signature index.
%                By default, this is OFF, the information
%                is neither collected nor printed.
switch f5statistics;
off1 'f5statistics;

% NOT exported and should NOT be used directly.
% Currently, f5modular is not available as an option mainly for two reasons:
%  1. f5modular is not compatible with f5integers;
%  2. computation with f5modular can be slower for some examples.
%
% f5modular - If f5 should use modular algorithms during computation.
%             Is set OFF by default, so all arithmetic operations
%             take place in the original coefficient domain.
%
%             Currently, f5modular ON assumes there are no parameters
%             in input coefficients.
switch f5modular;
off1 'f5modular;

% NOT exported and should NOT be used directly.
% f5certify - If f5 should certify the correctness of result
%             during modular computation (when f5modular is ON).
%             Is OFF dy default.
switch f5certify;
off1 'f5certify;

% Assertions should be OFF in production.
load!-package 'assert;
off1 'assert;
off1 'assert_procedures;
off1 'assert_inline_procedures;
off1 'assertinstall;
off1 'evalassert;

% For string manipulations in table printing;
% For `sfto_kernelp`;
load!-package 'rltools;

% The only function in the interface
put('f5, 'psopfn, 'f5_groebner);

% The only function in the interface
put('elimination, 'psopfn, 'f5_elimination);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STRUCTS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interface implemented in f5poly.red
procedure f5_isPolynomial(x); eqcar(x, 'p);
struct Polynomial checked by f5_isPolynomial;
struct Terms checked by listp;
struct Term checked by listp;
struct Coeffs checked by listp;
% Coeff can be either an Integer or a Standard Quotient.
struct Coeff;

% interface implemented in f5lp.red
procedure f5_isLabeledPolynomial(x); eqcar(x, 'lp);
procedure f5_isSignature(x); eqcar(x, 'sgn);
struct LabeledPolynomial checked by f5_isLabeledPolynomial;
struct Signature checked by f5_isSignature;

% interface implemented in f5primes.red
procedure f5_isPrimetracker(x); eqcar(x, 'pt);
struct Primetracker checked by f5_isPrimetracker;

% interface implemented in f5core.red
procedure f5_isBasistracker(x); eqcar(x, 'bt);
procedure f5_isCriticalPair(x); (null x) or eqcar(x, 'cp);
procedure f5_isRewriteRule(x); eqcar(x, 'rr);
struct Basistracker checked by f5_isBasistracker;
struct CriticalPair checked by f5_isCriticalPair;
struct RewriteRule checked by f5_isRewriteRule;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The main function that parses input arguments and calls the f5 routine
%
% First, if the format of the input is correct (see the format in the header of this file),
% each expression in the input is converted to a `Polynomial` object.
%   These Polynomials are then passed to the `core_constructModule` function.
% The `core_constructModule` cleans these polynomials a bit, converts each of
% them to a `LabeledPolynomial` object, and outputs the list of LabeledPolynomials.
%   These LabeledPolynomials are then passed to `core_groebner1` (if f5modular is OFF, default)
% or to `mod_groebnerModular1` (if f5modular is ON). Both of those functions
% return a list that contains a Groebner basis, and each generator in the basis
% is represented as a `LabeledPolynomial`.
%   Finally, each item in the Groebner basis list is converted to a Standard Form,
% and the resulting list is returned.
asserted procedure f5_groebner(u: List): List;
   begin scalar inputBasis, inputBasisSf, properIdeal, f, vars, ord, outputModule,
                saveTorder, w;
      if null u or not (listp u) then
         f5_argumentError();
      inputBasis := reval pop u;
      if not (listp inputBasis) or not (pop inputBasis eq 'list) or null inputBasis then
         f5_argumentError();
      properIdeal := t; while properIdeal and not null inputBasis do <<
         f := numr simp pop inputBasis;
         if numberp f and not null f then
            properIdeal := nil
         else if not null f then  % This line is for Gleb
            push(f, inputBasisSf)
      >>;
      if not properIdeal then
         return {'list, 1};
      inputBasis := reversip inputBasisSf;
      if null inputBasis then
         % This is a bit unclear mathematically, but we go with the design decisions of the groebner
         % package
         return {'list, 0};
      saveTorder := if not null u then <<
         % variables and sort mode are specified in f5 call
         vars := reval pop u;
         if not (listp vars) or not (pop vars eq 'list) then
            f5_argumentError();
         for each w in vars do
           if not sfto_kernelp(w) then
              f5_argumentError();
         ord := pop u;
         poly_initRing({vars, ord})
      >> else if not null cdr global!-dipvars!* then <<
         % variables <> {} and sort mode are specified using torder
         poly_initRing(nil)
      >> else <<
         % variables = {} and sort mode are specified using torder. Take variables from inputBasis.
         for each f in inputBasis do
            vars := union(vars, kernels f);
         vars := sort(vars, 'ordp);
         poly_initRing({vars})
      >>;
      w := errorset({'f5_groebner1, mkquote inputBasis}, t, !*backtrace);
      torder cdr saveTorder;
      if errorp w then
         return nil;
      outputModule := car w;
      return 'list . outputModule
   end;

asserted procedure f5_groebner1(inputBasis: List): List;
   begin scalar inputModule, outputModule;
      % convert input expressions to `Polynomial`s
      inputBasis := for each f in inputBasis collect
         poly_f2poly f;
      % construct the basis of the module,
      % elements of inputModule are `LabeledPolynomial`s
      inputModule := core_constructModule(inputBasis);
      % call the main groebner routine
      if !*f5modular then
         outputModule := mod_groebnerModular1(inputModule)
      else
         outputModule := core_groebner1(inputModule);
      % convert `LabeledPolynomial`s back to expressions
      outputModule := for each f in outputModule collect
         poly_2a lp_eval f;
      return outputModule
   end;

% Argument error
asserted procedure f5_argumentError();
   rederr "usage: f5(polynomials: List, vars: List, order: Any). For example,

          > f5({x*y + 1, y*z + 1}, {x, y, z}, lex);

          Or, using torder:

          > torder({x, y, z}, lex);
          > f5({x*y + 1, y*z + 1});
          ";


asserted procedure f5_elimination(u: List): List;
   begin scalar inputBasis, inputBasisSf, properIdeal, f, vars, ord, outputModule,
                saveTorder, w, presentVars, output;
      if null u or not (listp u) then
         f5_argumentError();
      inputGens := pop u;
      inputBasis := reval inputGens;
      if not (listp inputBasis) or not (pop inputBasis eq 'list) or null inputBasis then
         f5_argumentError();
      properIdeal := t; while properIdeal and not null inputBasis do <<
         f := numr simp pop inputBasis;
         if numberp f and not null f then
            properIdeal := nil
         else if not null f then  % This line is for Gleb
            push(f, inputBasisSf)
      >>;
      if not properIdeal then
         return {'list, 1};
      inputBasis := reversip inputBasisSf;
      if null inputBasis then
         % This is a bit unclear mathematically, but we go with the design decisions of the groebner
         % package
         return {'list, 0};
      saveTorder := <<
         % variables and sort mode are specified in f5 call
         vars := reval pop u;
         if not (listp vars) or not (pop vars eq 'list) then
            f5_argumentError();
         for each w in vars do
           if not sfto_kernelp(w) then
               f5_argumentError();
         for each f in inputBasis do
            presentVars := union(presentVars, kernels f);
         % make sure `vars` are last
         presentVars := nconc(setdiff(presentVars, vars), vars);
         poly_initRing({presentVars, 'gradlexgradlex, length(presentVars) - length(vars)})
      >>;
      % w := errorset({'f5_groebner1, mkquote inputBasis}, t, !*backtrace);
      % torder cdr saveTorder;
      % if errorp w then
      %    return nil;
      % outputModule := car w;
      outputModule := f5_groebner({inputGens});
      % outputModule := inputBasis;
      for each f in outputModule do <<
         if null setdiff(kernels (numr simp f), setdiff(presentVars, vars)) then
            push(f, output)
      >>;
      return 'list . output
   end;

asserted procedure f5_intersection(u: List): List;
   begin scalar inputBasis, inputBasisSf, properIdeal, f, vars, ord, outputModule,
                saveTorder, w, presentVars, output;
      if null u or not (listp u) then
         f5_argumentError();
      inputGens := pop u;
      inputBasis := reval inputGens;
      if not (listp inputBasis) or not (pop inputBasis eq 'list) or null inputBasis then
         f5_argumentError();
      properIdeal := t; while properIdeal and not null inputBasis do <<
         f := numr simp pop inputBasis;
         if numberp f and not null f then
            properIdeal := nil
         else if not null f then  % This line is for Gleb
            push(f, inputBasisSf)
      >>;
      if not properIdeal then
         return {'list, 1};
      inputBasis := reversip inputBasisSf;
      if null inputBasis then
         % This is a bit unclear mathematically, but we go with the design decisions of the groebner
         % package
         return {'list, 0};
      saveTorder := <<
         % variables and sort mode are specified in f5 call
         vars := reval pop u;
         if not (listp vars) or not (pop vars eq 'list) then
            f5_argumentError();
         for each w in vars do
           if not sfto_kernelp(w) then
               f5_argumentError();
         for each f in inputBasis do
            presentVars := union(presentVars, kernels f);
         % make sure `vars` are last
         presentVars := nconc(setdiff(presentVars, vars), vars)
      >>;
      % w := errorset({'f5_groebner1, mkquote inputBasis}, t, !*backtrace);
      % torder cdr saveTorder;
      % if errorp w then
      %    return nil;
      % outputModule := car w;
      outputModule := f5_groebner({inputGens});
      % outputModule := inputBasis;
      for each f in outputModule do <<
         if null setdiff(kernels (numr simp f), setdiff(presentVars, vars)) then
            push(f, output)
      >>;
      return 'list . output
   end;

endmodule;  % end of module f5

end;  % of file
