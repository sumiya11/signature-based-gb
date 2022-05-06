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

% f5modular - If f5 should use modular algorithms during computation.
%             Is set OFF by default, so all arithmetic operations
%             take place in the original coefficient domain.
switch f5modular;
off1 'f5modular;

% f5certify - If f5 should certify the correctness of result
%             during modular computation. Is set ON by default,
%             meaning that the output is guaranteed to be correct.
%             Otherwise, the algorithm is randomized and may
%             produce incorrect answer with a small probability
switch f5certify;
off1 'f5certify;

% f5fullreduce - If the resulting basis should be fully interreduced.
%                If this is ON, each generator in the output basis is
%                in the normal form with respect to other generators.
%                Otherwise, only head terms in the basis are reduced.
%                Is OFF by default.
switch f5fullreduce;
off1 'f5fullreduce;

% f5integers - If to perform computations assuming ring arithmetic
%             for coefficients. If set ON, polynomials are not divided
%             by the leading coefficient during reductions.
switch f5integers;
off1 'f5integers;

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

endmodule;

end;  % of file
