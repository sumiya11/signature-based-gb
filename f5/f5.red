
% Main module of the f5 package
module f5;

% The f5 module provides Groebner basis computation algorithm `f5` implementation
% based on the Fougere's F5 algorithm

% The interface consists of the `f5` function for computing Groebner bases.
%
% - The signature of `f5` is `f5(polynomials, variables, ordering)` where
%   . `polynomials` is the list of ideal generators;
%   . `variables` is the list of kernels to compute basis w.r.t.,
%     kernels not present in `variables` are treated as constants;
%   . `ordering` is the monomial ordering to compute the basis in,
%      possible options are `lex`, `revgradlex`;
%
% - The algorithm uses modular computations and, thus, is randomized.
%   Obtained result will be correct with high probability.
%   If *guaranteed correctness* is needed, use `on f5certify`;
%
% - If using modular computation is not desirable, set `off f5modular`.

create!-package('(f5 f5lp f5poly f5core f5primes f5mod), nil);


% If f5 should certify the correctness of result.
% False by default, meaning that the algorithm is randomized
% and may output incorrect answer with a small probability (~1/2^22)
switch f5certify;
off1 'f5certify;

% If f5 should use modular arithmetic.
% True by default
switch f5modular;
off1 'f5modular;

% If the output basis should be interreduced.
% If true, the basis is unique.
switch f5fullreduce;
on1 'f5fullreduce;

load!-package 'assert;
off1 'assert;

% The only function in the interface
put('f5, 'psopfn, 'f5_groebner);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interface implemented in f5poly.red
struct Polynomial;

% interface implemented in f5lp.red
struct LabeledPolynomial;

% interface implemented in f5primes.red
struct Primetracker;

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

      % initialize ground polynomial ring
      poly_initRing(variables, sortMode);

      inputBasis := for each f in inputBasis collect
         poly_f2poly numr simp f;
      % construct module basis from input polynomials
      inputModule := core_constructModule(inputBasis);

      if !*f5modular then
         outputModule := mod_groebnerModular1(inputModule)
      else
         outputModule := core_groebner1(inputModule);

      outputModule := 'list . for each f in outputModule collect
                        poly_poly2a lp_eval f;
      return outputModule
   end;

% Argument error
asserted procedure f5_argumentError();
   rederr "usage: f5(polynomials: List, variables: List, sortmode: Id)";

endmodule;

% f5({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

end;  % of file
