module f5param;
% The module provides parametric assumptions manipulating interface

% The helper module to keep track of assumptions made in parametric F5 computations.
% Provides functions `param_add*` for recording different kinds of assumptions 
% on polynomial coefficients containing parameters.
% Additionaly, provides function `param_clearAssumptions`, which clears all recorded assumptions,
% and function `param_dumpAssumptions`, which returns all recorded assumptions in a list.
% Currently, 4 types of assumptions are recorded in 4 separate lists:
%  . Assumptions implicitly made from the input (param_assumptionsInput);
%     For example, (1/a)*x in the input will sometimes produce an assumption in case `a` is a parameter.
%  . Assumptions on non-vanishing of leading terms in S-polynomial construction (param_assumptionsSpol);
%     That is, if an S-polynomial of f1,f2 is formed, assumptions on non-vanishing of leading coefficients
%     of f1 and f2 are made.
%  . Assumptions on non-vanishing of leading coefficients in polynomial reductions (param_assumptionsRed);
%  . Assumptions obtained when making polynomials monic/primitive depending on f5fractionfree (param_assumptionsNormalize);
%        (active if f5parametricNormalize is ON)
fluid '(param_assumptionsSpol!* param_assumptionsRed!*
        param_assumptionsInput!* param_assumptionsNormalize!*);

% Returns true if the given assumption `poly <> 0` is "interesting", 
% e.g., the following holds
%  1) `poly` is not a literal constant 
asserted procedure param_isAssumptionInteresting(poly: SQ): Boolean;
   if domainp poly then
      nil
   else
      t;

% Tries to "standardize" the given assumption `poly <> 0`
% and returns a standardized variant
asserted procedure param_standardize(poly: SQ): SQ;
   quotfx(numr poly, lc numr poly);

% Add `poly <> 0` to S-polynomial assumptions 
asserted procedure param_addAssumptionSpol(poly: SQ): List;
   if param_isAssumptionInteresting(poly) then
      param_assumptionsSpol!* := lto_insert(param_standardize(poly), param_assumptionsSpol!*);

% Adds `poly <> 0` to Reduction assumptions 
asserted procedure param_addAssumptionRed(poly: SQ): List;
   if param_isAssumptionInteresting(poly) then
      param_assumptionsRed!*  := lto_insert(param_standardize(poly), param_assumptionsRed!*);

% Adds `poly <> 0` to Input assumptions 
asserted procedure param_addAssumptionInput(poly: SQ): List;
   if param_isAssumptionInteresting(poly) then
      param_assumptionsInput!*  := lto_insert(param_standardize(poly), param_assumptionsInput!*);

% Adds `poly <> 0` to Normalize assumptions 
asserted procedure param_addAssumptionNormalize(poly: SQ): List;
   if param_isAssumptionInteresting(poly) then
      param_assumptionsNormalize!*  := lto_insert(param_standardize(poly), param_assumptionsNormalize!*);

% Clears all assumptions
asserted procedure param_clearAssumptions();
<<
   param_assumptionsSpol!* := nil;
   param_assumptionsRed!*  := nil;
   param_assumptionsNormalize!* := nil;
   param_assumptionsInput!* := nil
>>;

% Dumps (returns) all recorded assumptions.
% `param_dumpAssumptions` is a symbolic counterpart of the algebraic function `f5dumpAssumptions`;
% For example,
%  
%  > f5({a*x1*x2^2 + a*x1*x3^2 - b*x1 + a, a*x1^2*x2 + a*x2*x3^2 - b*x2 + a, a*x1^2*x3 + a*x2^2*x3 - b*x3 + a}, {x1,x2,x3}, revgradlex);
%  > f5dumpAssumptions();
%
%  Out:
%                            3               4  8      5  2    3  
%  {input={1},reductions={2*a ,a},spolys={2*a ,a *b,2*a ,a ,2*a ,a},normalize={}}
%
asserted procedure param_dumpAssumptions(u: List);
<< 
   l1 := 'list . for each f in param_assumptionsInput!* collect prepf f;
   l2 := 'list . for each f in param_assumptionsRed!* collect prepf f;
   l3 := 'list . for each f in param_assumptionsSpol!* collect prepf f;
   l4 := 'list . for each f in param_assumptionsNormalize!* collect prepf f;
   {
      'list,
      '(equal Input (reval l1)),
      '(equal Reductions (reval l2)),
      '(equal Spolys (reval l3)),
      '(equal Normalize (reval l4))
   }
>>;

endmodule;  % end of module f5param

end;  % of file