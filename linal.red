
% Computes the normal form of f w.r.t. polynomials r_i where i in Gprev
%
% Two cases are possible:
%  1. At least one reduction step happened,
%     then we return the t flag together with the normal form itself;
%  2. No reductions were performed, then nil flag and f itself are returned.
%
% During reductions, the topReduce flag controls the reduction type.
% If topReduce is set, only top-reductions happen. Otherwise,
% the polynomial is fully reduced.
asserted procedure core_normalForm(f: Polynomial, Gprev: List,
                                   r: Basistracker,
                                   topReduce: Boolean): DottedPair;
   begin scalar reducers, updated, reducer, reduced, f, g;
      % while polynomial f gets updated by reduction steps,
      % scan the list Gprev in search for possible reducers
      % prin2t {"####################"};
      % prin2t {"Normal form called with Gprev = ", Gprev};
      % prin2t {for each g in Gprev collect
      %         if not poly_iszero!?(lp_eval(core_getPoly(r, g))) then 
      %            poly_leadTerm(lp_eval(core_getPoly(r, g)))
      %         else 
      %            nil
      %      };
      topReduce := t; % for now, always true
      updated := t;
      while updated do <<
         updated := nil;
         reducers := Gprev;
         while reducers and (not poly_iszero!?(f)) do <<
            g := pop(reducers);
            reducer := lp_eval(core_getPoly(r, g));
            reduced . f := poly_tryTopReductionStep(f, reducer);
            % prin2t {"reducing by ", g, "...", "reduced? ", reduced};
            if reduced then
               push(g, reducers);
            updated := reduced or updated
         >>
      >>;
      return updated . f
   end;

load_package sparse;

sparse A(1, 1);

% in the real world, the iteration number is unknown,
% and is roughly bounded by 50 
iters := 2;
for i := 1:iters do <<
    % on each iteration:
    %   0. check if Ax = b is solvable (just check the rank)
    %   1. add several new rows (~99% zero entries)
    %   2. add one new column (~90% zero entries)
    %   3. check if solvable
    x := for j := 1:(car length(A)) do collect j;

    newrows := i;
    A := spextend(A, newrows, 0, 0);
    A(i, 1) := 1
>>;

end;