module f5radical;
% The module with the ideal radical implementation based on f5.

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
         prin2t {presentVars, vars};
         poly_initRing({presentVars, 'gradlexgradlex, length(presentVars) - length(vars)})
      >>;
      prin2t {"rr ", saveTorder};
      % w := errorset({'f5_groebner1, mkquote inputBasis}, t, !*backtrace);
      % torder cdr saveTorder;
      % if errorp w then
      %    return nil;
      % outputModule := car w;
      prin2t {"rad0", inputBasis};
      % outputModule := groebner(inputBasis);
      outputModule := inputBasis;
      prin2t {"o ", outputModule};
      for each f in outputModule do <<
         prin2t {f, kernels f,  setdiff(presentVars, vars), setdiff(kernels f, setdiff(presentVars, vars))};
         if null setdiff(kernels f, setdiff(presentVars, vars)) then
            push(f, output)
      >>;
      prin2t {"ororo", output};
      prin2t for each f in output collect lp_eval f;
      output := for each f in output collect
         poly_2a lp_eval f;
      return 'list . output
   end;

endmodule;  % end of module f5radical

end;  % of file