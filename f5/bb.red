module bb;

off1 'allfac;

load!-package 'dp;

put('buchberger, 'psopfn, 'bb_buchberger);

asserted procedure bb_buchberger(u: List): List;
   begin scalar inputBasis, variables, sortMode, outputBasis;
      if null u then
         bb_error();
      inputBasis := reval pop u;
      if not (pop inputBasis eq 'list) then
         bb_error();
      variables := reval pop u;
      if not (pop variables eq 'list) then
         bb_error();
      sortMode := pop u;
      if not null u then
         bb_error();
      dip_init(variables, sortMode);
      inputBasis := for each f in inputBasis collect
         dip_f2dip numr simp f;
      outputBasis := bb_buchberger1(inputBasis);
      outputBasis := 'list . for each f in outputBasis collect
         dip_2a f;
      return outputBasis
   end;

asserted procedure bb_error(): Void;
   rederr "usage: buchberger(polynomials: List, variables: List, sortmode: Id)";

asserted procedure bb_buchberger1(inputBasis: List): List;
   begin
      return inputBasis
   end;

endmodule;

tr bb_buchberger;
tr bb_buchberger1;

buchberger({3 * x1^2 * x2^2 - x2*x3^7}, {x1, x2, x3}, lex);

buchberger({3 * x1^2 * x2^2 - x2*x3^7}, {x1, x2, x3}, revgradlex);

end;  % of file
