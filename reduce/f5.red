
in "lp.red";

module f5;

load_package assert;
on assert;

put('f5, 'psopfn, 'f5_groebner);

asserted procedure f5_groebner(u: List): List;
   begin scalar inputBasis, variables, sortMode, outputBasis;
      if null u then
         f5_error();
      inputBasis := reval pop u;
      if not (pop inputBasis eq 'list) then
         f5_error();
      variables := reval pop u;
      if not (pop variables eq 'list) then
         f5_error();
      sortMode := pop u;
      if not null u then
         f5_error();
      dip_init(variables, sortMode);
      inputBasis := for each f in inputBasis collect
         dip_f2dip numr simp f;
      outputBasis := f5_groebner1(inputBasis);
      outputBasis := 'list . for each f in outputBasis collect
         dip_2a lp_evaluation f;
      return outputBasis
   end;

% Void return type failes assert check
asserted procedure f5_error();
   rederr "usage: buchberger(polynomials: List, variables: List, sortmode: Id)";

asserted procedure f5_selectNext(spolys: List): List;
   begin scalar i, idx, sgn, elem;
      i := 1;
      s := lp_signature(car spolys); % nonempty
      while spolys do <<
         p . spolys := spolys;
         if (i equal 1) or lp_potCompareSignatures(s, sgn) then <<
            idx  := i;
            elem := p;
            sgn  := s
         >>;
         i := i + 1
      >>;

      return idx . elem;
   end;

asserted procedure f5_groebner1(inputBasis: List): List;
   begin scalar basis, spolys, known_syz;
   integer i, ii;
      % form output list..
      basis  := lp_constructModule(inputBasis);

      % and initial s-polynomials
      spolys := nil;
      for each pi in basis do <<
         for each pj in basis do <<
            if not (pi equal pj) then <<
               msi . msj := lp_spolyMultSignatures(pi, pj);

               % do no add redundant pairs
               if lp_potCompareSignatures(msi, msj) then
                  spolys := nconc(spolys, { lp_spoly(pi, pj) })
            >>
         >>
      >>;

      ii := 1;
      while spolys do <<
         i . p := f5_selectNext(spolys);

         % inplace ?
         spolys := remove(spolys, i);

         p_nf := lp_regularNormalForm(p, basis);

         if not lp_isSyzygy(p_nf) then <<
            if (not lp_isSingularlyTopReducible(p_nf, basis)) then <<
               for each gg in basis do <<
                  msi . msj := lp_spolyMultSignatures(p_nf, gg);
                  if not (msi equal msj) then
                     spolys := nconc(spolys, { lp_spoly(p_nf, gg) })
               >>;
               basis := nconc(basis, { p_nf })
            >>
         >>;

         ii := ii + 1
         % if ii > 5 then spolys := nil;
      >>;

      return basis
   end;


trst f5_groebner;
trst f5_groebner1;

trst f5_selectNext;


endmodule;

f5({x1 + x2 + x3, x1*x2 + x2*x3 + x1*x3, x1*x2*x3 - 1}, {x1, x2, x3}, lex);

end;  % of file



