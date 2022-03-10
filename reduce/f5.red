module f5;

off1 'allfac;

load!-package 'dp;

put('groebner, 'psopfn, 'f5_groebner);

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
         dip_2a mod_evaluation f;
      return outputBasis
   end;

asserted procedure f5_error(): Void;
   rederr "usage: buchberger(polynomials: List, variables: List, sortmode: Id)";

procedure f5_selectNext(spolys);
   begin scalar i, idx, sgn, elem;
      i := 1;
      while spolys do <<
         p := car spolys;
         spolys := cdr spolys;
         s := mod_signature(p);
         if (i equal 1) or mod_potCompareSignatures(s, sgn) then <<
            idx  := i;
            elem := p;
            sgn  := s
         >>;
         i := i + 1
      >>;

      print {"IN SELECT idx, elem, sgn", idx, elem, " ---> ", sgn};

      return idx . elem;
   end;

procedure f5_groebner1(inputBasis);
   begin scalar basis, spolys, known_syz;
   integer i, ii;

      print {"input", inputBasis};

      % form output list..
      basis  := mod_constructModule(inputBasis);

      print {"basis", basis};
      print {"firrst", car basis};

      % and initial s-polynomials
      spolys := nil;
      for each pi in basis do <<
         for each pj in basis do <<
            if not (pi equal pj) then <<
               msij := mod_spolyMultSignatures(pi, pj);
               msi  := car msij;
               msj  := cdr msij;

               % do no add redundant pairs
               if mod_potCompareSignatures(msi, msj) then
                  spolys := nconc(spolys, { mod_spoly(pi, pj) })
            >>
         >>
      >>;

      print {"GENERATED SPOLYS", spolys};
      print {"MAIN CYCLE START"};

      ii := 1;
      while spolys do <<
         ip := f5_selectNext(spolys);
         i  := car ip;
         p  := cdr ip;

         % inplace ?
         spolys := remove(spolys, i);

         print "------------------";
         print {"signature ", mod_signature(p)};
         print {"iteration ", ii, p};
         print {"basis = ", basis};
         print {"spolys = ", spolys};
         print {"left spolys", length(spolys)};

         p_nf := mod_regularNormalForm(p, basis);
         print {"normal form", p_nf};

         if not mod_isSyzygy(p_nf) then <<
            print {"Not a syzygy"};
            if (not mod_isSingularlyTopReducible(p_nf, basis)) then <<
               print {"Not top reducible"};
               for each gg in basis do <<
                  print {"!!!", p_nf, gg};
                  msij := mod_spolyMultSignatures(p_nf, gg);
                  msi  := car msij;
                  msj  := cdr msij;

                  print {"???"};

                  if not (msi equal msj) then
                     spolys := nconc(spolys, { mod_spoly(p_nf, gg) })
               >>;
               basis := nconc(basis, { p_nf })
            >>
         >>;

         ii := ii + 1
         % if ii > 5 then spolys := nil;
      >>;

      return basis
   end;

endmodule;

tr f5_groebner;
tr f5_groebner1;

in "mod.red";

% groebner({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);
% groebner({x1*x2 + 1, x2*x3 + 1}, {x1, x2, x3}, lex);
groebner({x1 + x2 + x3, x1*x2 + x2*x3 + x1*x3, x1*x2*x3 - 1}, {x1, x2, x3}, lex);
groebner({10*x1*x2^2 - 11*x1 + 10, 10*x1^2*x2 - 11*x2 + 10}, {x1, x2}, revgradlex);

end;  % of file
