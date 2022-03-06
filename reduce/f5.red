
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
         dip_2a f;
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
   begin scalar basis, spolys;
   integer i;

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
               print {"pi, pj", pi, pj};
               mij := mod_spolyCofactors(pi, pj);
               mi := car mij;
               mj := cdr mij;

               si := mod_signature(pi);
               sj := mod_signature(pj);

               msi := mod_multSignature(si, mj);
               msj := mod_multSignature(sj, mi);
               print {"mi, mj", mi, mj};

               print {"si, sj, msi, msj", si, sj, msi, msj};
               print {"passing pi, pj", pi, pj};

               % do no add redundant pairs
               if mod_potCompareSignatures(msi, msj) then
                  spolys := nconc(spolys, { mod_spoly(pi, pj) })
            >>
         >>
      >>;

      print {"GENERATED SPOLYS", spolys};
      print {"MAIN CYCLE START"};

      i := 1;
      while spolys do <<
         ip := f5_selectNext(spolys);
         i  := car ip;
         p  := cdr ip;

         % inplace ?
         spolys := remove(spolys, i);

         print "------------------";
         print {"signature ", mod_signature(p)};
         print {"iteration ", i, p};
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
                  mi, mj := mod_spolyCofactors(p_nf, gg);
                  si := mod_signature(p_nf);
                  sj := mod_signature(gg);

                  msi := mod_multSignature(si, mi);
                  msj := mod_multSignature(sj, mj);

                  if not (msi equal msj) then
                     spolys := nconc(spolys, { mod_spoly(p_nf, gg) })
               >>;
               basis := nconc(basis, p_nf)
            >>
         >>;

         i := i + 1;
         % if i > 5 then spolys := nil;
      >>;

      return basis
   end;

endmodule;

tr f5_groebner;
tr f5_groebner1;

in "mod.red";

groebner({8*x1 + x2, x1*x2 + 4}, {x1, x2}, lex);

end;  % of file


