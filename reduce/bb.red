
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

% compute Spoly of p1 and p2,
% nil is returned if zeroed
procedure Spoly(p1, p2);
   begin scalar e1, e2, c1, c2, elcm, mult1, mult2, ans;
      e1 := dip_evlmon(p1);
      e2 := dip_evlmon(p2);

      elcm := ev_lcm(e1, e2);

      c1 := dip_lbc(p1);
      c2 := dip_lbc(p2);

      mult1 := ev_dif(elcm, e2);
      mult2 := ev_dif(elcm, e1);

      % dip_ilcomb(p1,c1,t1,p2,c2,t2);
      % Compute p1*c1^t1+p2*c2^t2.
      ans := dip_ilcomb(p2, c1, mult1, p1, bc_neg(c2), mult2);

      return ans
   end;

% find an element from basis with lead that divides lead of f
procedure FindTopReducer(f, basis);
   begin scalar ans;
      if null basis then
         ans := nil
      else
      begin scalar ef, g, eg, flag, res;
         g := car basis;
         eg := dip_evlmon(g);
         ef := dip_evlmon(f);
         print {"lead check", ef, eg};
         flag := ev_divides!?(eg, ef);   % order?
         print {"flag", flag};
         if not flag then
            ans := FindTopReducer(f, cdr basis)
         else
            ans := g
      end;
      return ans
   end;

procedure TopReduce(f, g);
   begin
      return Spoly(f, g)
   end;

procedure NormalForm(f, basis);
   begin scalar flag;
      flag := t;
      while flag do
      begin scalar g;
         g := FindTopReducer(f, basis);
         print {"reducer", g};
         if null g then
            flag := nil
         else
         begin
            f := TopReduce(f, g);
            print {"reduced", f};
            if null f then
               flag := nil;
         end;
      end;
      return f
   end;

asserted procedure bb_buchberger1(inputBasis: List): List;
   begin scalar basis, spolys; integer i;
      % form output list and initial s-polynomials
      basis := inputBasis;    % non-copy assign
      spolys := for each i in inputBasis join
                  for each j in inputBasis collect
                     if not (i equal j) then Spoly(i, j); % not all are assigned
      % some duplicates in spolys, but it's okay

      while spolys do
      begin scalar p, p_nf;
         p := car spolys;
         spolys := cdr spolys;
         print "------------------";
         print {"iteration ", i, p};
         print {"basis = ", basis, "spolys = ", spolys};
         if not null p then
         begin
            p_nf := NormalForm(p, basis);
            print {"normal form", p_nf};
            if not null p_nf then
            begin scalar new_spolys;
               basis := nconc(basis, {p_nf});
               new_spolys := for each g in basis collect Spoly(g, p_nf);
               spolys := nconc(spolys, new_spolys);
            end;
         end;
         i := i + 1;
         % if i > 5 then spolys := nil;
      end;

      return basis
   end;

endmodule;

tr bb_buchberger;
tr bb_buchberger1;

% buchberger({3 * x1^2 * x2^2 - x2*x3^7}, {x1, x2, x3}, lex);

% buchberger({3 * x1^2 * x2^2 - x2*x3^7}, {x1, x2, x3}, revgradlex);

% buchberger({x1^2 * x2 - 2 * x2, x2^2 + 1}, {x1, x2}, lex);

buchberger({x1 + x2, x1*x2 + 1}, {x1, x2}, lex);

end;  % of file


