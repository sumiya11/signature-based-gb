
module mod;

%  The Module provides module over polynomial ring element interface --
%  a tuple of polynomials for usage in f5 algorithms.
%  Tuple `p` is stored internally as a 3-item list:
%     {evaluation of `p`, signature index of `p`, signature monomial of `p`}
%
% The interface provides, in particular,
% functions mod_evaluation and mod_signature, that return
% first and {second, third} items of internal list respectively
%

off1 'allfac;

load!-package 'dp;

% for a given set of polynomials F construct
% a standard basis of unit vectors corresponding to F
% in internal representation
procedure mod_constructModule(F);
   begin scalar i, f, elem, outputModule;
      outputModule := nil;
      i := 1;
      while F do <<
         ff := car F;
         elem := mod_poly2unit(ff, i);
         outputModule := nconc(outputModule, {elem});
         i := i + 1;
         F := cdr F
      >>;
      return outputModule
   end;

% initialize module element (0, .., f, 0, .., 0),
% where f is at i-th position
procedure mod_poly2unit(f, i);
   begin scalar moduleElement, lead;
      moduleElement := f . i . ev_zero() . nil;
      return moduleElement
   end;

% return the signature of module element f
procedure mod_signature(f);
   cdr f;

% return the evaluation of module element f
procedure mod_evaluation(f);
   car f;

% multiply signature s by monomial vector e
procedure mod_multSignature(s, e);
   begin scalar i, ti;
      i := car s;
      ti := cadr s;
      print ti;
      print e;
      print {"mutliply", ti, e};
      return i . ev_sum(ti, e) . nil
   end;

% compare signatures s1 vs. s2 with
% position over term extension
procedure mod_potCompareSignatures(s1, s2);
   begin integer idx1, idx2;
   scalar ev1, ev2;
      print {"inside compare", s1, s2};
      idx1 := car s1;
      idx2 := car s2;
      ev1  := cadr s1;
      ev2  := cadr s2;
      return if idx1 equal idx2 then ev_comp(ev1, ev2)
         else (idx1 < idx2)
   end;

% minimal signature based on pot comparison strategy
procedure mod_potMinSignature(s1, s2);
   << if mod_potCompareSignatures(s1, s2) then s1 else s2 >>;

% maximal signature based on pot comparison strategy
procedure mod_potMaxSignature(s1, s2);
   << if mod_potCompareSignatures(s1, s2) then s2 else s1 >>;

% compare module elems f1 vs. f2 with
% position over term extension
procedure mod_potCompare(f1, f2);
   begin scalar s1, s2;
      s1   := mod_signature(f1);
      s2   := mod_signature(f2);
      mod_potCompareSignatures(s1, s2)
   end;

% check if f is a syzygy
procedure mod_isSyzygy(f);
   mod_evaluation(f) equal nil;

% extract an array of exponent vectors of f
procedure mod_extractDpExponents(f);
   begin scalar ans;
      while f do <<
         ans := nconc(ans, {car f});
         f   := cddr f >>;
      return ans
   end;

% try to reduce module elem f with module elem g regurarly
% returns success flag with the result of reduction
procedure mod_tryRegularReduce1(f, g);
   begin scalar fexps, fexp, fcoef, gcoeflead, gexplead,
                  m, evlf, evlg, sgnf, sgng, idxf, idxg, tf, tg,
                  anspoly, anssgn;
   integer nterms, i;
      evlf := mod_evaluation(f);
      evlg := mod_evaluation(g);
      sgnf := mod_signature(f);
      sgng := mod_signature(g);
      idxf := car sgnf;
      idxg := car sgng;
      tf   := cadr sgnf;
      tg   := cadr sgng;

      gexplead  := dip_evlmon(evlg);
      gcoeflead := dip_lbc(evlg);

      nterms := dip_length(evlf);

      print {"2 try reduce ", f , " with ", g};

      anspoly := f;
      anssgn := sgnf;

      for i := 1:nterms do <<
         fexp  := car evlf;
         fcoef := cadr evlf;
         evlf  := cddr evlf;
         % check that evaluations can be reduced
         % ORDER of arguments !!!!!
         flag := ev_divides!?(gexplead, fexp);
         if flag then <<
            m := ev_dif(fexp, gexplead);
            % check that reduction is regular
            if not (tf equal mod_multSignature(tg, m)) then <<
               anspoly := dip_ilcomb(evlf, bc_neg(gcoeflead), ev_zero(), evlg, fcoef, m);
               anssgn  := sgnf
            >>
         >>
      >>;

      print {"2 result ", flag, anspoly};

      return flag . anspoly . anssgn . nil
   end;

% try to reduce module elem f with module G regurarly
% returns success flag with the result of reduction
procedure mod_tryRegularReduce(f, G);
   begin scalar ans, gg;
      print {"1 try reduce ", f, "  with  ", G};
      if null G then
         ans := nil . f
      else <<
         gg := car G;
         ans := mod_tryRegularReduce1(f, gg);
         print {"ans", ans}
      >>;

      print {"1 result ", ans, car ans};

      return if not (car ans) then
         if not (null G) then mod_tryRegularReduce(f, cdr G) else ans
      else ans
   end;

procedure mod_regularNormalForm(f, G);
   begin scalar res;
      state := t . f;
      while car state do <<
         state := mod_tryRegularReduce(cdr state, G);
         print {"state" , state}
      >>;
      return cdr state
   end;

procedure mod_isSingularlyTopReducible(f, G);
   begin scalar evlf, evlg, tf, tg, leadf, leadg, m, ans;
      evlf  := mod_evaluation(f);
      leadf := dip_evlmon(evlf);
      tf := mod_signature(f);
      tg := mod_signature(g);
      for each gg in G do <<
         evlg := mod_evaluation(gg);
         leadg := dip_evlmon(evlg);
         if ev_divides!?(leadg, leadf) then <<
            m := ev_dif(leadf, leadg);
            if tf equal mod_multSignature(tg, m) then ans := t
         >>
      >>;
      return ans
   end;

% compute cofactors of Spoly of p1 and p2
procedure mod_spolyCofactors(f, g);
   begin scalar p1, p2, e1, e2, elcm, mult1, mult2;
      p1 := mod_evaluation(f);
      p2 := mod_evaluation(g);

      e1 := dip_evlmon(p1);
      e2 := dip_evlmon(p2);

      elcm := ev_lcm(e1, e2);

      mult1 := ev_dif(elcm, e2);
      mult2 := ev_dif(elcm, e1);

      return mult1 . mult2
   end;

% compute Spoly of f and g, where f and g are tuples in internal representation
procedure mod_spoly(f, g);
   begin scalar p1, p2, e1, e2, c1, c2, elcm, mult1, mult2, ans,
                  sgn1, sgn2, sgn11, sgn22, anssgn;

      % here, we may that signatures do not cancel each other,
      % so the updated signature can be computed easily

      p1 := mod_evaluation(f);
      p2 := mod_evaluation(g);

      e1 := dip_evlmon(p1);
      e2 := dip_evlmon(p2);

      elcm := ev_lcm(e1, e2);

      c1 := dip_lbc(p1);
      c2 := dip_lbc(p2);

      mult1 := ev_dif(elcm, e2);
      mult2 := ev_dif(elcm, e1);

      ans := dip_ilcomb(p2, c1, mult1, p1, bc_neg(c2), mult2);

      sgn1 := mod_signature(f);
      sgn2 := mod_signature(g);

      sgn11 := mod_multSignature(sgn1, mult1);
      sgn22 := mod_multSignature(sgn2, mult2);

      print {"sgns", sgn11, sgn22};

      anssgn := mod_potMaxSignature(sgn11, sgn22);

      return ans . car anssgn . cdr anssgn;
   end;


endmodule;


end;  % of file


