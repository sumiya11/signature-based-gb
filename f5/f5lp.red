module f5lp;

%  The LP module provides Labeled Polynomial interface --
%  a tuple of polynomials for usage in f5-style algorithms.
%  Polynomial `p` is stored internally as a 4-item list:
%     {'lp, evaluation of `p`, signature index of `p`, signature monomial of `p`}
%
% The interface provides, in particular,
% functions lp_evaluation and lp_signature, that return
% second and {third, fouth} items of internal list respectively
%

struct LabeledPolynomial;

% instantiate LabeledPolynomial from Polynomial and leading index
asserted inline procedure lp_LabeledPolynomial1(
                              poly: Polynomial,
                              idx: Integer): LabeledPolynomial;
   lp_LabeledPolynomial2(poly, { idx, poly_zeroExponent() });

% instantiate LabeledPolynomial from Polynomial and signature
asserted inline procedure lp_LabeledPolynomial2(
                              poly: Polynomial,
                              sgn: List): LabeledPolynomial;
   'lbl . poly . sgn;

asserted inline procedure lp_getPoly(lp: LabeledPolynomial): Polynomial;
   cadr lp;

asserted inline procedure lp_getLabel(lp: LabeledPolynomial): List;
   cddr lp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return the signature of module element f
asserted inline procedure lp_signature(f: LabeledPolynomial): List;
   lp_getLabel(f);

% return the evaluation of module element f
asserted inline procedure lp_evaluation(f: LabeledPolynomial): Polynomial;
   lp_getPoly(f);

% for a given set of polynomials F construct
% a standard basis of unit vectors corresponding to F
% in internal representation
asserted procedure lp_constructModule(F: List): List;
   begin scalar f, elem, outputModule;
   integer i;
      outputModule := nil;
      i := 1;
      while F do <<
         ff := car F;
         elem := lp_LabeledPolynomial1(ff, i);
         outputModule := nconc(outputModule, {elem});
         i := i + 1;
         F := cdr F
      >>;
      return outputModule
   end;

% multiply signature sgn by monomial exponent vector ev
asserted inline procedure lp_multSignature(sgn: List, ev: List);
   begin scalar idx, term;
      idx . term := sgn;
      return { idx, poly_sumExponents(term, ev) };
   end;

% compare signatures sgn1 vs. sgn2 with the
% position over term extension
asserted procedure lp_potCompareSignatures(sgn1: List, sgn2: List): Boolean;
   begin integer idx1, idx2;
   scalar ev1, ev2;
      idx1 . ev1 := sgn1;
      idx2 . ev2 := sgn2;
      return if idx1 equal idx2 then poly_cmpExponents(ev1, ev2)
         else (idx1 > idx2)
   end;

% minimal signature based on pot comparison strategy
asserted inline procedure lp_potMinSignature(sgn1: List, sgn2: List): List;
   << if lp_potCompareSignatures(sgn1, sgn2) then sgn1 else sgn2 >>;

% maximal signature based on pot comparison strategy
asserted inline procedure lp_potMaxSignature(sgn1: List, sgn2: List): List;
   << if lp_potCompareSignatures(sgn1, sgn2) then sgn2 else sgn1 >>;

% compare LabeledPolynomials as module elems lp1 vs. lp2 with
% position over term extension
asserted procedure lp_potCompareLP(lp1: LabeledPolynomial, lp2: LabeledPolynomial): Boolean;
   begin scalar s1, s2;
      sgn1   := lp_signature(lp1);
      sgn2   := lp_signature(lp2);
      lp_potCompareSignatures(sgn1, sgn2)
   end;

% checks if lp represents a syzygy
asserted inline procedure lp_isSyzygy(lp: LabeledPolynomial): Boolean;
   lp_evaluation(lp) equal nil;

% extract an array of exponent vectors of f
% TODO: remove this
asserted procedure lp_extractDpExponents(f);
   begin scalar ans;
      while f do <<
         ans := nconc(ans, {car f});
         f   := cddr f >>;
      return ans
   end;

% try to reduce module elem f with module elem g regularly,
% returns success flag with the result of reduction
asserted procedure lp_tryRegularReduce1(f: LabeledPolynomial, g: LabeledPolynomial): List;
   begin scalar fexps, fexp, fcoef, gcoeflead, gexplead,
                  m, evlf, evlg, sgnf, sgng, idxf, idxg, tf, tg,
                  anspoly, anssgn, tmpflag;
   integer nterms, i;
      evlf := lp_evaluation(f);
      evlg := lp_evaluation(g);
      sgnf := lp_signature(f);
      sgng := lp_signature(g);

      gexplead  := poly_leadExponent(evlg);
      gcoeflead := poly_leadCoeff(evlg);

      nterms := poly_length(evlf);

      anspoly := nil;
      anssgn := nil;

      flag := nil;

      tmpevlf := evlf;

      for i := 1:nterms do <<
         % TODO: rewrite this
         fexp    := car tmpevlf;
         fcoef   := cadr tmpevlf;
         tmpevlf := cddr tmpevlf;
         % check that evaluations can be reduced
         % ORDER of arguments !!!!!
         tmpflag := poly_dividesExponents(gexplead, fexp);
         if (not flag) and tmpflag then <<
            flag := t;
            m := poly_difExponents(fexp, gexplead);
            % check that reduction is regular
            if not (sgnf equal lp_multSignature(sgng, m)) then <<
               % new polynomial instance
               % TODO: abstract this a bit
               anspoly := dip_ilcomb(evlf, bc_neg(gcoeflead), ev_zero(), evlg, fcoef, m);
               % new signature instance
               newexp  := for each x in cadr sgnf collect x;
               anssgn  := car sgnf . newexp . nil
            >>
         >>
      >>;

      % TODO: Construction of a LabeledPolynomial
      return flag . lp_LabeledPolynomial2(anspoly, anssgn)
   end;

% try to reduce module elem f with module G regurarly
% returns success flag with the result of reduction
asserted procedure lp_tryRegularReduce(f: LabeledPolynomial, G: List): List;
   begin scalar ans, gg;
      if null G then
         ans := nil . f
      else if null (car f) then
         ans := nil . f
      else <<
         gg := car G;
         ans := lp_tryRegularReduce1(f, gg);
      >>;

      return if not (car ans) then
         if not (null G) then lp_tryRegularReduce(f, cdr G) else ans
      else ans
   end;

asserted procedure lp_regularNormalForm(
                              f: LabeledPolynomial,
                              G: List): LabeledPolynomial;
   begin scalar res;
      state := t . f;
      while car state do <<
         state := lp_tryRegularReduce(cdr state, G);
      >>;
      return cdr state
   end;

asserted procedure lp_isSingularlyTopReducible(f: LabeledPolynomial, G: List): Boolean;
   begin scalar evlf, evlg, tf, tg, leadf, leadg, m, ans;
      evlf  := lp_evaluation(f);
      leadf := poly_leadExponent(evlf);
      tf := lp_signature(f);
      tg := lp_signature(g);
      for each gg in G do <<
         evlg := lp_evaluation(gg);
         leadg := poly_leadExponent(evlg);
         if poly_dividesExponents(leadg, leadf) then <<
            m := difExponents(leadf, leadg);
            if tf equal lp_multSignature(tg, m) then ans := t
         >>
      >>;
      return ans
   end;

% compute cofactors of Spoly of p1 and p2
asserted procedure lp_spolyCofactors(f: LabeledPolynomial, g: LabeledPolynomial): Boolean;
   begin scalar p1, p2, e1, e2, elcm, mult1, mult2;
      p1 := lp_evaluation(f);
      p2 := lp_evaluation(g);

      e1 := poly_leadExponent(p1);
      e2 := poly_leadExponent(p2);

      elcm := poly_lcmExponents(e1, e2);

      mult1 := poly_difExponents(elcm, e2);
      mult2 := poly_difExponents(elcm, e1);

      return mult1 . mult2
   end;

% compute signatures of pi and pj
% after multiplication by cofactors of Spolynomial of pi and pj
asserted procedure lp_spolyMultSignatures(p1: LabeledPolynomial, pj: LabeledPolynomial): List;
   begin scalar mij, mi, mj, si, sj, msi, msj;
      mi . mj := lp_spolyCofactors(p1, pj);

      si := lp_signature(p1);
      sj := lp_signature(pj);

      return lp_multSignature(si, mj) . lp_multSignature(sj, mi)
   end;

% compute Spoly of f and g, where f and g are tuples in internal representation
asserted procedure lp_spoly(f: LabeledPolynomial, g: LabeledPolynomial): LabeledPolynomial;
   begin scalar p1, p2, e1, e2, c1, c2, elcm, mult1, mult2, ans,
                  sgn1, sgn2, sgn11, sgn22, anssgn;

      % here, we may that signatures do not cancel each other,
      % so the updated signature can be computed easily

      p1 := lp_evaluation(f);
      p2 := lp_evaluation(g);

      e1 := poly_leadExponent(p1);
      e2 := poly_leadExponent(p2);

      elcm := poly_lcmExponents(e1, e2);

      c1 := poly_leadCoeff(p1);
      c2 := poly_leadCoeff(p2);

      mult1 := poly_difExponents(elcm, e2);
      mult2 := poly_difExponents(elcm, e1);

      % TODO: abstract this a bit
      ans := dip_ilcomb(p2, c1, mult1, p1, bc_neg(c2), mult2);

      sgn1 := lp_signature(f);
      sgn2 := lp_signature(g);

      sgn11 := lp_multSignature(sgn1, mult2);
      sgn22 := lp_multSignature(sgn2, mult1);

      anssgn := lp_potMaxSignature(sgn11, sgn22);

      return lp_LabeledPolynomial2(ans, anssgn)
   end;


% trst lp_LabeledPolynomial;
% trst lp_spoly;
% trst lp_spolyMultSignatures;
% trst lp_spolyCofactors;
% trst lp_multSignature;

% trst lp_tryRegularReduce1;
% trst lp_tryRegularReduce;

% untrst lp_LabeledPolynomial;

endmodule;


end;  % of file
