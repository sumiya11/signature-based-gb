% F5 algorithm implementation module
%
% The module provides Groebner basis computation routine `f5` based on
% Fougeres F5 algorithm
module f5;

create!-package('(f5 f5lp f5poly), nil);

load_package assert;
on1 'assert;

put('f5, 'psopfn, 'f5_groebner);

struct Polynomial;
struct LabeledPolynomial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure f5_groebner(u: List): List;
   begin scalar inputBasis, variables, sortMode, isModular, outputBasis;
      if null u then
         f5_error();
      inputBasis := reval pop u;
      if not (pop inputBasis eq 'list) then
         f5_error();
      variables := reval pop u;
      if not (pop variables eq 'list) then
         f5_error();
      sortMode := pop u;
      isModular := pop u;
      if not null u then
         f5_error();
      poly_initRing(variables, sortMode);
      inputBasis := for each f in inputBasis collect
         poly_f2poly numr simp f;

      !*f5modular := isModular;
      if !*f5modular then
         outputBasis := f5_groebnerModular1(inputBasis)
      else
         outputBasis := f5_groebner1(inputBasis);

      outputBasis := 'list . for each f in outputBasis collect
                        poly_poly2a f;
      return outputBasis
   end;

% Void return type fails assert check
asserted procedure f5_error();
   rederr "usage: buchberger(polynomials: List, variables: List, sortmode: Id)";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MODULAR F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure f5_correctnessCheck(reconstructedBasis: List): List;
   t;

asserted procedure f5_modularReduction(inputBasis: List): List;
   begin scalar prime, ans, poly;
      ans := nil;
      prime := current!-modulus;
      while inputBasis do <<
         poly . inputBasis := inputBasis;
         poly := poly_reduceCoeffs(poly, prime);
         ans := poly . ans
      >>;
      return reversip(ans)
   end;

asserted procedure f5_rationalReconstruction(inputBasis: List): List;
   begin scalar prime;
      ans := nil;
      prime := current!-modulus;
      while inputBasis do <<
         poly . inputBasis := inputBasis;
         poly := poly_reconstructCoeffs(poly, prime);
         ans := poly . ans
      >>;
      return reversip(ans)
   end;

% TODO
asserted procedure f5_scaleDenominatorsInplace(inputBasis: List): List;
   begin scalar prime;
      ans := nil;
      prime := current!-modulus;
      while inputBasis do <<
         poly . inputBasis := inputBasis;
         poly := poly_reconstructCoeffs(poly, prime);
         ans := poly . ans
      >>;
      return reversip(ans)
   end;

asserted procedure f5_groebnerModular1(inputBasis: List): List;
   begin scalar reducedBasis, computedBasis, reconstructedBasis,
               correctness;
      correctness := nil;
      prime := initial_prime!*;

      % hmm ?
      % f5_scaleDenominatorsInplace(inputBasis);
      % now all denominators are 1

      % so, while the correctness check failes
      while not correctness do <<
         % select next prime
         prime := primes_nextLuckyPrime(inputBasis, prime);

         reducedBasis := f5_modularReduction(inputBasis);

         computedBasis := f5_groebner1(reducedBasis);

         % CRT (prevprime, prime) --> (prevprime*prime)

         % reconstruct modulo prevprime*prime
         reconstructedBasis := f5_rationalReconstruction(computedBasis);

         correctness := f5_correctnessCheck(reconstructedBasis);
      >>;

      return reconstructedBasis
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERIC F5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for a given set of polynomials F construct
% a standard basis of unit vectors corresponding to F
asserted procedure f5_constructModule(F: List): List;
   begin scalar f, elem, outputModule;
   integer i;
      outputModule := nil;
      i := 1;
      while F do <<
         ff := car F;
         elem := lp_LabeledPolynomial1(ff, i);
         outputModule := elem . outputModule;
         i := i + 1;
         F := cdr F
      >>;
      return outputModule
   end;

% for a given set of LabeledPolynomials F construct
% the set of all principal syzygies
asserted procedure f5_constructSyzygies(F: List): List;
   begin scalar syzygies;
   integer i;
      syzygies := nil;
      n := length(F);
      for i := 1:n do <<
         FF := cdr F;
         for j := i+1:n do <<
            p1 := car F;
            p2 := car FF;
            syzygies := lp_principalSyzygy(p1, p2) . syzygies;
            FF := cdr FF
         >>;
         F := cdr F
      >>;

      return syzygies
   end;

% for a given set of LabeledPolynomials F construct
% the set of all S-polynomials
asserted procedure f5_constructSpolys(F: List): List;
   begin scalar spolys;
      spolys := nil;

      n := length(F);
      for i := 1:n do <<
         FF := cdr F;
         for j := i+1:n do <<
            p1 := car F;
            p2 := car FF;
            if not poly_disjLead!?(lp_evaluation(p1), lp_evaluation(p2)) then
               spolys := lp_spoly(p1, p2) . spolys;
            FF := cdr FF
         >>;
         F := cdr F
      >>;

      return spolys
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure f5_syzygyCriterion(p: LabeledPolynomial, syzygies: List);
   begin scalar flag;
      flag := nil;
      if null syzygies then
         flag := nil
      else <<
         syz := lp_signature(car syzygies);
         sgn := lp_signature(p);

         if lp_signatureDivides(syz, sgn) then
            flag := t
         else
            flag := f5_syzygyCriterion(p, cdr syzygies)
      >>;
      return flag
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure f5_selectNext(spolys: List): List;
   begin scalar i, idx, sgn, elem;
      i := 1;
      s := lp_signature(car spolys); % nonempty
      while spolys do <<
         p . spolys := spolys;
         if (i equal 1) or lp_potCmpSignature(s, sgn) then <<
            idx  := i;
            elem := p;
            sgn  := s
         >>;
         i := i + 1
      >>;

      return idx . elem
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% given a Groebner basis
% returns the unique Groebner basis of the corresponding ideal
%
% Output contracts:
%  . the basis is interreduced
%  . the basis is normalized by leading coefficients
%  . the basis is sorted by leading terms increasingly
asserted procedure f5_standardizeOutput(basis: List): List;
   begin scalar iter, reducers, tobereduced, reduced, poly, output;
      sortedbasis := sort(basis, 'poly_cmpLPLead);

      output := nil;
      iter   := sortedbasis;
      reducers  := cdr sortedbasis;
      while iter do <<
         tobereduced . iter := iter;
         reduced := lp_normalForm(tobereduced, reducers, nil);
         if not (lp_isSyzygy reduced) then <<
            poly     := lp_evaluation(reduced);
            poly     := poly_normalize(poly);
            output   := poly . output;
            reducers := append(reducers, { reduced })
         >>;
         reducers := cdr reducers
      >>;

      return sort(output, 'poly_cmpPolyLead)
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asserted procedure f5_groebner1(inputBasis: List): List;
   begin scalar basis, spolys, known_syz;
   integer i, ii;
      % form initial basis list
      basis    := f5_constructModule(inputBasis);

      % form principal syzygies
      syzygies := f5_constructSyzygies(basis);

      % construct S-polynomials
      spolys   := f5_constructSpolys(basis);

      ii := 1;
      while spolys do <<
         % select one S-poly
         i . p := f5_selectNext(spolys);

         % inplace ?
         spolys := remove(spolys, i);

         % if syzygy criterion then discard
         if f5_syzygyCriterion(p, syzygies) then
            % pass
         else <<
            % signature-safe reduction
            p_nf := lp_normalForm(p, basis, t);

            if lp_isSyzygy(p_nf) then
               syzygies := p_nf . syzygies
            else if (not lp_isSingularlyTopReducible(p_nf, basis)) then <<
               for each gg in basis do <<
                  msi . msj := lp_spolyMultSignatures(p_nf, gg);
                  if not (msi equal msj) then
                     spolys := lp_spoly(p_nf, gg) . spolys
               >>;
               basis := p_nf . basis
            >>

         >>;

         ii := ii + 1;
         if ii > 1 then spolys := nil
      >>;

      basis := f5_standardizeOutput(basis);

      return basis
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trst f5_groebner;

% trst f5_selectNext;
% trst f5_constructSyzygies;
% trst f5_constructSpolys;
% trst f5_syzygyCriterion;
% trst f5_standardizeOutput;

% trst f5_groebner1;
% trst f5_modularReduction;
% trst f5_rationalReconstruction;

endmodule;

% f5({x1 + x2, x1*x2 + 1}, {x1, x2}, lex, nil);

end;  % of file
