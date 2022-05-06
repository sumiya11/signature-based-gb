module f5lp;
% The module implements LabeledPolynomial

% The lp module provides the Labeled Polynomial interface --
% a special polynomial type to be used in the f5-style algorithms.
% The LabeledPolynomial object `p` is represented as a 3-item list:
%   {'lp, evaluation of `p`, signature `p`}
% Where Evaluation of `p` is a `Polynomial` object (defined in f5poly file),
%        Signature of `p` if a `Signature` object.
% The `Signature` struct is described in the following.
%
% The LabeledPolynomial interface provides, in particular, procedures
%   . lp_eval(x) - returns evaluation of LabeledPolynomial x, a Polynomial
%   . lp_sgn(x)  - returns signature of LabeledPolynomial x, a Signature

% The Signature object `sgn` is a 3-item list:
%   {'sgn, index of `sgn`, term of `sgn`}
% Where index of `sgn` is an Integer,
%       term of `sgn` is a Term (defined in f5poly file)
%
% The signature interface has procedures for accessing the index and the term:
%   lp_indexSgn and lp_termSgn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct a Signature object from index and term
asserted inline procedure lp_Signature(idx: Integer, st: Term);
  {'sgn, idx, st};

% Instantiates LabeledPolynomial object from the given Polynomial `poly` and
% places garbage in the Signature position
%
% Currently, used only in computing plain polynomial normal form
asserted inline procedure lp_LabeledPolynomial0(
                              poly: Polynomial): LabeledPolynomial;
  lp_LabeledPolynomial1(poly, 0);

% instantiates LabeledPolynomial from Polynomial and the Signature index
asserted inline procedure lp_LabeledPolynomial1(
                              poly: Polynomial,
                              idx: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly, lp_Signature(idx, poly_identityTerm()));

% instantiates LabeledPolynomial from Polynomial and its signature
asserted inline procedure lp_LabeledPolynomial2(
                              poly: Polynomial,
                              sgn: Signature): LabeledPolynomial;
  {'lp, poly, sgn};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% LP & SIGNATURE INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return the signature of f
asserted inline procedure lp_sgn(f: LabeledPolynomial): Signature;
  caddr f;

% Return the evaluation of f
asserted inline procedure lp_eval(f: LabeledPolynomial): Polynomial;
  cadr f;

% set the evaluation of lp to Polynomial ev
asserted inline procedure lp_setEval(lp: LabeledPolynomial, ev: Polynomial);
  cadr lp := ev;

% set the signature of lp to Signature ev
asserted inline procedure lp_setSgn(lp: LabeledPolynomial, s: Signature);
  caddr lp := s;

% return the index of signature s
asserted inline procedure lp_indexSgn(s: Signature): Integer;
  cadr s;

% return the term of signature s
asserted inline procedure lp_termSgn(s: Signature): Term;
  caddr s;

% Signatures s1 == s2 ?
asserted inline procedure lp_eqSgn(s1: Signature, s2: Signature);
  lp_indexSgn(s1) #= lp_indexSgn(s2) and
    poly_eqTerm!?(lp_termSgn(s1), lp_termSgn(s2));

% Zero LabeledPolynomial is represented as
%   {'lp, zero Polynomial, any Signature}
% Check if LabeledPolynomial is zero
asserted inline procedure lp_iszero!?(lp: LabeledPolynomial);
  poly_iszero!?(lp_eval(lp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% SIGNATURE MANIPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multiply signature sgn by the term ev
asserted inline procedure lp_mulSgn(sgn: Signature, ev: Term);
  lp_Signature(lp_indexSgn(sgn), poly_mulTerm(lp_termSgn(sgn), ev));

% compare signatures s1 and s2 with the
% (reversed) Position over term order extension
asserted procedure lp_cmpSgn(s1: Signature, s2: Signature);
  if lp_indexSgn(s1) #= lp_indexSgn(s2) then
    poly_cmpTerm(lp_termSgn(s1), lp_termSgn(s2))
  else
    lp_indexSgn(s1) #< lp_indexSgn(s2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LP COEFFICIENTS MANIPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mainly, these fall back to simple polynomial coefficients operations
% from f5poly

% Normalize the evaluation of `f`
asserted inline procedure lp_normalize(f: LabeledPolynomial): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_normalize(lp_eval(f)), lp_sgn(f));

% Scale denominators of the evaluation of `f`
asserted inline procedure lp_scaleDenominators(f: LabeledPolynomial): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_scaleDenominators(lp_eval(f)), lp_sgn(f));

% Reduce coefficients of the evaluation of `f` modulo `prime`
asserted inline procedure lp_reduceCoeffs(f: LabeledPolynomial,
                                          prime: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_reduceCoeffs(lp_eval(f), prime), lp_sgn(f));

% Reconstruct coefficients of the evaluation of `f` modulo `prime`
asserted inline procedure lp_reconstructCoeffs(f: LabeledPolynomial,
                                          prime: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_reconstructCoeffs(lp_eval(f), prime), lp_sgn(f));

% Given LPs (polyaccum mod modulo) and (polycomp mod prime)
% construct a new LP with the evaluation equal
% to the modular reconstruction of evaluations of polyaccum and polycomp
asserted inline procedure lp_crtCoeffs(polyaccum: Polynomial, modulo: Integer,
                      polycomp: Polynomial, prime: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(
      poly_crtCoeffs(lp_eval(polyaccum), modulo, lp_eval(polycomp), prime),
      lp_sgn(polyaccum)
  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LabeledPolynomial sorting ad-hoc

% return true if lead(eval(lp1)) < lead(eval(lp2))
asserted procedure lp_cmpLPLead(lp1: LabeledPolynomial, lp2: LabeledPolynomial);
  poly_cmpPolyLead(lp_eval(lp1), lp_eval(lp2));

% return true if lead(eval(lp2)) < lead(eval(lp1))
asserted procedure lp_cmpLPLeadReverse(lp1, lp2);
  lp_cmpLPLead(lp2, lp1);

% return true if total_degree(lead(lp1)) < total_degree(lead(lp2))
asserted procedure lp_leadTotalDegreeCmp(lp1: LabeledPolynomial,
                                          lp2: LabeledPolynomial);
  poly_leadTotalDegreeCmp(lp_eval(lp1), lp_eval(lp2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endmodule;

end;  % of file
