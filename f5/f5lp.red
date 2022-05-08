module f5lp;
% The module implements LabeledPolynomial interface.

% The lp module provides the Labeled Polynomial interface --
% a special polynomial type to be used in the f5-style algorithm.
% Essentially, LabeledPolynomial object consists of two parts:
% a polynomial part, and a signature part.
%
% The LabeledPolynomial object `p` is represented as a 3-item list:
%   {'lp, evaluation of `p`, signature `p`}
% Where Evaluation of `p` is a `Polynomial` object (defined in f5poly file),
%        Signature of `p` is a `Signature` object (described somewhere below).
%
% LabeledPolynomial interface provides, in particular, procedures
%   . lp_eval(x) - returns the evaluation of LabeledPolynomial x
%   . lp_sgn(x)  - returns the signature of LabeledPolynomial x

% The Signature object `sgn` is a 3-item list:
%   {'sgn, index of `sgn`, term of `sgn`}
% Where the index of `sgn` is an Integer, and
%       the term of `sgn` is a Term (defined in f5poly file)
%
% The signature interface provides procedures for accessing
% the index and the term respectively:
%   lp_indexSgn and lp_termSgn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constructs a Signature from the given index and term
asserted inline procedure lp_Signature(idx: Integer, st: Term): Signature;
  {'sgn, idx, st};

% Instantiates LabeledPolynomial object from the given Polynomial `poly` and
% places garbage in the Signature position.
% Currently, used only in computing plain polynomial normal form
asserted inline procedure lp_LabeledPolynomial0(
                              poly: Polynomial): LabeledPolynomial;
  lp_LabeledPolynomial1(poly, 0);

% Instantiates LabeledPolynomial from Polynomial and the Signature index
asserted inline procedure lp_LabeledPolynomial1(
                              poly: Polynomial,
                              idx: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly, lp_Signature(idx, poly_identityTerm()));

% Instantiates LabeledPolynomial from Polynomial and its signature
asserted inline procedure lp_LabeledPolynomial2(
                              poly: Polynomial,
                              sgn: Signature): LabeledPolynomial;
  {'lp, poly, sgn};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% LP & SIGNATURE INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns the signature of f
asserted inline procedure lp_sgn(f: LabeledPolynomial): Signature;
  caddr f;

% Returns the evaluation of f
asserted inline procedure lp_eval(f: LabeledPolynomial): Polynomial;
  cadr f;

% Sets the evaluation of lp to ev
asserted inline procedure lp_setEval(lp: LabeledPolynomial, ev: Polynomial);
  cadr lp := ev;

% Sets the signature of lp to ev
asserted inline procedure lp_setSgn(lp: LabeledPolynomial, s: Signature);
  caddr lp := s;

% Returns the index of signature s
asserted inline procedure lp_indexSgn(s: Signature): Integer;
  cadr s;

% Returns the term of signature s
asserted inline procedure lp_termSgn(s: Signature): Term;
  caddr s;

% Checks if signatures s1 and s2 are equal
asserted inline procedure lp_eqSgn(s1: Signature, s2: Signature): Boolean;
  lp_indexSgn(s1) #= lp_indexSgn(s2) and
    poly_eqTerm!?(lp_termSgn(s1), lp_termSgn(s2));

% Checks if `lp` is zero as a LabeledPolynomial.
% Zero LabeledPolynomial is represented as
%   {'lp, zero `Polynomial`, any `Signature`}
% so we just check if the polynomial is zero.
asserted inline procedure lp_iszero!?(lp: LabeledPolynomial): Boolean;
  poly_iszero!?(lp_eval(lp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% SIGNATURE MANIPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multiplies signature `sgn` by term `ev`
asserted inline procedure lp_mulSgn(sgn: Signature, ev: Term): Signature;
  lp_Signature(lp_indexSgn(sgn), poly_mulTerm(lp_termSgn(sgn), ev));

% Compares signatures `s1` and `s2` with the
% (reversed) Position over term order extension
asserted procedure lp_cmpSgn(s1: Signature, s2: Signature): Boolean;
  if lp_indexSgn(s1) #= lp_indexSgn(s2) then
    poly_cmpTerm(lp_termSgn(s1), lp_termSgn(s2))
  else
    lp_indexSgn(s1) #< lp_indexSgn(s2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LP COEFFICIENTS MANIPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mainly, these functions fall back to simple polynomial coefficients
% manipulations from f5poly

% Normalizes the evaluation of `f`
asserted inline procedure lp_normalize(f: LabeledPolynomial): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_normalize(lp_eval(f)), lp_sgn(f));

% Scales coefficients of the evaluation of `f` by their common denominator
asserted inline procedure lp_scaleDenominators(f: LabeledPolynomial): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_scaleDenominators(lp_eval(f)), lp_sgn(f));

% Reduces coefficients of the evaluation of `f` modulo `prime`
asserted inline procedure lp_reduceCoeffs(f: LabeledPolynomial,
                                          prime: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_reduceCoeffs(lp_eval(f), prime), lp_sgn(f));

% Reconstructs coefficients of the evaluation of `f` modulo `prime`
asserted inline procedure lp_reconstructCoeffs(f: LabeledPolynomial,
                                               prime: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(poly_reconstructCoeffs(lp_eval(f), prime), lp_sgn(f));

% Given two LPs (polyaccum and polycomp) and two Integers (modulo and prime)
% constructs a new LP with the evaluation equal to the modular reconstruction
% of evaluations of polyaccum and polycomp modulo modulo and prime respectively.
asserted inline procedure lp_crtCoeffs(polyaccum: Polynomial, modulo: Integer,
                      polycomp: Polynomial, prime: Integer): LabeledPolynomial;
  lp_LabeledPolynomial2(
      poly_crtCoeffs(lp_eval(polyaccum), modulo, lp_eval(polycomp), prime),
      lp_sgn(polyaccum)
  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparators for sorting LabeledPolynomials

% Returns t if lead(eval(lp1)) < lead(eval(lp2))
asserted procedure lp_cmpLPLead(lp1: LabeledPolynomial,
                                  lp2: LabeledPolynomial): Boolean;
  poly_cmpPolyLead(lp_eval(lp1), lp_eval(lp2));

% Returns t if lead(eval(lp2)) < lead(eval(lp1))
asserted procedure lp_cmpLPLeadReverse(lp1: LabeledPolynomial,
                                        lp2: LabeledPolynomial): Boolean;
  lp_cmpLPLead(lp2, lp1);

% Returns t if total_degree(lead(eval(lp1))) < total_degree(lead(eval(lp2)))
asserted procedure lp_leadTotalDegreeCmp(lp1: LabeledPolynomial,
                                          lp2: LabeledPolynomial): Boolean;
  poly_leadTotalDegreeCmp(lp_eval(lp1), lp_eval(lp2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endmodule;  % end of module f5lp

end;  % of file
