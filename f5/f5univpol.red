module f5univpol;

% SparseVector - sparsely stored vector
% 
% SparseVector object is stored as
%  {'sv, indices, coeffs} 
% where `indices` is a list of integers - nonzero vector entries indices, 
% and `coeffs` is a list of corresponding nonzero coefficients.
%
% Coeffs are guaranteed to be nonzero, indices are guaranteed to be sorted.
%
% Zero SparseVector is represented as {'sv, nil, nil}

% Creates sparse vector with the given indices and coeffs
asserted procedure univ_SparseVector(indices: List, coeffs: List): SparseVector;
   {'sv, indices, coeffs};

% Returns the vector indices
asserted procedure univ_getIndicesVector(vect: SparseVector): List;
   cadr vect; 

% Returns the vector coeffs
asserted procedure univ_getCoeffsVector(vect: SparseVector): List;
   caddr vect; 

% Sets the vector indices
asserted procedure univ_setIndicesVector(vect: SparseVector, newindices: List): List;
   cadr vect := newindices; 

% Sets the vector coeffs
asserted procedure univ_setCoeffsVector(vect: SparseVector, newcoeffs: List): List;
   caddr vect := newcoeffs; 

% Returns the pivot index - the first indices entry
asserted procedure univ_pivotIndexVector(vect: SparseVector): Integer;
   car univ_getIndicesVector(vect); 

% Returns the pivot coeff - the first coeffs entry
asserted procedure univ_pivotCoeffVector(vect: SparseVector): Coeff;
   car univ_getCoeffsVector(vect); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MacaulayMatrix - matrix of shape
%  
%  | left | right |
%  | rows | rows  |
%
% where each column is labeled with a polynomial term; 
% the number of rows in both part is the same
% We perform elementary row operations on the _left part_ of the matrix,
% and apply the same operations to its _right part_.   
% 
% The matrix object is represented as 
%  {'mm, nrows, nleftcols, nrightcols, leftRows, rightRows, term2index, index2term}
% where
%  . nrows - number of rows in the matrix (same for the left and right parts),
%  . nleftcols / nrightcols - number of columns in the left / right part,
%  . leftRows / rightRows - list of left / right rows, each entry is a SparseVector object 
%  . term2index - hashtable that maps term to its index amidst the column labels in the _left part_
%  . index2term - vector that maps column label index to the corresponding term in the _left part_ 
%        x,    <-> 1
%        x^2   <-> 2
%              ...
%        xy^2  <-> 12  (indices are assigned arbitrarily, not in any particular term order)
%
% Important condition on leftRows: 
% each i-th row in leftRows is reduced w.r.t. {1..i-1} rows from leftRows.
% So, the left part of the matrix is in some sort of a normal reduced form 
% (except that left rows are not sorted by the pivot index, in general)

% Creates an empty `MacaulayMatrix`
asserted procedure univ_MacaulayMatrix(): MacaulayMatrix;
   {'mm, 0, 0, 0, nil, nil, mkhash(2^4, 'equal, 2), mkvect(0)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert polynomial `left` to a dense row, with the column labels as in the left matrix part
%  & Convert polynomial `right` to a dense row,  -//-
asserted procedure univ_convertToDenseRows(mmatrix: MacaulayMatrix, 
                                    left: Polynomial, right: Polynomial): DottedPair;
   begin scalar nrightcols, leftRow, rightRow;
      nrightcols := cadddr mmatrix;
      
      %              1  x  x^2 ...   x^n
      % rightRow is [0, 0, 0, ...,0, 1]
      rightRow := mkvect(nrightcols);
      for i := 0:nrightcols do
         putv(rightRow, i, poly_zeroCoeff());
      putv(rightRow, nrightcols, poly_oneCoeff());
      cadddr mmatrix := nrightcols + 1;

      nleftcols := caddr mmatrix;
      term2index := caddr cddddr mmatrix;   % any alternative ?
      index2term := cadddr cddddr mmatrix;
      for each term in poly_getTerms(left) do <<
         if not gethash(term, term2index) then <<
            puthash(term, term2index, nleftcols);
            if nleftcols >= length(index2term) then <<
               newindex2term := mkvect(length(index2term) * 2);
               for i := 0:nleftcols-1 do
                  putv(newindex2term, i, getv(index2term, i));
                cadddr cddddr mmatrix := newindex2term
            >>;
            index2term := cadddr cddddr mmatrix;
            putv(index2term, nleftcols, term);
            nleftcols := nleftcols + 1
         >>
      >>;
      caddr mmatrix := nleftcols;

      leftRow := mkvect(nleftcols);
      for i := 0:nleftcols do
         putv(leftRow, i, poly_zeroCoeff());
      terms   := poly_getTerms(left);
      coeffs  := poly_getCoeffs(left);
      while terms do <<
         term  := pop(terms);
         coeff := pop(coeffs);
         idx := gethash(term, term2index);
         prin2t {"HASH", idx};
         putv(leftRow, idx, coeff)
      >>;

      return leftRow . rightRow
   end;

trst univ_convertToDenseRows;

% Reduce `left` with `leftReducer` and `right` with `rightReducer` inplace.
%  left  = left  - C*leftReducer
%  right = right - C*rightReducer
% where `C` is chosen so that one (or more) coefficients in `left` is cancelled
% 
% ASSUMING values are SQ
asserted procedure univ_reduceDenseWithSparse(left: Vector, right: Vector, 
                                                leftReducer: SparseVector, rightReducer: SparseVector);
   begin scalar reducerCols, reducerVals, leftMul, col, newcf;
      reducerCols := univ_getIndicesVector(leftReducer);
      reducerVals := univ_getCoeffsVector(leftReducer);
      leftMul := negsq(getv(left, car reducerCols));
      for each col in reducerCols do <<
         newcf := addsq(multsq(pop(reducerVals), leftMul), getv(left, col));
         putv(left, col, newcf)
      >>;
      reducerCols := univ_getIndicesVector(rightReducer);
      reducerVals := univ_getCoeffsVector(rightReducer);
      for each col in reducerCols do <<
         newcf := addsq(multsq(pop(reducerVals), leftMul), getv(right, col));
         putv(right, col, newcf)
      >>
   end;

trst univ_reduceDenseWithSparse;

asserted procedure univ_reduceRows(mmatrix: MacaulayMatrix, left: Vector, right: Vector): Boolean;
   begin scalar leftRows, rightRows, lpivot, lrow, rrow, reduced;
      leftRows  := cadr cdddr mmatrix; 
      rightRows := caddr cdddr mmatrix;
      for each lrow in leftRows do <<
         lpivot := univ_pivotIndexVector(lrow);
         if not poly_iszeroCoeff!?(getv(left, lpivot)) then <<
            rrow := pop(rightRows);
            univ_reduceDenseWithSparse(left, right, lrow, rrow)
         >>
      >>;
      reduced := t;
      for i := 0:length(left)-1 do
         if not poly_iszeroCoeff!?(getv(left, i)) then
            reduced := nil;
      return reduced
   end;

trst univ_reduceRows;

% Divides coefficients in `left` and `right` by `C`
%  right = right / C 
%  left  = left  / C
% where `C` is the leading coefficient in `left`,
% mutating `left` and `right` inplace
asserted procedure univ_normalizeRows(left: SparseVector, right: SparseVector);
   begin scalar leftCoeffs, rightCoeffs;
      leftCoeffs  := univ_getCoeffsVector(left);
      rightCoeffs := univ_getCoeffsVector(right);
      leftMult := denr(car leftCoeffs) ./ numr(car leftCoeffs);
      leftCoeffs  := for each c in leftCoeffs collect 
                        multsq(leftMult, c);
      rightCoeffs := for each c in rightCoeffs collect 
                        multsq(leftMult, c);
      univ_setCoeffsVector(left, leftCoeffs);
      univ_setCoeffsVector(right, rightCoeffs)
   end;

trst univ_normalizeRows;

asserted procedure univ_addMatrixRow(mmatrix: MacaulayMatrix, left: SparseVector, right: SparseVector);
   begin scalar leftRows, rightRows;
      cadr mmatrix := cadr mmatrix + 1;

      leftRows  := cadr cdddr mmatrix; 
      rightRows := caddr cdddr mmatrix;

      % of course, this should be done better; but not a problem for now
      leftRows := reversip(leftRows);
      rightRows := reversip(rightRows);

      push(leftRow, leftRows);
      push(rightRow, rightRows);

      leftRows := reversip(leftRows);
      rightRows := reversip(rightRows);

      cadr cdddr mmatrix := leftRows;
      caddr cdddr mmatrix := rightRows;

      return leftRow . rightRow
   end;

trst univ_addMatrixRow;

% Given a dense vector row, produce a SparseVector from it
% e.g., [1, 0, 0, 2] --> {'sv, {0, 3}, {1, 2}} 
asserted procedure univ_extractSparseRow(row: Vector): SparseVector;
   begin scalar cols, coeffs;
      for i := 0:length(row)-1 do << 
         if not poly_iszeroCoeff!?(getv(row, i)) then <<
            push(i, cols);
            push(getv(row, i), coeffs)
         >>
      >>;
      return univ_SparseVector(reversip(cols), reversip(coeffs))
   end;

trst univ_extractSparseRow;

asserted procedure univ_findLinearRelation(mmatrix: MacaulayMatrix, left: Polynomial, right: Polynomial): DottedPair;
   begin scalar leftRow, rightRow, reduced; 
      leftRow . rightRow := univ_convertToDenseRows(mmatrix, left, right);
      reduced := univ_reduceRows(mmatrix, leftRow, rightRow);
      prin2t {"reduced", reduced, leftRow, rightRow};
      leftRow  := univ_extractSparseRow(leftRow);
      rightRow := univ_extractSparseRow(rightRow);
      if not reduced then <<
         univ_normalizeRows(leftRow, rightRow);
         univ_addMatrixRow(mmatrix, leftRow, rightRow)
      >>;
      return reduced . rightRow
   end;

trst univ_findLinearRelation;

%  
asserted procedure univ_extractGenerator(mmatrix: MacaulayMatrix, relation: SparseVector, 
                                             idx: Integer): Polynomial;
   begin scalar terms, coeffs;
      coeffs := univ_getCoeffsVector(relation);
      terms := for each x in univ_getIndicesVector(relation) 
                  collect poly_ithVariable(idx, x);
      return poly_Polynomial(reversip(terms), reversip(coeffs))
   end;

trst univ_extractGenerator;

% Computes the univariate polynomial of the ideal generated by G in index i.
% G - a list of Polynomials, a Groebner basis of zero dimensional ideal in n variables.
% i - an Integer, a variable index
asserted procedure univ_univpol1(G: List, i: Integer): Polynomial;
   begin scalar existsRelation, relation, tobereduced, reduced, mmatrix, flag;
         integer n;
      mmatrix := univ_MacaulayMatrix();
      prin2t {"MMatrix ", mmatrix};
      while not existsRelation do <<
         prin2t {"iteration ", n, poly_ithVariable(i, n)};
         tobereduced := poly_Polynomial({poly_ithVariable(i, n)}, {poly_oneCoeff()});
         flag . reduced := core_normalFormReducers(tobereduced, G, t);
         prin2t {tobereduced, " --> ", reduced};
         existsRelation . relation := univ_findLinearRelation(mmatrix, reduced, tobereduced); 
         if existsRelation then
            relation := univ_extractGenerator(mmatrix, relation, i);
         n := n + 1
      >>;
      return relation
   end;

trst univ_univpol1;

asserted procedure squarefree(f: Polynomial): Polynomial;
   begin scalar f;
      f := 1;
      for each fact in factorize(poly_2a(f)) do
         f := f * fact;
      return poly_f2poly(f)
   end;


% Computes the univariate polynomial of the ideal generated by G in index i.
% G - a list of Polynomials, a Groebner basis of zero dimensional ideal in n variables.
% i - an Integer, a variable index
asserted procedure f5_zerodimradical1(G: List): List;
   begin scalar fi, rad;
         integer n;
      for i := 1:length(global!-dipvars!*)-1 do <<
         fi := univ_univpol1(G, i);
         push(squarefree(fi), rad)
      >>;
      return rad
   end;

endmodule;  % end of module f5univpol

end;  % of file
