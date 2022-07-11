module f5univpol;

% SparseVector - sparsely stored vector
% 
% SparseVector object is stored as
%  {'sv, indices, values} 
% where `indices` is a list of integers - nonzero vector entries indices, 
% and `values` is a list of corresponding nonzero coefficients.
%
% Values are guaranteed to be nonzero, indices are guaranteed to be sorted.
%
% Zero SparseVector is represented as {'sv, nil, nil}

% Creates sparse vector with the given indices and values
asserted procedure univ_SparseVector(indices, values): SparseVector;
   {'sv, indices, values};

% Return the vector indices
asserted procedure univ_getIndicesVector(vect: SparseVector): List;
   cadr vect; 

% Return the vector values
asserted procedure univ_getValuesVector(vect: SparseVector): List;
   caddr vect; 

% Sets the vector indices
asserted procedure univ_setIndicesVector(vect: SparseVector, newindices: List): List;
   cadr vect := newindices; 

% Sets the vector values
asserted procedure univ_setValuesVector(vect: SparseVector, newvals: List): List;
   caddr vect := newvals; 

% Return the pivot index - the first indices entry
asserted procedure univ_pivotIndexVector(vect: SparseVector): Integer;
   car univ_getIndicesVector(vect); 

% Return the pivot value - the first values entry
asserted procedure univ_pivotValueVector(vect: SparseVector): Coeff;
   car univ_getValuesVector(vect); 

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
   {'mm, 0, 0, 0, nil, nil, mkhash(2^4, 'eq, 2), mkvect(0)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert polynomial `left` to a dense row, with the column labels as in the left matrix part
%  & Convert polynomial `right` to a dense row,  -//-
asserted procedure univ_convertToDenseRows(matrix: MacaulayMatrix, left: Polynomial, right: Polynomial): DottedPair;
   begin scalar nrightcols, leftRow, rightRow;
      nrightcols := cadddr matrix;
      
      %              1  x  x^2 ...   x^n
      % rightRow is [0, 0, 0, ...,0, 1]
      rightRow := mkvect(nrightcols);
      for i := 0:nrightcols do
         putv(rightRow, i, 0);
      putv(rightRow, nrightcols, 1);
      cadddr matrix := nrightcols + 1;

      nleftcols := caddr matrix;
      term2index := caddr cddddr matrix;   % any alternative ?
      index2term := cadddr cddddr matrix;
      for term in poly_getTerms(left) do <<
         if not gethash(term, term2index) then <<
            puthash(term, term2index, nleftcols);
            if nleftcols >= length(index2term) then <<
               newindex2term := mkvect(length(index2term) * 2);
               for i := 0:nleftcols do
                  putv(newindex2term, i, getv(index2term, i));
               index2term := newindex2term
            >>;
            nleftcols := nleftcols + 1;
            putv(index2term, nleftcols, term);
         >>
      >>;
      caddr matrix := nleftcols;

      leftRow := mkvect(nleftcols);
      for i := 0:nleftcols do
         putv(leftRow, i, 0);
      terms   := poly_getTerms(left);
      coeffs  := poly_getCoeffs(left);
      while terms do <<
         term  := pop(terms);
         coeff := pop(coeffs);
         idx := gethash(term, term2index);
         putv(leftRow, idx, coeff)
      >>;

      return leftRow . rightRow
   end;

trst univ_convertToRowDense;

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
      reducerVals := univ_getValuesVector(leftReducer);
      leftMul := negsq(getv(left, car reducerVals));
      for each col in reducerCols do <<
         newcf := addsq(pop(reducerVals), multsq(getv(left, col), leftMul));
         putv(left, col, newcf);
      >>;
      reducerCols := univ_getIndicesVector(rightReducer);
      reducerVals := univ_getValuesVector(rightReducer);
      for each col in reducerCols do <<
         newcf := addsq(pop(reducerVals), multsq(getv(right, col), leftMul));  % not a mistake
         putv(right, col, newcf);
      >>
   end;

trst univ_reduceDenseWithSparse;

asserted procedure univ_reduceRows(matrix: MacaulayMatrix, left: Vector, right: Vector): Boolean;
   begin scalar leftRows, rightRows, lpivot, lrow, rrow, reduced;
      leftRows  := cadr cdddr matrix; 
      rightRows := caddr cdddr matrix;
      for each lrow in leftRows do <<
         lpivot := univ_pivotIndexVector(lrow);
         if not (getv(left, lpivot) = 0) then <<
            rrow := pop(rightRows);
            univ_reduceDenseWithSparse(left, right, lrow, rrow);
         >>
      >>;
      reduced := t;
      for i := 0:length(left) do <<
         if not (getv(left, i) = 0) then
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
      leftCoeffs  := univ_getValuesVector(left);
      rightCoeffs := univ_getValuesVector(right);
      leftMult := denr(car leftCoeffs) ./ numr(car leftCoeffs);
      leftCoeffs  := for each c in leftCoeffs collect 
                        multsq(leftMult, c);
      rightCoeffs := for each c in rightCoeffs collect 
                        multsq(leftMult, c);
      univ_setIndicesVector(left, leftCoeffs);
      univ_setIndicesVector(right, rightCoeffs)
   end;

trst univ_normalizeRows;

asserted procedure univ_addMatrixRow(matrix: MacaulayMatrix, left: SparseVector, right: SparseVector);
   begin scalar leftRows, rightRows;
      leftRows  := cadr cdddr matrix; 
      rightRows := caddr cdddr matrix;

      % of course, this should be done better; but not a problem for now
      leftRows := reversip(leftRows);
      rightRows := reversip(rightRows);

      push(leftRow, leftRows);
      push(rightRow, rightRows);

      leftRows := reversip(leftRows);
      rightRows := reversip(rightRows);

      cadr cdddr matrix := leftRows;
      caddr cdddr matrix := rightRows;

      return leftRow . rightRow
   end;

trst univ_addMatrixRow;

% Given a dense vector row, produce a SparseVector from it
% e.g., [1, 0, 0, 2] --> {'sv, {0, 3}, {1, 2}} 
asserted procedure univ_extractSparseRow(row: Vector): SparseVector;
   begin scalar cols, coeffs;
      for i := 0:length(row) do << 
         if not (getv(row, i) = 0) then <<
            push(i, cols);
            push(coffs, getv(row, i))
         >>
      >>;
      return univ_SparseVector(reversip(cols), reversip(coeffs))
   end;

trst univ_extractSparseRow;

asserted procedure univ_findLinearRelation(matrix: MacaulayMatrix, left: Polynomial, right: Polynomial): DottedPair;
   begin scalar leftRow, rightRow, reduced; 
      leftRow . rightRow := univ_convertToDenseRows(matrix, left, right);
      reduced := univ_reduceRows(matrix, leftRow, rightRow);
      leftRow  := univ_extractSparseRow(leftRow);
      rightRow := univ_extractSparseRow(rightRow);
      if not reduced then <<
         univ_normalizeRows(leftRow, rightRow);
         univ_addMatrixRow(matrix, leftRow, rightRow)
      >>;
      return reduced . rightRow
   end;

trst univ_findLinearRelation;

%  
asserted procedure univ_extractGenerator(relation: SparseVector, idx: Integer): Polynomial;
   begin scalar terms, coeffs;
      coeffs := univ_getValuesVector(relation);
      terms := for each x in univ_getIndicesVector(relation) 
               collect poly_ithVariable(idx, x);
      return poly_Polynomial(terms, coeffs)
   end;

% Computes the univariate polynomial of the ideal generated by G in index i.
% G - a list of Polynomials, a Groebner basis of zero dimensional ideal in n variables.
% i - an Integer, a variable index
asserted procedure univ_univpol1(G: List, i: Integer): Polynomial;
   begin scalar existsRelation, relation, tobereduced, reduced, matrix;
         integer n;
      matrix := univ_MacaulayMatrix();
      while not existsRelation do <<
         tobereduced := poly_Polynomial({poly_ithVariable(i, n)}, {poly_oneCoeff()});
         reduced := core_normalFormReducers(tobereduced, G, t);
         existsRelation . relation := univ_findLinearRelation(matrix, reduced, tobereduced); 
         if existsRelation then
            relation := univ_extractGenerator(relation, i);
         n := n + 1
      >>;
      return relation
   end;

trst univ_univpol1;

endmodule;  % end of module f5univpol

end;  % of file
