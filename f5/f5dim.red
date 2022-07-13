module f5dim;

asserted procedure dim_maxIndependent1(HT: List, k: Integer, U: List, M: List): List;
   begin scalar x;
      newM := M;
      for i := k:n do <<
         if U + {Xi} does not intersect with HT then
            newM = dim_maxIndependent1(HT, i + 1, U + {Xi}, newM)
      >>;
      if U is not contained in any V from newM then
         newM = newM + {U}
      return newM
   end;

asserted procedure dim_maxIndependent1(G: List): List;
   begin scalar x;
      headTerms := for each f in G collect poly_leadTerm(f);
      M := dim_dimensionRecursive(headTerms, 1, nil, nil);
      d := 0;
      U := nil;
      for each independentSet in M do <<
         if length(independentSet) > d then <<
            d := length(independentSet);
            U := independentSet
         >>
      >>;
      return U
   end;

endmodule;  % end of module f5dim

end;  % of file
