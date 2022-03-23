module f5poly;

load!-package 'dp;

off1 'allfac;

% return exponent vector of zeros
asserted procedure poly_zeroExponent(): List;
  ev_zero();

% return sum of exponent vectors e1, e2
asserted procedure poly_sumExponents(e1: List, e2: List): List;
  ev_sum(e1, e2);

% return difference of exponent vectors e1, e2
asserted procedure poly_difExponents(e1: List, e2: List): List;
  ev_dif(e1, e2);

% return lcm of exponent vectors e1, e2
asserted procedure poly_lcmExponents(e1: List, e2: List): List;
  ev_lcm(e1, e2);

% comparator for exponent vectors e1, e2
asserted procedure poly_cmpExponents(e1: List, e2: List): List;
  ev_comp(e1, e2);

% check if exponent e1 divides e2
asserted procedure poly_dividesExponents(e1: List, e2: List): List;
  ev_divides!?(e1, e2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% returns leading exponent of poly
asserted procedure poly_leadExponent(poly: Polynomial): List;
  dip_evlmon(poly);

% returns leading coefficient of poly
asserted procedure poly_leadCoeff(poly: Polynomial); % so, what does this return exactly
  dip_lbc(poly);

% returns length of poly, i.e., the number of terms
asserted procedure poly_length(poly: Polynomial): Integer;
  dip_length(poly);


endmodule;


end;  % of file
