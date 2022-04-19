
module f5hash;

% struct MonomHashtable
%     {'hashtable, table, values, filled}
%
%   values - an array, maps INDEX --> KEY . EXPONENT
%   table  - a htable, maps KEY --> HASH --> INDEX
%
%   An EXPONENT is mapped to KEY first and stored in values at INDEX.
%   Then KEY->INDEX is inserted into table.
%   Finally, a single INDEX is returned to user to represent an EXPONENT.
%

struct MonomHashtable;

asserted procedure hash_MonomHashtable();
  begin scalar ht, values;
        integer filled;

    table  := mkhash(8, 'equal, nil);
    values := mkvect(8);
    filled := 0;

    return {'hashtable, table, values, filled}
  end;

asserted inline procedure hash_getTable(ht);
  cadr ht;

asserted inline procedure hash_getValues(ht);
  caddr ht;

asserted inline procedure hash_getFilled(ht);
  cadddr ht;

asserted inline procedure hash_incFilled(ht);
  cadddr ht := cadddr ht #+ 1;

asserted inline procedure hash_getExpByIndex(ht, i);
  cadr getv(hash_getValues(ht), i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for EXPONENT ex return a KEY number (hash)
asserted procedure hash_hashExp(ht, ex: List);
  for each i in ex sum(i);

% Insert EXPONENT vector `ex` into hashtable `ht`
% and return a unique integer identifier INDEX
%
asserted procedure hash_insertExp(ht, ex: List);
  begin scalar table, values;
        integer hashE, vidx;
    table  := hash_getTable(ht);
    values := hash_getValues(ht);

    % calculate KEY (hash) of new exponent
    hashE := hash_hashExp(ht, ex);

    % vidx points to values array, or nil otherwise
    vidx  := gethash(hashE, table);

    % TODO: this is not really correct in terms of uniqueness

    % if vidx is nil, then exponent ex is not in hashtable yet
    % Add vidx to hashtable
    if null vidx then <<
      vidx := hash_getFilled(ht);
      hash_incFilled(ht);

      % table:  KEY --> INDEX 1..n
      puthash(hashE, table, vidx);

      % values: INDEX 1..n --> (KEY . EXPONENT)
      putv(values, vidx, {hashE, ex})
    >>;

    return vidx
  end;

% Check if monom i1 divides monom i2
% assuming exponent vectors are stored in hashtable `ht`
asserted procedure hash_isMonomDiv(ht, i1: Integer, i2: Integer);
  begin scalar e1, e2;
    e1 := hash_getExpByIndex(ht, i1);
    e2 := hash_getExpByIndex(ht, i2);

    return hash_isExpDiv(e1, e2)
  end;

% return the product of monomials i1 and i2
asserted procedure hash_monomProd(ht, i1: Integer, i2: Integer);
  begin scalar e1, e2, es;
    e1 := hash_getExpByIndex(ht, i1);
    e2 := hash_getExpByIndex(ht, i2);

    es := hash_sumExp(e1, e2);

    return hash_insertExp(ht, es)
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if exponent e1 divides e2
asserted procedure hash_isExpDiv(e1: List, e2: List);
	if null e1 then
		t
	else if car e1 > car e2 then
		nil
	else
		hash_isExpDiv(cdr e1, cdr e2);

% return sum of exponent vectors e1, e2
asserted procedure hash_sumExp(e1: List, e2: List): List;
	if null e1 then
		nil
	else
		(car e1 + car e2) . hash_sumExp(cdr e1, cdr e2);

endmodule;

% trst hash_insertExp;

end;  % of file
