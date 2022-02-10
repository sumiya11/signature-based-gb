lex: hm
off1 'allfac;

put('initialize, 'psopfn, 'ht_initialize);
put('insert, 'psopfn, 'ht_insert);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computes the hash of x
symbolic procedure hash(x);
    begin integer h;
        if pairp x then
            h := length(x)     % if list then length
        else if fixp x then
            h := x             % if integer then identity
        else
            hashing_error();   % otherwise not supported
        return h
    end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize hashtable, u has only one parameter - initial table size
symbolic procedure ht_initialize(u);
    begin integer n;
        print(u);
    % here u is a list with 1 number
        if null u then
            input_error();

        n := reval pop u;

        % hashtable internally would be a list (allocated, load, values, table)
        hashtable := ht_initialize1(n);

        % convert to algebraic list
        hashtable := 'list . hashtable;

        return hashtable
    end;

symbolic procedure ht_initialize1(n);
    begin scalar values, table;
    integer load, allocated;
        values := mkvect(n);  % to store objects
        table  := mkvect(n);
        allocated := n;
        return {allocated, load, values, table}
    end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symbolic procedure ht_insert(u);
    begin scalar hashtable, newval, newhash, sym;

        if length u neq 2 then
            input_error();

        hashtable := reval pop u;
        newval    := reval pop u;

        print("In insert");
        print(hashtable);
        print(newval);

        sym := cdr hashtable;

        newhash   := ht_insert1(sym, newval);

        cdr hashtable := sym;

        return newhash
    end;

symbolic procedure ht_insert1(hashtable, newval);
    begin scalar h, allocated, values, table;
        allocated   := car hashtable;
        valuestable := cddr hashtable;

        print("internal");
        print(valuestable);

        h := hash(newval);
        h := mod(h, allocated);

        putv(cadr valuestable, h, newval);
        putv(car valuestable, h, h);

        %cadddr hashtable := table;
        %caddr hashtable  := values;

        return h
    end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symbolic procedure input_error();
    begin
        rederr "bad input";
    end;

symbolic procedure hashing_error();
    begin
        rederr "hashing for this type is not supported, sorry";
    end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


endmodule;

% Alex: hm?
tr ht_initialize;
tr ht_initialize1;
tr ht_insert;
tr ht_insert1;

end;


