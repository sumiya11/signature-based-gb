
module lists;

put('ReverseList, 'psopfn, 'lists_reverse);
put('MergeLists, 'psopfn, 'lists_merge);


asserted procedure ReverseError(): Void;
    rederr "bad arguments for reverse";

asserted procedure MergeError(): Void;
    rederr "bad arguments for merge";


% reverse u (almost) inplace
asserted procedure lists_reverse(u: List): List;
    begin scalar input, current, next, result, output;
        % input checking
        if null u then
            ReverseError();
        if length u neq 1 then
            ReverseError();

        % convert to symbolic
        input := cdr reval pop u;

        % reverse
        result  := nil;
        current := input;
        next    := cdr current;
        while next do
            begin;
                next := cdr current;
                print {"next, current, result", next, current, result};
                cdr current := result;
                result  := current;
                current := next;
            end;

        % convert to export
        output := 'list . result;

        return output
    end;

% merge two lists from u
asserted procedure lists_merge(u: List): List;
    begin scalar left, right, result, output;
        if null u then
            MergeError();
        if length u neq 2 then
            MergeError();

        left  := cdr reval pop u;
        right := cdr reval pop u;

        result := nil;
        while left and right do
            begin scalar iteml, itemr;
                iteml := car left;
                itemr := car right;
                print {"iteml, itemr, result", iteml, itemr, result};
                if iteml < itemr
                    then <<result := nconc(result, {iteml}); left  := cdr left>>
                    else <<result := nconc(result, {itemr}); right := cdr right>>;
        end;

        if left  then result := nconc(result, left);
        if right then result := nconc(result, right);

        output := 'list . result;

        return output
    end;

endmodule;

tr lists_reverse;
tr lists_merge;

end;  % of file

