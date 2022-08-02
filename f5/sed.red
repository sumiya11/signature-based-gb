dir := lisp lastcar lto_stringsplit(get!-current!-directory(), {'!/});
num1 := lisp caddr assoc ('(core_normalform), profile_alist!*);
num2 := lisp caddr assoc ('(core_topReductionF5), profile_alist!*);
total := lisp (time() - profile_time!*);
share dir, num1, num2, total;
str := lisp lto_sconcat({"echo '", dir, ",", num1, ",", num2, ",",total , "' >> ~/signature-based-gb/f5/putin.csv"});
share str;
lisp system str;
end;