dir := lisp lastcar lto_stringsplit(get!-current!-directory(), {'!/});
num1 := lisp caddr assoc ('(core_normalform), profile_alist!*);
num2 := lisp caddr assoc ('(core_topReductionF5), profile_alist!*);
num3 := lisp caddr assoc ('(core_interreduceInput), profile_alist!*);
total := lisp (time() - profile_time!*);
share dir, num1, num2, num3, total;
str := lisp lto_sconcat({"echo '", dir, ",", num1, ",", num2, ",", num3, ",",total , "' >> ~/signature-based-gb/f5/putin2.csv"});
share str;
lisp system str;
end;