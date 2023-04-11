function value = eval_rt_check(sirt, x, poly)
%{
The output should be same as that from eval_pdf function
The only modification here is to force the output which has input outside 
the range to be 0
This requires one step of check before passing inputs to eval_pdf
%}

L = poly.domain(1);
U = poly.domain(2);
rg_flag = all((x <= U & x >= L), 1); % range flag, telling which column is in the range
ind = find(rg_flag);
value = NaN(size(x, 1), size(x, 2));
if ~isempty(ind)
    value(:, ind) = eval_rt(sirt, x(:, ind));
end
end


