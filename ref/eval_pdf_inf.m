function value = eval_pdf_inf(sirt, x)
%{
The output should be same as that from eval_pdf function
The only modification here is to force the output which has input outside 
the range to be 0
This requires one step of check before passing inputs to eval_pdf
%}


flag = any(isnan(x)); % range flag, telling which column is in the range
ind = find(~flag);
value = zeros(1, size(x, 2));
if ~isempty(ind)
    value(ind) = eval_pdf(sirt, x(:, ind));
end
end


