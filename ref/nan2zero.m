function A = nan2zero(A)
% input is an array. output should be an array of same size with NaN
% replaced by zeros.
A(isnan(A)) = 0;
A(isinf(A)) = 0;
end

