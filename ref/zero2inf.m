function y = zero2inf(x)
% used in dirt-ftt to eliminate the points out of bound
y = x;
ind = x == 0;
y(ind) = inf;
end


