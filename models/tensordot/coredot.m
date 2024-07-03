function y = coredot(core1, core2)
    % compute the dot product of two tensor cores
    % output is a matrix of size r1(core1)r1(core2) * r2(core1)r2(core2)
    [r1, n1, r2] = size(core1);
    [s1, n2, s2] = size(core2);
    if n1 ~= n2
        error("size of two tensor cores should be equivalent.")
    end
    c1 = reshape(permute(core1, [1,3,2]), r1*r2, n1, 1);
    c2 = reshape(permute(core2, [2,1,3]), 1, n2, s1*s2);
    y = sum(c1.*c2, 2);
    y = reshape(y, r1, r2, s1, s2);
    y = reshape(permute(y, [3,1,4,2]), r1*s1, r2*s2);
end