function v = ttdot(cores1, cores2)
    % compute the inner product of two tensor trains represented by cores
    % output is a scalar

    % check the dimension and size of inputs
    if length(cores1) ~= length(cores2)
        error("Two input tensor trains should have same length.")
    end
    d = length(cores1);
    % s1 = zeros(1, d);
    % s2 = zeros(1, d);
    % for k = 1:d
    %     s1(k) = size(cores1{k}, 2);
    %     s2(k) = size(cores2{k}, 2);
    % end
    % if any(s1~=s2)
    %     error("Two input tensor trains should have same size.")
    % end

    v = coredot(cores1{1}, cores2{1});
    for k = 2:d
        v = v*coredot(cores1{k}, cores2{k});
    end
end