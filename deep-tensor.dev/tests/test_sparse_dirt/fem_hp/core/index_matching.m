function ind_matching = index_matching(data1, data2, tol)
%index_matching
%
% find corresponding index in data1 for each element of data2
%
% Tiangang Cui, 17/Oct/2016

n1 = size(data1, 1);
n2 = size(data2, 1);

ind_matching = zeros(n2,1);
for i = 1:n2
    dist = sum((repmat(data2(i,:), n1, 1) - data1).^2, 2).^0.5;
    ind = dist < tol;
    if sum(ind) > 1
        disp('Error: multiple matching')
    end
    if sum(ind) > 0
        ind_matching(i) = find(ind);
    end
end


end

