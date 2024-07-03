function [udata, ind_data, ind_udata] = robust_unique(data, tol)
%robust_unique
%
% find unique rows of a input matrix, 
% 
% inputs:
%   data: input data
%   tol:  tolorence that defines two values are the same
%
% outputs: 
%   udata:      unique rows
%   ind_data:   incides of data that correspond to udata
%   ind_udata:  incides of udata that correspond to data
%
% udata = data(ind_data,: ) 
% data = udata(ind_udata,:)
%
% Tiangang Cui, 17/Oct/2016

tmp = round(data/tol);
[~,ind_data,ind_udata] = unique(tmp, 'rows', 'stable');
udata = data(ind_data,:);


% compute the pairwise distance 

%{
n = size(data, 1);

ind_udata = zeros(n,1);
active = true(n, 1);
j = 0;
for i = 1:n
    if active(i)
        j = j+1;
        ii  = i:n;
        dist = sum((repmat(data(i,:), (n-i+1), 1) - data(ii,:)).^2, 2).^0.5;
        ind = ii(dist < tol);
        ind_udata(ind) = j;
        active(ind) = false;
    end
end

[~, ind_data, ~] = unique(ind_udata, 'stable');

%[~, ind_data, t_ind_udata] = unique(ind_udata, 'stable');
%sum(abs(t_ind_udata - ind_udata))

udata = data(ind_data, :);
%}

%{
tmp = round(data/tol);
[~,ind_data1,ind_udata1] = unique(tmp, 'rows', 'stable');
udata1 = data(ind_data1,:);

norm(udata-udata1)
norm(ind_data-ind_data1)
norm(ind_udata-ind_udata1)
%}

end
