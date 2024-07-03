function [alpha, alphas] = multi_index(nd, p)
%multi_index
%
% Tiangang Cui, 29/Oct/2016

alphas      = cell(p+1,1);    % multi-index
alphas{1}   = zeros(1,nd);    % multi-index for length 0

if nd == 1  % dimension = 1
    for q = 1:p
        alphas{q+1} = q;
    end
else
    for q = 1:p
        s           = nchoosek(1:nd+q-1,nd-1);
        s1          = zeros(size(s,1),1);
        s2          = (nd+q)+s1;
        alphas{q+1} = flipud(diff([s1 s s2],1,2))-1;   % -1 due to MATLAB indexing
        if sum(alphas{q+1},2) ~= q*ones(nchoosek(nd+q-1,nd-1),1)
            error('The sum of each row has to be equal to q-th order');
        end
    end
end
alpha = cell2mat(alphas);
