function mllkd = minus_log_like_dirt(obs, d)
%
% minus log likelihood for a given observable model output 
% gradient w.r.t. observable model output
% Fisher information matrix at observable model output

switch obs.like
    case {'normal'}
        misfit = (d - obs.data)./obs.std;
        mllkd  = 0.5*misfit(:).^2; % minus log-likelihood
    case {'poisson'}
        mllkd = - obs.data(:).*log(d(:)) + d(:); 
    otherwise
        error('likelihood not implemented')
end

end