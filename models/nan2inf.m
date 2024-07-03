function x = nan2inf(x)
% input is a 1 \times N vector
% output is a vector of the same size
% replace all nan in x by -inf
x(isnan(x)) = inf;
end