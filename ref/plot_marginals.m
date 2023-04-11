function plot_marginals(data, varargin)
% This function plots the 1d marginal and 2d marginal densities estimated 
% from the samples. Run the function compute_marginals before calling this
% function.
%
% Inputs: 
%
% data:     data used for plotting
%
% Optional parameters:
%
% box:      a bounding box in the form of [left, bottom, width, height] for 
%           placing the plots. The default is [0.08,0.06,0.9,0.9]
%
% Example:
% 
% load samples.mat
% data = compute_marginals(samples, 'ngrid', 100, 'qlist', [0.9,0.7,0.5,0.3,0.1])
% plot_marginals(data)
%
% Tiangang Cui, August, 2019

defaultBox = [0.08, 0.06, 0.9, 0.9];

p = inputParser;
%
addRequired(p,'data');
addOptional(p,'box', defaultBox, @(x) isnumeric(x) && length(x) == 4);
parse(p,data,varargin{:});

box = p.Results.box;
d   = data.d;

% parameters for the subplots
% first create the bounding box [left, bottom, width, height]
left    = box(1);
bottom  = box(2);
width   = box(3);
height  = box(4);
size_x  = width/d;
size_y  = height/d;
fac     = 0.98;

figure('position', [100, 100, 800, 800])
for i = 1:d
    % place the one d marginal plot alond the diagonal
    
    subplot('position',[left+size_x*(i-1) bottom+size_y*(d-i) size_x*fac size_y*fac])
    plot(data.xi{i},data.fi{i});
    
    if i ~= 1
        set(gca,'YTick', [])
    end
    if i ~= d
        set(gca,'XTick', [])
    else
%         xlabel(['x_{' num2str(d) '}'])
    end
    
    % then create two d marginals
    for j = (i+1):d
        subplot('position',[left+size_x*(i-1) bottom+size_y*(d-j) size_x*fac size_y*fac])
        contour(data.xi{i}, data.xi{j}, data.fs{i,j}, data.qants{i,j})
        if j ~= d
            set(gca,'XTick', [])
        else
%             xlabel(['x_{' num2str(i) '}'])
        end
        if i ~= 1
            set(gca,'YTick', [])
        else
%             ylabel(['x_{' num2str(j) '}'])
        end
    end
%     set(gca,'XTicklabel',[], 'YTicklabel', [])
end

end
