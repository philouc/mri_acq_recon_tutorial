% function plot_colored_curve(s,col_values,bound,n_colors) 
%
% This function plots a curve with color depending on the value of a given
%vector (e.g., the l2 norm of the speed at a given time)
%
% Input: s: \in R^{nd} a n x d array representing the curve
%        color values: \in R^n list of norm values
%        bound: maximum theoretical value of the criterion
%        n_colors: number of colors
%
% Output: the curve with colored varying with respect to the the chosen
% criterion

function plot_colored_curve(s,col_values,bound,n_colors)
    if nargin<5
        n_colors=100;
    end
    H=jet(n_colors);
    color_values_rescaled=min(n_colors,floor(1+(n_colors-1)*(col_values/bound)));
    scatter(s(:,1),s(:,2),12,H(color_values_rescaled,:),'filled')
    colorbar('YTick',[0 1],'YTickLabel', {'0';'Bound'}, 'FontSize',20)
