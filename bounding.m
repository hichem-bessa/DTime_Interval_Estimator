function [y_, y_bar] = bounding(x)
%BOUNDING Decomposes a signal into its positive and residual components.
%
% [y_, y_bar] = bounding(x)
%
% Given a vector or matrix `x`, this function computes:
%   y_bar = max(x, 0)     : element-wise positive part of x
%   y_    = y_bar - x     : residual (negative or zero components flipped)
%
% This is commonly used in interval observer implementations for
% bounding uncertainties or signal deviations.
%
% Example:
%   x = [-2; 3];
%   [y_, y_bar] = bounding(x)
%   => y_bar = [0; 3], y_ = [2; 0]

y_bar = max(x, 0);
y_ = y_bar - x;

end
