%==========================================================================
%
% derivative  Numerical differentiation of data (i.e. arrays) over the 
% domain of the data or at specified points.
%
%   dy = derivative(x,y)
%   dy = derivative(x,y,x_star)
%
% See also diff, gradient, iderivative, derivest.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-08-27
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Basic_Numerical_Calculus.pdf
%
% REFERENCES:
%   [1] http://www.ohiouniversityfaculty.com/youngt/IntNumMeth/lecture27.pdf
%   [2] https://en.wikipedia.org/wiki/Finite_difference_method
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - ((N+1)×1 or 1×(N+1) double) independent variable data
%   y       - ((N+1)×1 or 1×(N+1) double) dependent variable data
%   x_star  - (OPTIONAL) (p×1 or 1×p double) points at which to 
%             differentiate
%
% -------
% OUTPUT:
% -------
%   dy      - ((N+1)×1, 1×(N+1), p×1, or 1×p double) derivative of y = f(x)
%             w.r.t. x evaluated at:
%               --> all the points in x (cumulative differentiation)
%               --> all the points in x_star (point differentiation)
%
% -----
% NOTE:
% -----
%   --> If "x_star" is not input, then "dy" stores the derivative of 
%       y = f(x) with with respect to x at the points in "x".
%   --> N = number of data points (i.e. length of "y" and "x")
%   --> p = number of points to differentiate at (i.e. length of "x_star")
%
% THIRD-PARTY FILE — not included in this repository.
% Download from: https://github.com/tamaskis/numerical_differentiation-MATLAB
% (also on MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/97267)
%
%==========================================================================
function [dy,x] = derivative(x,y,x_star)
    error('derivative.m is a third-party file. Download from: https://github.com/tamaskis/numerical_differentiation-MATLAB');
end
