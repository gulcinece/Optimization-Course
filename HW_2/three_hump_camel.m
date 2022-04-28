function [f, gradf, hessianf ] = three_hump_camel(x)

% Description:
% Calculates Three Hump Camel function, its gradient and Hessian matrix at
% given point x.
%
% Input:
%  x: given point
%
% Output:
%   f: Three Hump Camel function at point x
%   gradf: Gradient of Three Hump Camel fuction at point x
%   hessianf: Hessian matrix of Three Hump Camel function at point x
%
% Usage:
% three_hump_camel(x)

f        = 2*x(1)^2 - 1.05*x(1)^4 + (x(1)^6)/6 + x(1)*x(2) + x(2)^2;
gradf    = [4*x(1) - 4.2*x(1)^3 + x(1)^5 + x(2); x(1)+2*x(2)];
hessianf = [4 - 12.6*x(1)^2 + 5*x(1)^4, 1; 1, 2];
end