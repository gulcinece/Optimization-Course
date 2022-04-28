% CA2_Exercise1.m
%
% Author: Pelin Ciloglu
%
% Description:
% Primal-dual interior point method is applied for the different initial
% vectors and sigma values. Also, MATLAB built-in function 'linprog' is
% tested for the same parameters.
%
% The input:
%   A: This is the Ax=b matrix. 
%   b: Vector. This is the right hand side of Ax=b.
%   c: Vector. This is from minimize  J(x) = c'x. 
%   X0: initial vector
%   S0: initial vector
%   Y0: initial vector
%   Nmax: maximum number of iteration
%   p: a parameter which is greater than n (p>n)
%   sigma: real number in [0,1]
%
% Ouput:
%   criteria = c'*X
%   X:optimal value
%   Y:optimal value
%   S:optimal value
%
% Usage:
% [criteria, X, Y, S] = Interior_Point(A,b,C,X0,Y0,S0,Nmax,p,sigma)

%% Part a

% Problem 1
A = [1, 1, 1, 0; 2, 1, 0, 1];
C = [-1; 1; 0; 0];
b = [40; 60];

X0 = [0; 30; 10; 30];
% X0 = [15; 15; 5; 20];
% X0 = [30; 30; 5; 20];

Y0=[-1; -1];
S0 = [2; 3; 1; 1];

Nmax = 500;
tol = 1e-6;
[criteria, iter, X, Y, S] = Interior_Point(A,b,C,X0,Y0,S0,Nmax,tol)


%% Part c

format short;
options = optimset('LargeScale','off');
[X_lin,FVAL,EXITFLAG,OUTPUT]=linprog(C,[],[],A,b,zeros(size(C)),[],[],options)





