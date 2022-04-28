function [x, hist, hist_err, iters] = secant(f,x_0,x_1,tol,maxit)
% Input:
%   f : a user supplied function
%   x_0, x_1 : initial guess
%   tol : a positive real number (the stopping tolerance)
%   maxit: a positive integer specifying the max number
%   of iterations allowed.
% Output:
%   x : approximate solution to f(x) = 0
%   hist : an array (a vector) of the values of x_k
%   hist_err : an array (a vector) of the error, i.e., x^* - x_k
%   iter : the number of iterations taken

iters = 2;                               

% Initialize the vector holding the value of x.
hist(1) = x_0;  
hist(2) = x_1;
x = x_0;

% Compute value of function at x_0 and x_1
[fval0, ~] = feval(f,x_0);
[fval1, ~] = feval(f,x_1);

% Initialize the vector holding the error.
hist_err(1) = fval0;   
hist_err(2) = fval1; 

while (abs(fval1)>tol && (iters<maxit))
    
    % Update the solution
    x = x_1 - fval1*(x_1-x_0)/(fval1-fval0);
    
    % Compute value of function at current point for stopping criteria
    x_0 = x_1;  fval0 = fval1;
    x_1 = x;    [fval1, ~] = f(x_1);
    
    % Save the history
    hist(iters+1) = x;
    hist_err(iters+1) = fval1;
    
    % Increase the iteration number
    iters = iters+1;
    
end

iters = iters - 2;


end