function [x,hist,hist_err, iters] = Newton(f,x_0,tol,maxit)
% Input:
%   f : a user supplied function
%   x_0 : initial guess
%   tol : a positive real number (the stopping tolerance)
%   maxit: a positive integer specifying the max number
%   of iterations allowed.
% Output:
%   x : approximate solution to f(x) = 0
%   hist : an array (a vector) of the values of x_k
%   hist_err : an array (a vector) of the error, i.e., x^* - x_k
%   iter : the number of iterations taken

iters = 1;                               

% Initialize the vector holding the value of x.
hist(1) = x_0;                             
x = x_0;

% Compute value of function at x_0
[fval, df] = feval(f,x);

% Initialize the vector holding the error.
hist_err(1) = fval;                              

while (abs(fval)>tol && (iters<maxit))
    
    % Update the solution
    x = x - fval/df;
    
    % Compute value of function at current point for stopping criteria
    [fval, df] = f(x);
    
    % Save the history
    hist(iters+1) = x;
    hist_err(iters+1) = fval;
    
    % Increase the iteration number
    iters = iters+1;
    
end

iters = iters -1;


end