

function [X,NormDf,it] = SR1_inverse(fhandle,x0,tol,H0,maxit,alpha0,c,beta,amax)

% Description:
%  Obtains the iterations of a given function with given datas by using 
%  BFGS method 
%
% Input:
%  x0: initial point
%  tol: tolerance
%  gradf: gradient of the function
%  H0 : initial inverse Hessian matrix
%  maxit: maximum number of iterations
%  alpha0: initial step size
%  c: parameter
%  beta: parameter
%  amax: maximum number of iterations of Armijo function
%
% Output:
%  X: Iterations of fhandle up to convergence 
%  NormDf: Norm of Gradients of f(x)
%  it: Number of iterations
%
% Usage:
%  SR1_inverse(x0,tol,gradf,B0,maxit)



% Allocate initial point
x(:,1)=x0;

H=H0;

% Calculate function values of initial point
[~,fgrad] = feval(fhandle,x0);
df(:,1)=fgrad;

normdf(:,1)=norm(df(:,1));

i=1; 

while( i < maxit &&  norm(fgrad) > tol)

    % Compute the search direction
    p=-H*df(:,i); 
      
    % Compute step length
    alpha = armijo(fhandle,x(:,i),p,alpha0,c,beta,amax);
    
    % Update the point
    x(:,i+1)=x(:,i)+alpha*p;

    % Compute gradient
    [~,fgrad] = feval(fhandle,x(:,i+1));
    
    df(:,i+1)=fgrad; 
    
    % Compute norm of gradient
    normdf(:,i+1)=norm(df(:,i+1)); 
    
    s=x(:,i+1)-x(:,i);
    
    y=df(:,i+1)-df(:,i);
    
    %Inverse update formula
    H=H+(((s-H*y)*(s-H*y)')/((s-H*y)'*y)); 
    
    i=i+1;
end

X=x;           %Roots of f(x)
Df=df;         %Gradients of f(x) at roots 
NormDf=normdf; %Norm of Gradients of f(x)
it=i-1;        %Number of iterations


end