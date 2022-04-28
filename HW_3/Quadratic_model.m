% Function of Quadratic Model
function [m] = Quadratic_model(fhandle,x,p)

[f,~,~]=feval(fhandle,x);
[~,gradf,~]=feval(fhandle,x);
[~,~,Hess]=feval(fhandle,x);

m = f + gradf'*p + 0.5*p'*Hess*p;

end