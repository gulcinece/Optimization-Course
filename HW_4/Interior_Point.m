function [criteria,iter,X,Y,S] = Interior_Point(A,b,C,X0,Y0,S0,Nmax,tol)

% function [criteria , X] = Interior_Point(A,b,C,X0,Y0,S0,Nmax)
%    This function implments the basic interior point algorithm.
%    Consider min  c'*X    Primal problem
%             s.t. Ax = b 
%            
%             max  b'*Y     Dual problem
%                  A'*Y <=c iff A'*Y + S  = c
%                                        S >= 0
%
% The input:
%   A: This is the Ax=b matrix. 
%   b: Vector. This is the right hand side of Ax=b.
%   c: Vector. This is from minimize  J(x) = c'x. 
%   X0: initial vector
%   S0: initial slack vector
%   Y0: initial vector
%
% Ouput:
%   criteria = c'*X
%   X:optimal value

% Diagonal Matrices
X = diag(X0) ; 
S = diag(S0) ; 

% size of matrices
n=size(X,1);
m=size(Y0,1);
e=ones(n,1);

% Zeros matrix in DF
H=zeros(n,m);
H1=zeros(n,n);
H2=zeros(m,m);

I=eye(n);
sigma=0.001;
k = 1;
err = 1;
while(norm(err)>tol && k<=Nmax) 
    
    mu=(X0'*S0/n);
    
    F1 = zeros(size(S0,1),1);
    F2 = zeros(size(b,1),1);
    F3 = X*S*e - sigma*mu.*e;
    
    % F(x,\lambda,s)
    F = [ F1;F2;F3];
    
    % Jacobian of F
    DF = [H1  A'   I ;
          A   H2   H';
          S   H    X  ];
    % Search Direction
    P = -DF\F;
    
    dX = P(1:n); % 
    dY = P(n+1:n+m);
    dS = P(n+m+1:end);
  
    alpha = step_size(X0,dX,S0,dS,tol);   
    
    % Update all variables
    X0 = X0 + alpha*dX ; 
    S0 = S0 + alpha*dS ; 
    Y0 = Y0 + alpha*dY ; 
    
    err = alpha*dX;
    
    X = diag(X0) ; 
    S = diag(S0) ; 
    k = k + 1;
end
% optimal point
iter = k;
X=diag(X);
Y=Y0;
S=diag(S);
% optimal value
criteria = C'*X ; 

end

% One choice of step-size
function [alpha] = step_size(X0,dX,S0,dS,tol)
% function  [alpha] = step_size(X0,dX,S0,dS,n)
% This function calculates the suitable step size alpha.

% Initial step length 
alpha1 = -1/min(min(dX./X0),-1);
alpha2 = -1/min(min(dS./S0),-1);
alpha_min = min(alpha1, alpha2);
p = size(X0,1)*1000;
fun = @(var)p*log((X0+var*dX)'*(S0+var*dS))-sum(log((X0+var*dX).*(S0+var*dS)));
alpha = fminbnd(fun,0,alpha_min);
% eta = 0.999;
% alphax = -1/min(min(dX./X0),-1);
% alphax = min(1, eta * alphax);
% alphas = -1/min(min(dS./S0),-1);
% alphas = min(1, eta * alphas);
% alpha = min(alphax, alphas);

% alpha = 5;

X = X0 + alpha * dX;
S = S0 + alpha * dS;

% temp = norm(A*X-b);
% temp = X.*S;
% temp2 = A*X - b;

% While condition fails to hold
% while(any(temp < 1e-6*(X0'*S0)/n))
% any(temp2 < 1e-6*(A*X0-b))
% any(X<=1e-6) || any(S<=1e-6)
while(any(X<=tol) || any(S<=tol))
    
    % Decrease the alpha 
    alpha = alpha * 0.5;
    
    % Calculate new X and S
    X = X0 + alpha * dX;
    S = S0 + alpha * dS;
    
%     temp = X.*S;
%     temp2 = A*X-b;
end

end

