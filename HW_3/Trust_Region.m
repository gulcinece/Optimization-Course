
function [X,Grad,it,p]= Trust_Region(fhandle,x0,delta0,tau1,tau2,maxit,tol)

% Description: The function implements trust region algorithm with given
%   inputs and returns iteration points,norms of the gradients at these
%   points and iteration number.

% Inputs: 
%   fhandle: function handle
%   x0: initial guess
%   delta0: initial trust region radius
%   tau1,tau2: parameters
%   maxit: maximum iteration number
%   tol: tolerance for stopping criteria

% Outputs:
%   X: iteration points from initial value to convergence value
%   Grad: Norms of the gradient vector at iteration points
%   it: iteration number

% Usage: 
%   Trust_Region(fhandle,x0,delta0,tau1,tau2,maxit,tol)

% Allocate initial point
x(:,1)=x0;
% Calculate function values of initial point
[~,gradf(:,1)]=feval(fhandle,x(:,1));
Grad(:,1)=norm(gradf(:,1));

i=1;
delta=delta0;

while(Grad(:,i)>tol && i<maxit)
    
    [~,~,Hess]=feval(fhandle,x(:,i));
    
    % Compute the search direction minimizing the model function m(p)
    p(:,i)=-Hess\gradf(:,i); 
   
    if norm(p(:,i))<=delta  % if true then accept directon
        
        % Function value at current point
        [f,~,~]=feval(fhandle,x(:,i)); 
        % Function value at next iteration point
        [fnext,~,~]=feval(fhandle,(x(:,i)+p(:,i))); 
        % Model function value at next iteration piont
        m=Quadratic_model(fhandle,x(:,i),p(:,i)); 
        % Actual reduction over predicted reduction 
        rho=(f-fnext)/(f-m); 
        
        if rho<=tau1          % If true then iteration is unsuccessful
            x(:,i+1)=x(:,i);
            gradf(:,i+1)=gradf(:,i);
            Grad(:,i+1)=norm(gradf(:,i));
            delta=delta/2; 
            
        else                    % successful step
            x(:,i+1)=x(:,i)+p(:,i);  % update the point
            [~,gradf(:,i+1)]=feval(fhandle,x(:,i+1));
            Grad(:,i+1)= norm(gradf(:,i+1));   
            
            if (tau1<rho) && (rho<tau2)     % If true then do not change the radius
                delta=delta;
            end
            
            if rho>=tau2                 % If true then increase the radius
                delta=delta*2;
            end
            
        end
        
        
    else       % If norm(p)>delta we can't use newton direction
       
        % Find lamda satisfying norm(inverse(H+lamda*I)*g)=delta
        fun=@(lamda) norm((Hess+lamda*eye(size(Hess)))\gradf(:,i))-delta; 
        lamda=fzero(fun,[0 100]);
        
        % Find search direction using the lamda
        p(:,i)=-inv(Hess+lamda*eye(size(Hess)))*gradf(:,i);
        % Compute function value at current point
        [f,~,~]=feval(fhandle,x(:,i)); 
        
        % Compute function value at next iteration point
        [fnext,~,~]=feval(fhandle,(x(:,i)+p(:,i))); 
        % Compute model function m(p)
        m = Quadratic_model(fhandle,x(:,i),p(:,i));  
        % Find actual reduction over predicted reduction
        rho=(f-fnext)/(f-m); 
        
        if rho<=tau1
            x(:,i+1)=x(:,i);
            gradf(:,i+1)=gradf(:,i);
            Grad(:,i+1)=norm(gradf(:,i));
            delta=delta/2;
            
        else
            x(:,i+1)=x(:,i)+p(:,i);
            [~,gradf(:,i+1)]=feval(fhandle,x(:,i+1));
            Grad(:,i+1)= norm(gradf(:,i+1));
            if (tau1<rho) && (rho<tau2)
                delta=delta;
            end
            if rho>=tau2
                delta=delta*2;
            end
        end
    
    end
     i=i+1;
end

X=x;
it=i-1;
end



    