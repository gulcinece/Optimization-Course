
function [f,df] = f(x)
    f  = (x^2+1)*(x-1);
   df  = 3*x^2-2*x+1;
end

