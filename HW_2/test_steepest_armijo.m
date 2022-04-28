
x=[-0.5,1];
tol=[1.0e-3,1.0e-6,1.0e-9];
maxit=10000;
alpha0 = 1; 
c = 1.0e-4; 
beta = 0.5; 
amax = 100; 

[X1,Grad1,it1] = steepest_descent_armijo(@three_hump_camel,x,tol(1),maxit, alpha0, c, beta, amax);
[X2,Grad2,it2] = steepest_descent_armijo(@three_hump_camel,x,tol(2),maxit, alpha0, c, beta, amax);
[X3,Grad3,it3] = steepest_descent_armijo(@three_hump_camel,x,tol(3),maxit, alpha0, c, beta, amax);

figure(1)
hold on
semilogy(0:it1-1,Grad1,'-o')
semilogy(0:it2-1,Grad2,'-s')
semilogy(0:it3-1,Grad3,'-d')
hold off
xlabel('iterations','fontsize',18)
ylabel('Norm of Gradient','fontsize',18)
title('Steepest Descent','fontsize',18)