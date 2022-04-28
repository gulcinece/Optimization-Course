

x=[-0.5,1];
tol=[1.0e-3,1.0e-6,1.0e-9];
maxit=10000;

[X1,Grad1,it1] = steepest_descent(@three_hump_camel,x,tol(1),maxit);
[X2,Grad2,it2] = steepest_descent(@three_hump_camel,x,tol(2),maxit);
[X3,Grad3,it3] = steepest_descent(@three_hump_camel,x,tol(3),maxit);

figure(1)
hold on
semilogy(0:it1-1,Grad1,'-o')
semilogy(0:it2-1,Grad2,'-s')
semilogy(0:it3-1,Grad3,'-d')
hold off
xlabel('iterations','fontsize',18)
ylabel('Norm of Gradient','fontsize',18)
title('Steepest Descent','fontsize',18)