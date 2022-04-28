

x0     =[0,-1];
tau1   = 1/4;
tau2   = 3/4;
tol    = 1.0e-2;
maxit  = 10000;
delta0 =1;
 
[XX,Grad,it,p] = Trust_Region(@rosenbrock_hes,x0,delta0,tau1,tau2,maxit,tol);
figure(1)
plot(0:it,Grad,'--*')
xlabel('Number of iterations','fontsize',18)
ylabel('Norm of Gradient','fontsize',18)
title('Trust Region','fontsize',18)

figure(2)
p_x=linspace(-1.5,1.5,it-1);
p_x=p_x';
p_y=linspace(-1.5,1.5,it-1);
p_y=p_y';
[X,Y]=meshgrid(p_x,p_y);
m=zeros(it-1);
for i=1:400
 m(i) = Quadratic_model(@rosenbrock_hes,[0,-1],[X(i);Y(i)]);
end
contour(X,Y,m,'b','ShowText','on');
hold on
plot([0 0], [p(1,1), p(2,1)])
hold off
xlabel('p_1','fontsize',18)
ylabel('p_2','fontsize',18)
title('Contour Lines of Quadratic Model','fontsize',18)


