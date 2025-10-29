clc;clear;close all ;
epsilon = 1e-4; itmax =100;mu=1e-4;eta=0.1;
x0 = [10;10];   % case 1
[xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,nH,CHN,IFLAG,k] = BFGS(@Rosenbrock,x0,epsilon,mu,eta,itmax);



disp("xmin convert to")
disp(xmin)
% Display the minimum function value found
disp("fmin value is")
disp(fmin)
% Check if the optimization was successful
if IFLAG == 0
    disp('Optimization converged successfully.');
else
    disp('Optimization did not converge.');
end
figure;
% Plot the convergence of the function values
semilogy(1:k, Fk(1:k), 'o-');
xlabel('Iteration');
ylabel('Function Value');
title('Convergence of Function Values');
grid on;