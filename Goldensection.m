function [alpha, fmin] = Goldensection(a,b,fun,g,sdir,epsilon,itmax,c1,c2)
% Golden Section Search on [a,b] for scalar fun(alpha)
    %check at lamda = 1 by wolf and Amijo's condition
    if fun(1) <= fun(0)+c1*sdir'*g(0) && ... 
        abs(sdir'*g(1)) <= -c2*sdir'*g(0)

        alpha = 1;
        fmin = fun(1) ;
        
    %find exact path : Lamda
    else
    phi = (1 + sqrt(5)) / 2;
    rho = phi - 1;
    x1 = b - rho*(b - a);
    x2 = a + rho*(b - a);
    f1 = fun(x1);
    f2 = fun(x2);
    for k = 1:itmax
        if f1 > f2
            a = x1;
            x1 = x2;  f1 = f2;
            x2 = a + rho*(b - a);
            f2 = fun(x2);
        else
            b = x2;
            x2 = x1;  f2 = f1;
            x1 = b - rho*(b - a);
            f1 = fun(x1);
        end
        if abs(b - a) < epsilon, break; end
    end

    alpha = (a + b) / 2;
    fmin  = fun(alpha);

    end

end
