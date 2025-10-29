function [f] = Rosenbrock(x, options)
    x1 = x(1); x2 = x(2);
    
    if options == 1

        f = 100*(x2 - x1^2)^2 + (1 - x1)^2;
        
    else
        f = [ -400*x1*(x2 - x1^2) - 2*(1 - x1);
               200*(x2 - x1^2) ];
    end
end
