function [m, b] = fo_taylor(func,x0)
m = Grad(func,x0);
b = -(func(x0) - dot(m,x0));
end

function g = Grad(fun, x0)
delta = x0 / 1000;             
for i = 1 : length(x0)
    if x0(i) == 0
        delta(i) = 1e-12;      
    end
    u = x0;                    
    u(i) = x0(i) + delta(i);
    f1 = feval(fun,u);     
    u(i) = x0(i) - delta(i);
    f2 = feval(fun,u);     
    g(i,1) = (f1 - f2) / (2 * delta(i)); 
end
end
