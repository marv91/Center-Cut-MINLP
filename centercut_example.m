func    = [-1,-1];
A       = [-2.7 -1]; 
b       = [-27];
g1      = @(x) (x(1)-8)^2 + (x(2)-10)^2 - 20; 
intcon  = 2;
lb      = [0, 0];
ub      = [20, 20];
nonlcon = {g1};

sol     = centercut(func, intcon, A, b, lb, ub, nonlcon)

