# Center-Cut-MINLP
This Matlab code implements the Center-Cut Algorithm presented by Jan Kronqvist, David Bernal, Andreas Lundell and Tapio Westerlund. The algorithm solves convex MINLP Problems. The original paper can be found here: http://www.optimization-online.org/DB_HTML/2018/02/6456.html


Usage:

```matlab
func    = [-1,-1];
A       = [-2.7 -1]; 
b       = [-27];
g1      = @(x) (x(1)-8)^2 + (x(2)-10)^2 - 20; 
intcon  = 2;
lb      = [0, 0];
ub      = [20, 20];
nonlcon = {g1};

sol     = centercut(func, intcon, A, b, lb, ub, nonlcon)
```
