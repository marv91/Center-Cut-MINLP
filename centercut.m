function [ret, time] = centercut(func, intcon, A_l, b_l, lb, ub, nonlcon, details)

time = 0;
 
if ~exist('details','var')
     details = 1;
end

size_A_l = size(A_l);
size_nonlcon = length(nonlcon);
dimension = length(ub);

A_l_2 = horzcat(A_l,zeros(size_A_l(1),1)); 
r_func = [zeros(length(func),1)
           -1];

c = @(x_) nonlinfun(nonlcon,x_);
ceq = [];
nonlinfcn = @(x_)deal(c(x_),ceq);       
       
x = [];
r = inf;
B_A = [];
B_b = [];
counter = 1;

while(r > 0.01)

    K = zeros(1, size_nonlcon);
    violated = 1;
    
if counter > 1
    
    violated = 0;
    for i = 1:size_nonlcon
        if(nonlcon{i}(x)>=0) 
            K(i) = 1; 
        end
        if(nonlcon{i}(x)>0) 
            violated = 1; 
        end
    end

end

if violated == 0

Aeq = zeros(dimension,dimension);
Aeq(intcon,intcon) = 1;
beq = zeros(dimension,1);
beq(intcon) = x(intcon);

obj_func = @(x) dot(x,func);

tic;
x_new = fmincon(obj_func,x,A_l,b_l,Aeq,beq,lb,ub,nonlinfcn,optimoptions('fmincon','Display', 'off'));
time = time + toc;

x = x_new;
val = dot(x,func);

B_A = vertcat(B_A,[func, norm(func,2)]);
B_b = vertcat(B_b,val);

for i = 1:size_nonlcon
        if(abs(nonlcon{i}(x))<1E-7)
            [m_i, b_i] = fo_taylor(nonlcon{i}, x);
            B_A = vertcat(B_A,[transpose(m_i), norm(m_i)]);
            B_b = vertcat(B_b,b_i);
        end
end

else    

for i = 1:size_nonlcon
 
        if(K(i)==1)
            [m_i, b_i] = fo_taylor(nonlcon{i}, x);
            B_A = vertcat(B_A,[transpose(m_i), norm(m_i)]);
            B_b = vertcat(B_b,b_i);
        end
    
end

end

tic;
sol = intlinprog(r_func,intcon,vertcat(B_A,A_l_2),vertcat(B_b,b_l),[],[],[lb 0],[ub max(ub)],[],optimoptions('intlinprog','Display', 'off'));
time = time + toc;

x = sol(1:dimension);
r = sol(dimension+1);

counter = counter + 1;
    
end

ret = x;


if details == 1
disp('Center-Cut Details');
disp(['----------------------------------------------------']);
disp(['number of iterations: ',num2str(counter),' iterations']);
disp(['time elapsed        : ',num2str(time),' seconds']);
disp(['----------------------------------------------------']);
end

end

function c = nonlinfun(nl,x)

NLEQ_n = length(nl);

for i = 1:NLEQ_n
    c(i) = nl{i}(x); 
end

end









