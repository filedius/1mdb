clear;
clc;
close all;

%% Question 2 (a)
k=0.5;
w=6;
V=200;
G=1/3600;
r=0.0002;
T0=280;
T1=290;
g=9.81;
d=0.5;
r=0.0002;


N0=30;
d0=0.5;

b0=g*(T1-T0)/T0;
c0=0;
p0=0;

pmeters=[b0; c0; p0];
% [t, dxdt]=ode45(@f,t, INPUTS)
tspan=1:3600; %s

[t, x1]=ode45(@function1, tspan, pmeters,[], N0, V, k, w, d0, G,r);

Ti=((T0/g).*x1(:,1))+T0;


figure;
subplot(3,1,1)
plot(t,x1(:,2))
xlabel('t')
ylabel('c')

subplot(3,1,2)
plot(t,x1(:,3))
xlabel('t')
ylabel('p')

subplot(3,1,3)
plot(t,Ti)
xlabel('t')
ylabel('T')

d=0.1:0.01:1;

p=Inf;

for i =d
    N=30;
    n=1;
    p(end)=Inf;

    while p(end)>0.01
        [t, x]=ode45(@function1, tspan, pmeters,[], N, V, k, w, i, G,r);
        p=x(:,3);
        N=N-1;
    end
    Ti=((T0/g).*x(end,1))+T0;
   

    fprintf('d= %g, N= %g, T= %g \n', i, N+1, Ti);
    if Ti<290
        break
    end
end

fprintf('End of Question 2(a) \n');
        
%% 2 (b)
c0=0;
p0=0;
for i=1:4
    param=[b0;c0;p0];
    [t1, xb1]=ode45(@function1, 1:3600, param, [], N0, V, k, w, d0, G,r);
    [t2, xb2]=ode45(@function1, 1:1800, [xb1(end,1),xb1(end,2),xb1(end,3)],[], 0, V, k, w, d0, 0,r);
    
    param=[xb2(end,1),xb2(end,2),xb2(end,3)];
end

 
for i =d
    N=30;
    n=1;
    p(end)=Inf;

    while p(end)>0.01
        for j=1:4
            param=[b0;c0;p0];
            [t1, xb1]=ode45(@function1, 1:3600, param, [], N, V, k, w, i, G,r);
            [t2, xb2]=ode45(@function1, 1:1800, [xb1(end,1),xb1(end,2),xb1(end,3)],[], 0, V, k, w, i, 0,r);
            param=[xb2(end,1),xb2(end,2),xb2(end,3)];
        end


        p=xb2(:,3);
        N=N-1;
    end
    Ti=((T0/g).*xb2(end,1))+T0;
   

    fprintf('d= %g, N= %g, T= %g \n', i, N+1, Ti);
    if Ti<290
        break
    end
end
    
    
    
    

%% Functions
function x= function1(t, pmeters, N, V, k, w, d, G,r )
b = pmeters(1); c=pmeters(2); p=pmeters(3);

dbdt=((0.0028*N +0.028)/V)-(((k/3)*w*sqrt(b*d^3))/V)*b;
dcdt=(G/V)-(((k/3)*w*sqrt(b*d^3))/V)*c;
dpdt=r*(N-1)*(1-p)*c;

x=[dbdt;dcdt;dpdt];
end

