clear;
clc;
close all;

%% 1. Outdoor Transmission of COVID-19
% Part (a)
% Assume no change in z

D= 1; % m^2 s^-1
t_a=linspace(0,10800,10800);
M=98/3600*4;
x_a1=-8;
y_a1=8;
x_a3=0;
y_a3=8;
x_a7=-4;
y_a7=4;
x_a13=1;
y_a13=1;
h=1.8;
C1=function1(x_a1, y_a1, 0 ,M , D, t_a,h);
C3=function1(x_a3, y_a3, 0 ,M , D, t_a,h);
C7=function1(x_a7, y_a7, 0 ,M , D, t_a,h);
C13=function1(x_a13, y_a13, 0 ,M , D, t_a,h);

figure;

loglog(t_a, C1, t_a, C3, t_a, C7, t_a, C13)

legend('Concentration at 1','Concentration at 3','Concentration at 7','Concentration at 13');
grid on

xlabel('log (time[s])');
ylabel('log (concentration [g/s^3])');
title('Temporal evolution of virus concentration');

%% Part (b)
x_b = linspace (-9, 9);
y_b = linspace (-9, 9);
z_b=0;
t_b=10800;

fxy=zeros(length(x_b), length (y_b));
for j=1:length(x_b)
    for k=1:length(y_b)
        fxy(j,k)=function1(x_b(j),y_b(k),0, M, D, t_b,h);
    end
end



figure;

contourf(x_b, y_b, fxy)
xlabel('x[m]')
ylabel('y[m]')
zlabel('concentration')
title('Spatial distribution of virus concentration at head height (b)')

colorbar

U=3;
fxy2=zeros(length(x_b), length (y_b));
for j=1:length(x_b)
    for k=1:length(y_b)
        fxy2(j,k)=function2(x_b(j),y_b(k),0, M, D, t_b, U);
    end
end


figure;


contourf(x_b, y_b, fxy2)
xlabel('x[m]')
ylabel('y[m]')
zlabel('concentration')
title('Spatial distribution of virus concentration at head height (b)')


colorbar
%% Part (c) 
% 1 - No wind condition
x_c = linspace (-9, 9);
y_c = linspace (-9, 9);
z_c=0;
t_c=10800;

fxy_c=zeros(length(x_c), length (y_c));
changes_c=zeros(length(x_c), length (y_c));
for j=1:length(x_c)
    for k=1:length(y_c)
        fxy_c(j,k)=function3(x_b(j),y_b(k),z_c, M, D, t_b,  1.8);
        changes_c(j,k)=((fxy_c(j,k)-fxy(j,k))/fxy(j,k))*100;
    end
end

figure;
subplot(2,1,1)
surf(x_c, y_c, fxy_c)
xlabel('x[m]')
ylabel('y[m]')
zlabel('Concentration')
colorbar
title('Spatial distribution of virus concentration (with disinfectant) at head height, no wind condition')

subplot(2,1,2)
contourf(x_c, y_c, fxy_c)
xlabel('x[m]')
ylabel('y[m]')
zlabel('Concentration')
colorbar

figure;

contourf(x_c, y_c, changes_c)
xlabel('x[m]')
ylabel('y[m]')
zlabel('changes %')
title('Percentage changes in concentration with disinfectant applied to all ground surfaces')
colorbar

% Wind condition

fxy_c2=zeros(length(x_c), length (y_c));
changes_c2=zeros(length(x_c), length (y_c));
for j=1:length(x_c)
    for k=1:length(y_c)
        fxy_c2(j,k)=function4(x_b(j),y_b(k),z_c, M, D, t_b, U, 1.8);
        changes_c2(j,k)=((fxy_c2(j,k)-fxy2(j,k))/fxy2(j,k))*100;
    end
end

figure;
subplot(2,1,1)
surf(x_c, y_c, fxy_c2)
xlabel('x[m]')
ylabel('y[m]')
zlabel('Concentration')
colorbar
title('Spatial distribution of virus concentration (with disinfectant) at head height,  wind condition')

subplot(2,1,2)
contourf(x_c, y_c, fxy_c2)
xlabel('x[m]')
ylabel('y[m]')
zlabel('Concentration')
colorbar

figure;

contourf(x_c, y_c, changes_c2)
xlabel('x[m]')
ylabel('y[m]')
zlabel('changes %')
title('Percentage changes in concentration with disinfectant applied to all ground surfaces, Wind condition')
colorbar
%% Part (d)

%(i) No wind
z_d=linspace(-1.8,3);
x_d=0;
y_d=15;
t_d=10800;

for i=1:length(z_d)
    fxy3(i)=function1(x_d, y_d, z_d(i), M, D, t_d,h);
end

% (ii) Wind, 3m/s from South West
for i=1:length(z_d)
    fxy4(i)=function2(x_d, y_d, z_d(i), M, D, t_d, U);
end 

figure

subplot(1,2,1)
plot(fxy3, z_d);
xlabel('concentration');
ylabel('height');
grid on
title('No wind condition')

subplot(1,2,2)
plot(fxy4, z_d);
xlabel('concentration');
ylabel('height');
grid on
title('3m/s wind from Southwest condition')






%% Functions
% No wind
function [C] = function1(x, y, z, M, D, t, h)
r= sqrt(x.^2+y.^2+z.^2);
if r<1
    r=NaN;
end
r_2=sqrt(x.^2+y.^2+(z-2*h).^2);
C= (M./(4.*pi.*D.*r)).*erfc(r./sqrt(4*D*t))+(M./(4.*pi.*D.*r_2)).*erfc(r_2./sqrt(4*D*t));

end

% Wind, 3m/s from South West
function [C] = function2(x, y,z, M, D, t, U,h ) % Assume at t_b we have reached steady state
r= sqrt(x.^2+y.^2+z.^2);

if r<1
    r=NaN;
end
C= (M./(4.*pi.*D.*r)).*exp(-(U/(2*D))*(r-x));
end

% For (c), No wind condition
function [C] = function3(x, y,z, M, D, t, h)
r= sqrt(x.^2+y.^2+z.^2);
r_2=sqrt(x.^2+y.^2+(z-2*h).^2);
if r<1
    r=NaN;
end
if r_2<1
    r_2=NaN;
end
C= (M./(4.*pi.*D.*r)).*erfc(r./sqrt(4*D*t))-(M./(4.*pi.*D.*r_2)).*erfc(r_2./sqrt(4*D*t));
end

% For (c), with absorber, Wind condition
function [C] = function4(x, y,z, M, D, t, U, h) % Assume at t_b we have reached steady state
r= sqrt(x.^2+y.^2+z.^2);
r_2= sqrt(x.^2+y.^2+(z-2*h).^2);
if r<1
    r=NaN;
end
if r_2<1
    r_2=NaN;
end
C= (M./(4.*pi.*D.*r)).*exp(-(U/(2*D))*(r-x))-(M./(4.*pi.*D.*r_2)).*exp(-(U/(2*D))*(r_2-x));
end