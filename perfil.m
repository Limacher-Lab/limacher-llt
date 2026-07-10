% This function is called by plot2_shape3D.m to plot blade shape.

function [xm,ym] = perfil(filename)
%-------------------------------
%XiYi = load('NACA0012.dat');
%XiYi = load('NACA644421.txt');
XiYi = readmatrix(filename);
%-------------------------------
    nx = length(XiYi(:,1));
    ny = length(XiYi(:,2));
%-------------------------------

x = XiYi(:,1);
y = XiYi(:,2);

for i = 1:nx-1
    Ax(i)  = (x(i+1)-x(i))*y(i);
    Ay(i)  = (y(i+1)-y(i))*x(i);
    xA(i) = (x(i+1)+x(i))/2*Ax(i);
    yA(i) = (y(i+1)+y(i))/2*Ay(i);
end

xm = sum(xA)/sum(Ax);
ym = sum(yA)/sum(Ay);

%plot(XiYi(:,1),XiYi(:,2),xm,ym,'r*');axis equal;grid on;