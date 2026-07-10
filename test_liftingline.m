% Test lifting-line code

clear; clc; close all;

forcefile = "xf-n0012-il-100000_truncated.csv";
wingfilename = "blade_elliptic.txt";

 
alfas = [-22:22];
[CL,CD,y,Gamma,v] = liftingline(wingfilename,forcefile,alfas);

figure();
plot(alfas,CL,'b.-',alfas,CD,'r.-');
data = readmatrix(forcefile);
hold on;
plot(data(:,1),data(:,2),'c.-');
grid on;


