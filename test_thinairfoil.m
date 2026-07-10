% Test thinairfoil function


clear; clc; close all;

% Inputs
airfoilname = 'NACA0012_profile.txt';
forcefile = 'xf-n0012-il-100000.csv';

figure();
h = plotcamberline(airfoilname);

%%%
alfa = -10:2:16;

[cl,cm] = thinairfoil(airfoilname,alfa);

forcedata = readmatrix(forcefile);

figure();
plot(alfa,cl,'b.-');
hold on;
plot(forcedata(:,1), forcedata(:,2),'k');
legend('thin-airfoil theory','Xfoil','location','northoutside');
xlabel('$\alpha$','interpreter','latex');
ylabel('$c_l$', 'interpreter','latex');
set(gcf,'color','w');
grid on;

figure();
plot(alfa,cm,'b.-');
hold on;
plot(forcedata(:,1), forcedata(:,5),'k');
legend('thin-airfoil theory','Xfoil','location','northoutside');
xlabel('$\alpha$','interpreter','latex');
ylabel('$c_m$', 'interpreter','latex');
set(gcf,'color','w');
grid on;