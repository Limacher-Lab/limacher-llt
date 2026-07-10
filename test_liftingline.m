% Test lifting-line code

clear; clc; close all;

forcefile = "xf-n0012-il-100000_truncated.csv";
wingfilename = "blade_elliptic.txt";

% [~,~,y,Gamma,v] = liftingline(wingfilename,forcefile,2);
% 
% figure();
% yyaxis left;
% plot(y,Gamma,'b.-');
% ylim([0,1.1*max(Gamma)]);
% yyaxis right; 
% plot(y,v,'r.-');
% ylim([-1.1*max(abs(v)), 1.1*max(abs(v))])
% leg = legend('$\Gamma$','$v$');
% xlabel('$y^*$','interpreter','latex');
% set(leg,'interpreter','latex','location','northoutside');
% grid on;
% set(gcf,'color','w');


%%

 
alfas = [-20:20];
[CL,CD,y,Gamma,v] = liftingline(wingfilename,forcefile,alfas);

figure();
plot(alfas,CL,'b.-',alfas,CD,'r.-');
data = readmatrix(forcefile);
hold on;
plot(data(:,1),data(:,2),'c.-');
grid on;



% figure();
% [~,b] = size(Gamma);
% for ii = 1:b
%     plot(y,Gamma,'k-.');
%     drawnow;
%     pause(1);
% end
% 
% 
% figure();
% [~,b] = size(v);
% for ii = 1:b
%     plot(y,v,'k-.');
%     drawnow;
%     pause(1);
% end



