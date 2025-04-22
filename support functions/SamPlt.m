function SamPlt(x)
% plot the samples
if ~iscell(x)
    figure
    theta = x;
    scatter(theta(1,:),theta(2,:),4,'filled',"MarkerFaceColor",[0.1,0.1,0.1],"LineWidth",2,"MarkerEdgeAlpha",0.5);
    xmin = min(theta(1,:),[],"all");
    xmax = max(theta(1,:),[],"all");
    ymin = min(theta(2,:),[],"all");
    ymax = max(theta(2,:),[],"all");
else
    Nite = length(x);
    figure
    xmin = +inf;
    xmax = -inf;
    ymin = +inf;
    ymax = -inf;
    for ite=1:Nite
        theta = x{ite};
        scatter(theta(1,:),theta(2,:),4,'filled',"LineWidth",2,"MarkerEdgeAlpha",0.5);
%         if ite<Nite
%             scatter(theta(1,:),theta(2,:),4,'filled',"LineWidth",2,"MarkerEdgeAlpha",0.2,"MarkerFaceColor",[0.8,0.8,0.8]);
%         else
%             scatter(theta(1,:),theta(2,:),4,'filled',"LineWidth",2,"MarkerEdgeAlpha",0.5);
%             Ns = size(theta,3);
%             color = linspace(0,1,Ns)';
%             for is = 1:Ns
%                 scatter(squeeze(theta(1,:,is)),squeeze(theta(2,:,is)),"Marker",".","MarkerEdgeColor",[color(is,:),0.1,0.1],"MarkerFaceColor",[color(is,:),0.1,0.1],"LineWidth",100,"MarkerFaceAlpha",1,"MarkerEdgeAlpha",1);
%                 hold on
%             end
%         end
        xmin = min(min(theta(1,:),[],"all"),xmin);
        xmax = max(max(theta(1,:),[],"all"),xmax);
        ymin = min(min(theta(2,:),[],"all"),ymin);
        ymax = max(max(theta(2,:),[],"all"),ymax);
        hold on
    end
end
ylim([ymin,ymax]);xlim([xmin,xmax]);
% ylim([80,200]);xlim([2,8]);
xlabel('$\theta_1$','Interpreter','latex');ylabel('$\theta_2$','Interpreter','latex');
set(gca,"Fontsize",20)
% legend('$i=0$','$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','$i=7$','$i=9$',"box","off",'interpreter','latex',"Location","northwest","FontSize",15)
% titlename = ['$','\rm{','CW-aMH','}\ \rm{','Example','\ ','d','}$'];
% titlename = ['$','\rm{','aCS','}\ \rm{','Example','\ ','c','}$'];
% titlename = ['$','\rm{','ES','}\ \rm{','Example','\ ','c','}$'];
% title(titlename,'interpreter','latex');
close all
end