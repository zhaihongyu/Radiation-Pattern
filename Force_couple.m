% Force couple representation
clc;
clear all;
close all;

h=figure;
set(h,'Position',[50 50 1000 1000])
% set(h,'Position',get(0,'Screensize'))
for i=1:3;
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        hold on;
        axis off
%         grid on;
% Plot 3-D coordinate system
        X=zeros(1,3);
        Y=zeros(1,3);
        Z=zeros(1,3);
        U=[1 0 0];
        V=[0 1 0];
        W=[0 0 1];
        quiver3(X,Y,Z,U,V,W,1.5);
        set(gca,'XDir','reverse','YDir','reverse','xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1])
%     Plot positive Force couple    
        FC_SPoint=[0 0 0];
        FC_SPoint(j)=FC_SPoint(j)+0.5;
        FC_SPoint_Rev=-FC_SPoint;
        FC_X=[FC_SPoint(1) FC_SPoint_Rev(1)];
        FC_Y=[FC_SPoint(2) FC_SPoint_Rev(2)];
        FC_Z=[FC_SPoint(3) FC_SPoint_Rev(3)];
        
        FC_UVW=zeros(3,2);
        FC_UVW(i,1)=1;
        FC_UVW(i,2)=-1;
        h1=quiver3(FC_X,FC_Y,FC_Z,FC_UVW(1,:),FC_UVW(2,:),FC_UVW(3,:),1);
        set(h1,'LineWidth',2.0,'color',[0 0 0]);
%         Plot dot line
        h2=plot3(FC_X,FC_Y,FC_Z,':k');
        set(h2,'LineWidth',2.0)
        
        xlabel('x');
        ylabel('y');
        zlabel('z');
        view(3);
        axis square;
    end
end