clc;
clear all;
close all;
theta=0:pi/45:pi;
theta_num=size(theta,2);
phi=0:pi/45:2*pi;
phi_num=size(phi,2);
[Theta,Phi]=meshgrid(theta,phi);
% P-wave radiation pattern
RP_P=zeros(phi_num,theta_num);
RP_Px=zeros(phi_num,theta_num);
RP_Py=zeros(phi_num,theta_num);
RP_Pz=zeros(phi_num,theta_num);
% S-wave radiation pattern
RP_S=zeros(phi_num,theta_num);
RP_Sx=zeros(phi_num,theta_num);
RP_Sy=zeros(phi_num,theta_num);
RP_Sz=zeros(phi_num,theta_num);
% Direction cosin vector
e_r=zeros(phi_num,theta_num,3);
e_theta=zeros(phi_num,theta_num,3);
e_phi=zeros(phi_num,theta_num,3);
% Different kinds of seismic moment tensor
M_ISO=1/3*[1 0 0;0 1 0;0 0 1];
M_LD=[1 0 0;0 0 0;0 0 0];
M_CLVD=sqrt(6)/6*[1 0 0;0 -2 0;0 0 1];
M_DC=sqrt(2)/2*[0 1 0;1 0 0;0 0 0];
% M_DC=sqrt(2)/2*[0 0 1;0 0 0;1 0 0];
M=zeros(3,3,4);
M(:,:,1)=M_ISO;
M(:,:,2)=M_LD;
M(:,:,3)=M_CLVD;
M(:,:,4)=M_DC;
delta=[1 0 0;0 1 0;0 0 1];
% Calculate the radiation pattern and plot it
for M_Num=1:4
    for j=1:theta_num
        for i=1:phi_num
            %         Y(i,j)=sin(Theta(i,j))*cos(Phi(i,j));
            %         X(i,j)=sin(Theta(i,j))*sin(Phi(i,j));
            %         Z(i,j)=cos(Theta(i,j));
            e_r(i,j,1)=sin(Theta(i,j))*cos(Phi(i,j));
            e_r(i,j,2)=sin(Theta(i,j))*sin(Phi(i,j));
            e_r(i,j,3)=cos(Theta(i,j));
            E_R=reshape(e_r(i,j,:),1,3);
            
            %             Calculate the P-wave radiation
            RP_P(i,j)=E_R*M(:,:,M_Num)*E_R';
            RP_Px(i,j)=e_r(i,j,1)*RP_P(i,j);
            RP_Py(i,j)=e_r(i,j,2)*RP_P(i,j);
            RP_Pz(i,j)=e_r(i,j,3)*RP_P(i,j);
            %             Calculate the S-wave radiation
            RP_Sx(i,j)=(delta(1,1)-e_r(i,j,1)*e_r(i,j,1))*e_r(i,j,1)*M(1,1,M_Num)+...
                (delta(1,1)-e_r(i,j,1)*e_r(i,j,1))*e_r(i,j,2)*M(1,2,M_Num)+...
                (delta(1,1)-e_r(i,j,1)*e_r(i,j,1))*e_r(i,j,3)*M(1,3,M_Num)+...
                (delta(1,2)-e_r(i,j,1)*e_r(i,j,2))*e_r(i,j,1)*M(2,1,M_Num)+...
                (delta(1,2)-e_r(i,j,1)*e_r(i,j,2))*e_r(i,j,2)*M(2,2,M_Num)+...
                (delta(1,2)-e_r(i,j,1)*e_r(i,j,2))*e_r(i,j,3)*M(2,3,M_Num)+...
                (delta(1,3)-e_r(i,j,1)*e_r(i,j,3))*e_r(i,j,1)*M(3,1,M_Num)+...
                (delta(1,3)-e_r(i,j,1)*e_r(i,j,3))*e_r(i,j,2)*M(3,2,M_Num)+...
                (delta(1,3)-e_r(i,j,1)*e_r(i,j,3))*e_r(i,j,3)*M(3,3,M_Num);
            
            RP_Sy(i,j)=(delta(2,1)-e_r(i,j,2)*e_r(i,j,1))*e_r(i,j,1)*M(1,1,M_Num)+...
                (delta(2,1)-e_r(i,j,2)*e_r(i,j,1))*e_r(i,j,2)*M(1,2,M_Num)+...
                (delta(2,1)-e_r(i,j,2)*e_r(i,j,1))*e_r(i,j,3)*M(1,3,M_Num)+...
                (delta(2,2)-e_r(i,j,2)*e_r(i,j,2))*e_r(i,j,1)*M(2,1,M_Num)+...
                (delta(2,2)-e_r(i,j,2)*e_r(i,j,2))*e_r(i,j,2)*M(2,2,M_Num)+...
                (delta(2,2)-e_r(i,j,2)*e_r(i,j,2))*e_r(i,j,3)*M(2,3,M_Num)+...
                (delta(2,3)-e_r(i,j,2)*e_r(i,j,3))*e_r(i,j,1)*M(3,1,M_Num)+...
                (delta(2,3)-e_r(i,j,2)*e_r(i,j,3))*e_r(i,j,2)*M(3,2,M_Num)+...
                (delta(2,3)-e_r(i,j,2)*e_r(i,j,3))*e_r(i,j,3)*M(3,3,M_Num);
            
            RP_Sz(i,j)=(delta(3,1)-e_r(i,j,3)*e_r(i,j,1))*e_r(i,j,1)*M(1,1,M_Num)+...
                (delta(3,1)-e_r(i,j,3)*e_r(i,j,1))*e_r(i,j,2)*M(1,2,M_Num)+...
                (delta(3,1)-e_r(i,j,3)*e_r(i,j,1))*e_r(i,j,3)*M(1,3,M_Num)+...
                (delta(3,2)-e_r(i,j,3)*e_r(i,j,2))*e_r(i,j,1)*M(2,1,M_Num)+...
                (delta(3,2)-e_r(i,j,3)*e_r(i,j,2))*e_r(i,j,2)*M(2,2,M_Num)+...
                (delta(3,2)-e_r(i,j,3)*e_r(i,j,2))*e_r(i,j,3)*M(2,3,M_Num)+...
                (delta(3,3)-e_r(i,j,3)*e_r(i,j,3))*e_r(i,j,1)*M(3,1,M_Num)+...
                (delta(3,3)-e_r(i,j,3)*e_r(i,j,3))*e_r(i,j,2)*M(3,2,M_Num)+...
                (delta(3,3)-e_r(i,j,3)*e_r(i,j,3))*e_r(i,j,3)*M(3,3,M_Num);
        end
    end
    Facecolor=[0.5 0.4 1];
    Edgecolor=[0 0 0];
    subplot(2,4,M_Num)
    view(3)
    hold on
    %     surf(RP_Px,RP_Py,RP_Pz,'edgecolor',Edgecolor);
    surf(RP_Px(:,1:23),RP_Py(:,1:23),RP_Pz(:,1:23),'edgecolor',Edgecolor);
    surf(RP_Px(:,23:46),RP_Py(:,23:46),RP_Pz(:,23:46),'edgecolor',Edgecolor);

    axis equal
    axis off
%     surf(RP_Px,RP_Py,RP_Pz,'facecolor',[0.5 0.4 1],'edgecolor',[0,0,0]);
    
    grid off
%     axis square
    
    subplot(2,4,4+M_Num)
    hold on
    %     surf(RP_Sx,RP_Sy,RP_Sz,'facecolor',Facecolor,'edgecolor',Edgecolor);
%     plot3(RP_Sx,RP_Sy,RP_Sz,'b.');
    surf(RP_Sx(:,1:23),RP_Sy(:,1:23),RP_Sz(:,1:23),'facecolor',Facecolor,'edgecolor',Edgecolor);
    surf(RP_Sx(:,23:46),RP_Sy(:,23:46),RP_Sz(:,23:46),'facecolor',Facecolor,'edgecolor',Edgecolor);
    axis([-1 1 -1 1 -1 1]);
    view(3)
    grid off
%     axis square
end




