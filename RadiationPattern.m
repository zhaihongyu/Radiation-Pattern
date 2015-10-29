clc;
clear all;
close all;
theta=0:pi/20:pi;
theta_num=size(theta,2);
phi=0:pi/20:2*pi;
phi_num=size(phi,2);
[Theta,Phi]=meshgrid(theta,phi);
% P-wave radiation pattern
RadiationPattern_PX=zeros(theta_num,phi_num);
RadiationPattern_PY=zeros(theta_num,phi_num);
RadiationPattern_PZ=zeros(theta_num,phi_num);
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
% RadiationPattern_P=zeros(3,theta_num*phi_num);
% Calculate the P wave radiation pattern and plot it
figure
for M_Idx=1:4
    
    for j=1:theta_num
        for i=1:phi_num
            
            X(i,j)=sin(Theta(i,j))*cos(Phi(i,j));
            Y(i,j)=sin(Theta(i,j))*sin(Phi(i,j));
            Z(i,j)=cos(Theta(i,j));
            %% First method (Correct!!!)
            %{
            Direction_Cosine=[X(i,j),Y(i,j),Z(i,j)];
            RadiationPattern_P=Direction_Cosine'*Direction_Cosine*M(:,:,M_Idx)*Direction_Cosine';
            RadiationPattern_PX(j,i)=RadiationPattern_P(1);
            RadiationPattern_PY(j,i)=RadiationPattern_P(2);
            RadiationPattern_PZ(j,i)=RadiationPattern_P(3);
            %}
            %% Second emthod (Correct!!!)
            Ele_Vector_CC=[X(i,j),0,0;0,Y(i,j),0;0,0,Z(i,j)];
            CC_to_SC=...
                [sin(Theta(i,j))*cos(Phi(i,j)),sin(Theta(i,j))*sin(Phi(i,j)),cos(Theta(i,j));
                cos(Theta(i,j))*cos(Phi(i,j)),cos(Theta(i,j))*sin(Phi(i,j)),-sin(Theta(i,j));
                -sin(Phi(i,j)),cos(Phi(i,j)),0;];
            %             Ele_Vector_SC=CC_to_SC*Ele_Vector_CC;
            RadiationPattern_P=CC_to_SC(1,:)*M(:,:,M_Idx)*CC_to_SC(1,:)';
            RadiationPattern_PX(j,i)=RadiationPattern_P*sin(Theta(i,j))*cos(Phi(i,j));
            RadiationPattern_PY(j,i)=RadiationPattern_P*sin(Theta(i,j))*sin(Phi(i,j));
            RadiationPattern_PZ(j,i)=RadiationPattern_P*cos(Theta(i,j));
            %}
        end
    end
    Facecolor=[0.5 0.4 1];
    Edgecolor=[0 0 0];
%     figure
    subplot(4,1,M_Idx)
    view(3)
    hold on
        surf(RadiationPattern_PX,RadiationPattern_PY,RadiationPattern_PZ,'edgecolor',Edgecolor,'facecolor',Facecolor);
%     surf(RP_Px(:,1:23),RP_Py(:,1:23),RP_Pz(:,1:23),'edgecolor',Edgecolor);
%     surf(RP_Px(:,23:46),RP_Py(:,23:46),RP_Pz(:,23:46),'edgecolor',Edgecolor);
    axis([-1 1 -1 1 -1 1]);
%     surf(RP_Px,RP_Py,RP_Pz,'facecolor',[0.5 0.4 1],'edgecolor',[0,0,0]);
    axis equal
    axis off
    grid off
end

% Calculate the S wave radiation pattern and plot it
MT_Vector=zeros(6,1);
RadiationPattern_SX=zeros(theta_num,phi_num);
RadiationPattern_SY=zeros(theta_num,phi_num);
RadiationPattern_SZ=zeros(theta_num,phi_num);

RadiationPattern_SVX=zeros(theta_num,phi_num);
RadiationPattern_SVY=zeros(theta_num,phi_num);
RadiationPattern_SVZ=zeros(theta_num,phi_num);

RadiationPattern_SHX=zeros(theta_num,phi_num);
RadiationPattern_SHY=zeros(theta_num,phi_num);
RadiationPattern_SHZ=zeros(theta_num,phi_num);
figure
for M_Idx=1:4
    %Transform the matrix into vector
    MT_Vector(1:3)=M(1,:,M_Idx);
    MT_Vector(4:5)=M(2,2:3,M_Idx);
    MT_Vector(6)=M(3,3,M_Idx);
    for j=1:theta_num
        for i=1:phi_num
            X(i,j)=sin(Theta(i,j))*cos(Phi(i,j));
            Y(i,j)=sin(Theta(i,j))*sin(Phi(i,j));
            Z(i,j)=cos(Theta(i,j));
            %% First method (BUG!!!)
            %
            Direction_Cosine=[X(i,j),Y(i,j),Z(i,j)];
            %     Generate the first line of S-wave matrix A
            S_A11=Direction_Cosine(1 )^3-Direction_Cosine(1 );
            S_A12=2*Direction_Cosine(1 )^2*Direction_Cosine(2 )-Direction_Cosine(2 );
            S_A13=2*Direction_Cosine(1 )^2*Direction_Cosine(3 )-Direction_Cosine(3 );
            S_A14=Direction_Cosine(1 )*Direction_Cosine(2 )^2;
            S_A15=2*Direction_Cosine(1 )*Direction_Cosine(2 )*Direction_Cosine(3 );
            S_A16=Direction_Cosine(1 )*Direction_Cosine(3 )^2;
            S_A(1,:)=-[S_A11 S_A12 S_A13 S_A14 S_A15 S_A16];
            %     Generate the second line of S-wave matrix A
            S_A21=Direction_Cosine(1 )^2*Direction_Cosine(2 );
            S_A22=2*Direction_Cosine(1 )*Direction_Cosine(2 )^2-Direction_Cosine(1 );
            S_A23=2*Direction_Cosine(1 )*Direction_Cosine(2 )*Direction_Cosine(3 );
            S_A24=Direction_Cosine(2 )^3-Direction_Cosine(2 );
            S_A25=2*Direction_Cosine(2 )^2*Direction_Cosine(3 )-Direction_Cosine(3 );
            S_A26=Direction_Cosine(2 )*Direction_Cosine(3 )^2;
            S_A(2,:)=-[S_A21 S_A22 S_A23 S_A24 S_A25 S_A26];
            %     Generate the second line of S-wave matrix A
            S_A31=Direction_Cosine(1 )^2*Direction_Cosine(3 );
            S_A32=2*Direction_Cosine(1 )*Direction_Cosine(2 )*Direction_Cosine(3 );
            S_A33=2*Direction_Cosine(1 )*Direction_Cosine(3 )^2-Direction_Cosine(1 );
            S_A34=Direction_Cosine(2 )^2*Direction_Cosine(3 );
            S_A35=2*Direction_Cosine(2 )*Direction_Cosine(3 )^2-Direction_Cosine(2 );
            S_A36=Direction_Cosine(3 )^3-Direction_Cosine(3 );
            S_A(3,:)=-[S_A31 S_A32 S_A33 S_A34 S_A35 S_A36];
            
            RadiationPattern_S=[S_A(1,:)*MT_Vector,S_A(2,:)*MT_Vector,S_A(3,:)*MT_Vector];
            RadiationPattern_SX(j,i)=Direction_Cosine(1)*norm(RadiationPattern_S);
            RadiationPattern_SY(j,i)=Direction_Cosine(2)*norm(RadiationPattern_S);
            RadiationPattern_SZ(j,i)=Direction_Cosine(3)*norm(RadiationPattern_S);
            %}
            
            %% Second method (Correct!!!)
            %
            CC_to_SC=...
                [sin(Theta(i,j))*cos(Phi(i,j)),sin(Theta(i,j))*sin(Phi(i,j)),cos(Theta(i,j));
                cos(Theta(i,j))*cos(Phi(i,j)),cos(Theta(i,j))*sin(Phi(i,j)),-sin(Theta(i,j));
                -sin(Phi(i,j)),cos(Phi(i,j)),0;];
            %SV Radiation Pattern
            RadiationPattern_SV=CC_to_SC(2,:)*M(:,:,M_Idx)*CC_to_SC(1,:)';
            RadiationPattern_SVX(j,i)=RadiationPattern_SV*sin(Theta(i,j))*cos(Phi(i,j));
            RadiationPattern_SVY(j,i)=RadiationPattern_SV*sin(Theta(i,j))*sin(Phi(i,j));
            RadiationPattern_SVZ(j,i)=RadiationPattern_SV*cos(Theta(i,j));
            %SH Radiation Pattern
            RadiationPattern_SH=CC_to_SC(3,:)*M(:,:,M_Idx)*CC_to_SC(1,:)';
            RadiationPattern_SHX(j,i)=RadiationPattern_SH*sin(Theta(i,j))*cos(Phi(i,j));
            RadiationPattern_SHY(j,i)=RadiationPattern_SH*sin(Theta(i,j))*sin(Phi(i,j));
            RadiationPattern_SHZ(j,i)=RadiationPattern_SH*cos(Theta(i,j));
            %S Radiation Pattern
            RadiationPattern_S=norm([RadiationPattern_SV,RadiationPattern_SH]);
            RadiationPattern_SX(j,i)=RadiationPattern_S*sin(Theta(i,j))*cos(Phi(i,j));
            RadiationPattern_SY(j,i)=RadiationPattern_S*sin(Theta(i,j))*sin(Phi(i,j));
            RadiationPattern_SZ(j,i)=RadiationPattern_S*cos(Theta(i,j));
            %}
        end
    end
    Facecolor=[0.5 0.4 1];
    Edgecolor=[0 0 0];
    %{
    figure
    subplot(1,3,1)
    view(3)
    hold on
    surf(RadiationPattern_SVX,RadiationPattern_SVY,RadiationPattern_SVZ,'edgecolor',Edgecolor,'facecolor',Facecolor);
%     surf(RadiationPattern_SVX(1:11,:),RadiationPattern_SVX(1:11,:),abs(RadiationPattern_SVX(1:11,:)),'edgecolor',Edgecolor);
    axis([-1 1 -1 1 -1 1]);
    axis equal
    subplot(1,3,2)
    view(3)
    hold on
    surf(RadiationPattern_SHX,RadiationPattern_SHY,RadiationPattern_SHZ,'edgecolor',Edgecolor,'facecolor',Facecolor);
%     surf(RadiationPattern_SVX(1:11,:),RadiationPattern_SVX(1:11,:),abs(RadiationPattern_SVX(1:11,:)),'edgecolor',Edgecolor);
    axis([-1 1 -1 1 -1 1]);
    axis equal
    %}
    subplot(4,1,M_Idx)
%     subplot(1,3,3)
    view(3)
    hold on
    surf(RadiationPattern_SX,RadiationPattern_SY,RadiationPattern_SZ,'edgecolor',Edgecolor,'facecolor',Facecolor);
%     surf(RadiationPattern_SVX(1:11,:),RadiationPattern_SVX(1:11,:),abs(RadiationPattern_SVX(1:11,:)),'edgecolor',Edgecolor);
    axis([-1 1 -1 1 -1 1]);
    axis equal
    axis off
    grid off
end


