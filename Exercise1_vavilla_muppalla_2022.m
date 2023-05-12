%%
                          % Excercise 1 %
                   % Precise Point Positioning %
                   % Advanced concepts of Positioning and Navigation %
                   % Vavilla Sai Viswanath %
                   % Kanthi Lochan Muppalla % 
%%
clc;
clear;
% Loading the data
%load 'Exercise11.mat'
load 'ACPN_PPP_data.mat'
%%
rec = [3800689.4257 882077.6724 5028791.5720];
XYZ = zeros(480,3);
for i = 1:480
    XYZ(i,:) = rec;
end
% Epochs = (249:249+4*60*2-1);
Epochs = (1:480);

% Code obs
IFCOMC = OMC.CODE;
%IFCOMC(isnan(IFCOMC)) = 0;
C_IFCOMC=IFCOMC(Epochs,:); %code obs for given epochs
% Phase obs
IFPOMC = OMC.PHAS;
%IFPOMC(isnan(IFPOMC)) = 0;
P_IFPOMC=IFPOMC(Epochs,:); % phase obs for given epochs

[r,c]=find(~isnan(C_IFCOMC));
All_sat = unique(c)';
sat = length(All_sat);


[cix,ciy] = find(isnan(C_IFCOMC)); 

Az = SATPOS.AZI(Epochs,:);
El = SATPOS.ELEV(Epochs,:);
% Azimuth
%Az(isnan(Az)) = 0;
% Elevation
%El(isnan(El)) = 0;
Az_ind = sub2ind(size(Az),cix,ciy);
Az(Az_ind) = NaN;
El_ind = sub2ind(size(El),cix,ciy);
El(El_ind) = NaN;

% Trop mapping func
mfw = MFW;
%mfw(isnan(mfw)) = 0;
M_fw = mfw(Epochs,:); % For given epochs
M_fw_ind = sub2ind(size(M_fw),cix,ciy);
M_fw(M_fw_ind) = NaN;
% Satellite position 
X = SATPOS.X(Epochs,:);
Y = SATPOS.Y(Epochs,:);
Z = SATPOS.Z(Epochs,:);
X_ind = sub2ind(size(X),cix,ciy);
X(X_ind) = NaN;
Y_ind = sub2ind(size(Y),cix,ciy);
Y(Y_ind) = NaN;
Z_ind = sub2ind(size(Z),cix,ciy);
Z(Z_ind) = NaN;

%removing NaN columns
C_IFCOMC = C_IFCOMC(:,~all(isnan(C_IFCOMC)));
P_IFPOMC = P_IFPOMC(:,~all(isnan(P_IFPOMC)));
Az = Az(:,~all(isnan(Az)));
El = El(:,~all(isnan(El)));
M_fw = M_fw(:,~all(isnan(M_fw)));
X = X(:,~all(isnan(X)));
Y = Y(:,~all(isnan(Y)));
Z = Z(:,~all(isnan(Z)));

%% Task 1.a
% Visible satellites
figure(1)
plot(r,c,'k*','Markersize',1);
title('Visible satellites for given epochs')
xlabel('Epoch');
ylabel('PRN');
% Calculation of DOP values 
for i=1:480
    A{i}=[-cos(Az(i,:)).*cos(El(i,:)); -sin(Az(i,:)).*cos(El(i,:)); -sin(El(i,:))]';
    A{i} = [ A{i} ones(sat,1)];
%     A{i} = [ A{i} ones(32,1)];
    A{i}(any(isnan(A{i})'),:)=[];
    Q{i}=inv(A{i}'*A{i});
    Q_NN(i) = Q{i}(1,1);
    Q_EE(i) = Q{i}(2,2);
    Q_UU(i) = Q{i}(3,3);
    Q_TT(i) = Q{i}(4,4);
    GDOP(i)=sqrt(Q_NN(i) + Q_EE(i) + Q_UU(i) + Q_TT(i));              % Geomteric dilution of precision
    PDOP(i)=sqrt(Q_NN(i) + Q_EE(i) + Q_UU(i));                     % Position dilution of precision
    HDOP(i)=sqrt(Q_NN(i) + Q_EE(i));                            % Horizantal dilution of precision
    VDOP(i)=sqrt(Q_UU(i));                                   % Vertical dilution of precision
end
%GDOP=[GDOP{:}];
%PDOP=[PDOP{:}];
%HDOP=[HDOP{:}];
%VDOP=[VDOP{:}];
% Plotting DOP values
figure(2)
plot(GDOP,'r')
hold on
plot(PDOP,'g')
hold on 
plot(HDOP,'b')
hold on
plot(VDOP,'k')
title('DOP Values'),
xlabel('Epochs')
ylabel('DOP')
legend('GDOP','PODP','HDOP','VDOP');
%%
% Geometric distance
d = sqrt((X-XYZ(1,1)).^2+(Y-XYZ(1,2)).^2+(Z-XYZ(1,3)).^2);
% Partial derivatives of position with respect to geometric distance
A_pos_x = (XYZ(1,1)-X)./d;
A_pos_y = (XYZ(1,2)-Y)./d;
A_pos_z = (XYZ(1,3)-Z)./d;
%% Task 1.b
%% Design matrix setup 
% initialization
A_static = [];                                            % Design matix for static mode
A_kinematic = [];                                         % Design matrix for kinematic mode
cell_sta = cell(1,480);
cell_kin = cell(1,480);
% for i = 1:480
%     as = cat(1, A_pos_x(i,:),A_pos_y(i,:),A_pos_z(i,:))'; % Concatenating partial derivatives
%     as(find(isnan(as(:,1))),:)=[];                        % Finding non NaN values
%     Clk_sta = zeros(size(as,1),480);
%     Clk_sta(:,i)= 1;                                      % Receiver clock error in static mode
%     Clk_kin = zeros(size(as,1),1920);                     % Reciver Clock error in kinematic mode will also include coordinates 
%     P_k=[as ones(size(as,1),1)]; 
%     Clk_kin(:,4*i-3:4*i)=P_k;                             % Receiver clock error in kinematic mode
%     t = M_fw(i,:)'; 
%     t(t == 0) =[]; 
%     Tcomp=[zeros(size(t,1),4)];                           % Tropospheric component
%     Tcomp(:,ceil(i/120))=t;
%     C_amb = zeros(size(as,1),sat);                        % Ambiguity for code observations
%     P_amb = zeros(size(C_amb));                           % Ambiguity for phase observations
%     vs = find((Az(i,:)~= 0));                             % Visible satellites
%     for j = 1:length(vs)
%         vs_amb = find(All_sat == vs(j));
%         P_amb(j,vs_amb )=1;                               % Phase ambiguity
%     end
%     cell_sta{i} = [C_amb Tcomp as Clk_sta ; 
%                    P_amb Tcomp as Clk_sta];
%     A_static = [A_static; cell_sta{i}];                   % Design matrix for static mode
%     cell_kin{i} = [C_amb Tcomp Clk_kin ; 
%                    P_amb Tcomp Clk_kin];
%     A_kinematic =[A_kinematic ; cell_kin{i}];             % Design matrix for kinematic mode
% end

for i = 1:480
    as = cat(1, A_pos_x(i,:),A_pos_y(i,:),A_pos_z(i,:))'; % Concatenating partial derivatives
    as(find(isnan(as(:,1))),:)=[];                        % Finding non NaN values
    clk = zeros(size(as,1),480);
    clk(:,i)= 1;                                      % Receiver clock error in static mode 
    as_kin = zeros(size(as,1),1440);
    as_kin(:,3*i-2:3*i)=as;                             % Receiver clock error in kinematic mode
    t = M_fw(i,:)'; 
    t(isnan(t)) =[]; 
    Tcomp=[zeros(size(t,1),4)];                           % Tropospheric component
    Tcomp(:,ceil(i/120))=t;
    C_amb = zeros(size(as,1),sat);                        % Ambiguity for code observations
    P_amb = zeros(size(C_amb));                           % Ambiguity for phase observations
    vs = find(~isnan(Az(i,:)));                             % Visible satellites
    for j = 1:length(vs)
        vs_amb = find(All_sat == vs(j));
        P_amb(j,vs_amb )=1;                               % Phase ambiguity
    end
    cell_sta{i} = [as clk Tcomp C_amb; 
                   as clk Tcomp P_amb];
    A_static = [A_static; cell_sta{i}];                   % Design matrix for static mode
    cell_kin{i} = [as_kin clk Tcomp C_amb; 
                   as_kin clk Tcomp P_amb];
    A_kinematic =[A_kinematic ; cell_kin{i}];             % Design matrix for kinematic mode
end
%% Task 1.c
%% Setting up stachastic model with elevation dependent weighting 
Sc=3; 
Sp= 0.003;
% p=[];                                                 
% for j=1:480
%     i = find(~isnan(X(j,:)));
%     p=[p (sin(El(j,i))).^2/Sc^2 (sin(El(j,i))).^2/Sp^2];   
% end
% p=diag(p);
% % stochastic model
% Sm =inv(p);
P_cell = cell(480,1);
for i=1:480
vis_sat = sum(~isnan(C_IFCOMC(i,:)));
vis_prn = (find(~isnan(C_IFCOMC(i,:))))';
P_temp = zeros(vis_sat*2,1);
     for m = 1:vis_sat
         P_temp(m) = sin(El(i,vis_prn(m))^2)./Sc^2;
         P_temp(m+vis_sat) = sin(El(i,vis_prn(m))^2)./Sp^2;
     end
     P_cell{i} = P_temp;
end
P = cell2mat(P_cell);
% Weight matrix
p = diag(P);

%% Task 1.d
%% computing coordinate solution with One step solution for static and kinematic mode

r_pos=XYZ(1:480,:); % receiver position
numberOfZeros1 = sum(C_IFCOMC(:)~=0);
numberOfZeros2 = sum(P_IFPOMC(:)~=0);
numberOfZeros3 = sum(~isnan(d(:)));
numberOfZeros4 = sum(~isnan(A_pos_x(:)));
numberOfZeros5 = sum(~isnan(A_pos_y(:)));
numberOfZeros6 = sum(~isnan(A_pos_z(:)));
numberOfZeros7 = sum(M_fw(:)~=0);


% L= [];                      % observation vector
% a_0= [];                    % geometric distance
% for i=1:1:480
%     Vs = find((C_IFCOMC(i,:) ~= 0));   % visible satellites
%     L= [ L C_IFCOMC(i,Vs) P_IFPOMC(i,Vs)];  % observation vector (code and phase combined)
%     a_0=[ a_0 d(i,Vs) d(i,Vs)];  % computing correction
% end
% delta_l= L-a_0;      % Observed minus computed vector(OMC)
% delta_l = delta_l'; 

delta_l= [];
for i = 1:480
    cd = (C_IFCOMC(~isnan(C_IFCOMC(i,:))));
    %ph = nonzeros(P_IFPOMC(i,:));
    ph = (P_IFPOMC(~isnan(P_IFPOMC(i,:))));
    delta_l = [delta_l; cd'; ph'];
end
%  One step LSA for static and kinematic modes

% cofactor matrix of the estimated paramters
Qxx_static = inv(A_static'*p*A_static);                % for static mode 
Qxx_kinematic = inv(A_kinematic'*p*A_kinematic);       % for kinematic mode
% deviations wrt apriori values
dx_cap_static = Qxx_static*A_static'*p*delta_l;        % for static mode
dx_cap_kinematic = Qxx_kinematic*A_kinematic'*p*delta_l;   % for kinematic mode
% residuals
v_cap_static = A_static*dx_cap_static-delta_l;           % for static mode
v_cap_kinematic = A_kinematic*dx_cap_kinematic-delta_l;  % for kinematic mode
% Aposteriori variance factor
Apos_var_static = v_cap_static'*p*v_cap_static/(size(A_static,1)-size(A_static,2));     % for static mode
Apos_var_kinematic = v_cap_kinematic'*p*v_cap_kinematic/(size(A_kinematic,1)-size(A_kinematic,2));   % for kinematic mode
% covariance matrix 
C_xx_Static = Apos_var_static*Qxx_static;            % for static mode
C_xx_kinematic = Apos_var_kinematic*Qxx_kinematic;   % for kinematic mode
% Coordinate solution
% for static mode
X_static = r_pos(1,1)+dx_cap_static(sat+5);
Y_static = r_pos(1,2)+dx_cap_static(sat+6);
Z_static = r_pos(1,3)+dx_cap_static(sat+7);
XYZ_static =[X_static Y_static Z_static];
% for kinematic mode
for i=1:1:480
X_kinematic(i)= r_pos(i,1)+dx_cap_kinematic(21+4*i-3);
Y_kinematic(i)= r_pos(i,2)+dx_cap_kinematic(21+4*i-2);
Z_kinematic(i)= r_pos(i,3)+dx_cap_kinematic(21+4*i-1);
end
XYZ_kinematic =[X_kinematic', Y_kinematic', Z_kinematic'];
%% task 1.e 
% dividing the design matrix into global and epoch wise parameters for
% better computation using
%% pre elimination and back substitution approach
A1_static = A_static(:,1:24);  % global parameters for static mode
A2_static = A_static(:,25:end); % epoch parameters for static mode
A1_kinematic = A_kinematic(:,1:21); %global parameters for kinematic mode
A2_kinematic = A_kinematic(:,(22:end)); % epoch parameters for kinematic mode
% computation of system of normal equations for static and kinematic mode
N11_static = A1_static'*p*A1_static;
N12_static = A1_static'*p*A2_static;
N11_kinematic=(A1_kinematic'*p*A1_kinematic)';
N21_kinematic=A2_kinematic'*p*A1_kinematic;
N21_static = N12_static';
N22_static = A2_static'*p*A2_static;
N12_kinematic = A1_kinematic'*p*A2_kinematic;
N22_kinematic = A2_kinematic'*p*A2_kinematic;
n1_static = A1_static'*p*delta_l;
n2_static = A2_static'*p*delta_l;
n1_kinematic = A1_kinematic'*p*delta_l;
n2_kinematic = A2_kinematic'*p*delta_l;
% pre elimination
N11_bar_static = N11_static-N12_static*inv(N22_static)*N21_static;
n1_bar_static = n1_static-N12_static*inv(N22_static)*n2_static;
N11_bar_kinematic = N11_kinematic-N12_kinematic*inv(N22_kinematic)*N21_kinematic;
n1_bar_kinematic = n1_kinematic-N12_kinematic*inv(N22_kinematic)*n2_kinematic;
Q11_static = inv(N11_bar_static);
Q11_kinematic = inv(N11_bar_kinematic);
x1_cap_static =Q11_static * n1_bar_static;
x1_cap_kinematic = Q11_kinematic * n1_bar_kinematic;
% %Back substitution
x2_cap_static = inv(N22_static) * ( n2_static - N21_static * x1_cap_static);
x2_cap_kinematic = inv(N22_kinematic) * ( n2_kinematic - N21_kinematic * x1_cap_kinematic);
v_static = A1_static*x1_cap_static + A2_static*x2_cap_static - delta_l;
v_kinematic = A1_kinematic*x1_cap_kinematic + A2_kinematic*x2_cap_kinematic - delta_l;
Q22_static  = inv(N22_static) + inv(N22_static) * N12_static' * Q11_static * N12_static * inv(N22_static);
Q22_kinematic =inv(N22_kinematic) + inv(N22_kinematic) * N12_kinematic' * Q11_kinematic * N12_kinematic * inv(N22_kinematic);
% coordinate solution
X_PE=XYZ(1,1)+x1_cap_static(22);
Y_PE=XYZ(1,2)+x1_cap_static(23);
Z_PE=XYZ(1,3)+x1_cap_static(24);
XYZ_PE =[X_PE,Y_PE,Z_PE];
for i=1:480
    X_BS(i) = XYZ(i,1) + x2_cap_kinematic(3*i-2);
    Y_BS(i) = XYZ(i,2) + x2_cap_kinematic(3*i-1);
    Z_BS(i) = XYZ(i,3) + x2_cap_kinematic(3*i-0);
end
XYZ_BS =[X_BS; Y_BS; Z_BS]';
%% Task 2
%% Transformation of coordinate solution from ECEF coordinate system to topocentric coordinate system
% and plotting the coordinate solution using the reference coordinates
% Reference coordinates
Xr = rec(1);   
Yr = rec(2);    
Zr = rec(3); 
% phi = 52° 17' 47.25"
% lambda = 10° 27' 50.2"
% converting dms to degrees
% phi    = dms2degrees([52 17 47.25]);
% lambda =  dms2degrees([10 27 50.2]);
wgs84 = wgs84Ellipsoid('meter');
[phi,lambda,h] = ecef2geodetic(wgs84,Xr,Yr,Zr);

% rotation matrix
Rot = [-sin(phi)*cos(lambda) -sin(phi)*sin(lambda) cos(phi);
     -sin(lambda)                 cos(lambda)           0;
     cos(phi)*cos(lambda)  cosd(phi)*sin(lambda)  sin(phi)];
% transformation by multiplying rotation matrix
tr = Rot * XYZ_kinematic';
% coordinate solution after transformation 
N_static = Rot * [X_static-Xr;Y_static-Yr;Z_static-Zr];
N_kinematic = Rot * [X_kinematic-Xr;Y_kinematic-Yr;Z_kinematic-Zr];
% calculating mean, standard deviation and root mean square error 
% of each component
% mean
mean_x = mean(N_kinematic(1,:))*100;
mean_y = mean(N_kinematic(2,:))*100;
mean_z = mean(N_kinematic(3,:))*100;
% standard deviation    
s_x  = std(N_kinematic(1,:))*100;
s_y  = std(N_kinematic(2,:))*100;
s_z  = std(N_kinematic(3,:))*100;
% root mean square
rms_x  = sqrt(mean(N_kinematic(1,:))^2 + std(N_kinematic(1,:))^2)*100;  
rms_y  = sqrt(mean(N_kinematic(2,:))^2 + std(N_kinematic(2,:))^2)*100;       
rms_z  = sqrt(mean(N_kinematic(3,:))^2 + std(N_kinematic(3,:))^2)*100;
% plots 
figure (3)
plot(Epochs, N_kinematic(1,:), 'b')
hold on
plot(Epochs, N_kinematic(2,:), 'k')
plot(Epochs, N_kinematic(3,:), 'm')
title('Coordinate solution in topocentric frame vs time')
xlabel('epoch');
ylabel('N,E,U [m]');
legend('N: mean = 1.0719 std = 0.6288 rms = 1.2427 ' ,  'E: mean = -0.0198 std = 1.0053 rms = 1.0055' ,'U: mean = -0.2381 std = 1.2793 rms = 1.3013');
figure (4)
plot(N_kinematic(1,:), N_kinematic(2,:), 'k*')
hold on
plot(N_static(1),N_static(2),'ro')
title('2D Scatter figure');
xlabel('East [m]');
ylabel('North [m]');
legend('Kinematic','Static');
%% Task 3
%% Formal standard deviation of coordinate solutiion in topocentric frame 
% in static and kinematic mode
% static mode
% using cofactor matrix we will find the standard deviation
Q_topo_static = Rot * Qxx_static(sat+5:sat+7 , sat+5:sat+7) * Rot';
NN_static = sqrt(Q_topo_static(1,1));
EE_static = sqrt(Q_topo_static(2,2));
UU_static = sqrt(Q_topo_static(3,3));
% for kinematic mode
NN_kinematic = [];
EE_kinematic = [];
UU_kinematic = [];
for i=22:4:4*480+21
    Q_ECEF_kinematic = Qxx_kinematic(i:i+2,i:i+2);     % cofactor matrix in ECEF reference frame
    Q_t_kinematic = Rot * Q_ECEF_kinematic * Rot';           % cofactor matrix in topocentric reference frame
    NN_kinematic_r = sqrt(Q_t_kinematic(1,1));
    EE_kinematic_r = sqrt(Q_t_kinematic(2,2));
    UU_kinematic_r = sqrt(Q_t_kinematic(3,3));
    NN_kinematic = [NN_kinematic ; NN_kinematic_r];
    EE_kinematic = [EE_kinematic ; EE_kinematic_r];
    UU_kinematic = [UU_kinematic ; UU_kinematic_r];
end
NEU_kin = [NN_kinematic EE_kinematic UU_kinematic];
%plot
figure (5)
plot( NEU_kin(:,1),'-r')
hold on
plot( NEU_kin(:,2),'-g')
hold on
plot( NEU_kin(:,3),'-b')
title('Formal Standard Deviation'),
xlabel('Epochs (30s)'),
ylabel('std(cm)'),
legend('N', 'E', 'U');
hold off;
%% Task 4
%% Plotting of other estimated parameters
% receiver clock error
figure (6)
plot(dx_cap_static(25:503)/(3*10^8)*10^6)
title('Receiver Clock error'),
xlabel('Epochs'),
ylabel('Clock [micro sec]'),
grid on;
%plotting trophospheric wet delay
ZWD_static=dx_cap_static(18:21);
ZWD_static_plot = [ZWD_static(1)*ones(120), ZWD_static(2)*ones(120), ZWD_static(3)*ones(120), ZWD_static(4)*ones(120)];
ZWD_kinematic=dx_cap_kinematic(18:21);
ZWD_kinematic_plot = [ZWD_kinematic(1)*ones(2*60), ZWD_kinematic(2)*ones(2*60), ZWD_kinematic(3)*ones(2*60), ZWD_kinematic(4)*ones(2*60)];
figure (7)
plot(Epochs, ZWD_static_plot);
hold on
plot(Epochs, ZWD_kinematic_plot);
title('Zenith Wet Delay(Static and Kinematic)');
xlabel('Epochs');
ylabel('ZTD [m]');
grid on ;
%% Task 5
%% plot the code and phase residuals with respect to satellite elevation
code_res_kinematic=zeros(480,32);
phase_res_kinematic=zeros(480,32);
k=0;
for i=1:480
    index=find((C_IFCOMC(i,:)~=0));
    num=size(index,2);
    code_res_kinematic(i,index) = v_cap_kinematic(k+1:k+num);      
    phase_res_kinemaatic(i,index) = v_cap_kinematic(k+num+1:k+num*2);  
    k=k+2*num;
end
figure(8)
plot(El/pi*180,code_res_kinematic)
xlim([10 90])
title('code residuals wrt elevation')
xlabel('Elevation(in degrees)')
ylabel('kinematic-code(m)')
figure(9)
plot(El/pi*180,phase_res_kinemaatic)
xlim([10 90])
title('phase residuals wrt elevation')
xlabel('Elevation(in degrees)')
ylabel('kinematic-phase(m)')

%% Task 6
%% Using Identical Standard deviation for code and phase observations
% code obs with identical standard deviation
C_isd=[];
for i=1:480
    I = El(i,:);
    I=(I(find(I~=0)));
    C_isd = [C_isd (sin(I))/Sc^2 (sin(I))/Sc^2];
end
 C_isd((C_isd)==0)=[];
 % weight matrix
 p_isd=diag(C_isd);
 % cofactor matrix
 Q_isd=inv(A_kinematic'*p_isd*A_kinematic);
dx_cap_kinematic_isd = Q_isd * A_kinematic'* p_isd * delta_l;
for i=1:480
    X_kinematic_isd(i)=XYZ(1,1)+dx_cap_kinematic_isd(4*i-3);
    Y_kinematic_isd(i)=XYZ(1,2)+dx_cap_kinematic_isd(4*i-2);
    Z_kinematic_isd(i)=XYZ(1,3)+dx_cap_kinematic_isd(4*i-1);
end
% coordinate solution with identical standard deviation 
N_isd = Rot * [X_kinematic_isd-Xr;Y_kinematic_isd-Yr;Z_kinematic_isd-Zr];
% mean
m_N_isd = mean(N_isd(1,:));
m_E_isd = mean(N_isd(2,:));
m_U_isd = mean(N_isd(3,:));     
% std
std_N_isd  = std(N_isd(1,:));
std_E_isd  = std(N_isd(2,:));
std_U_isd  = std(N_isd(3,:)); 
% rms
rms_N_isd  = sqrt(mean(N_isd(1,:))^2 + std(N_isd(1,:))^2); 
rms_E_isd  = sqrt(mean(N_isd(2,:))^2 + std(N_isd(2,:))^2);
rms_U_isd  = sqrt(mean(N_isd(3,:))^2 + std(N_isd(3,:))^2);     
figure (10)
plot(Epochs, N_kinematic(1,:), 'r')
hold on
plot(Epochs, N_kinematic(2,:), 'g')
hold on
plot(Epochs, N_kinematic(3,:), 'b')
grid on
title("Topocentric coordinates using Identical standard deviation in Kinematic mode");
xlabel("Epochs(30s)");
ylabel("Coords (m)");
legend('N', 'E','U');
hold off