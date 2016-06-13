% This script reproduces the results of the paper Badin and Williams 2010
% A model which incorporates from Mixed layer physics to compare its
% interactions with the overturning circulation. Slight improvement over
% the Marshall and Radko 2003 work, which assumed a buoyancy relaxation
% condition at the surface.

% parameters
close all
clear all

global f Lx Ly H hm cpa rhoa rhow cd a b c eps tau0 K k sbc Ah y z tau Tair uw Pa

f = -10^(-4);
Lx = 21*10^6;
Ly = 2*10^6;
H = 1000;
hm = 200;          % mixed layer depth
cpa = 1008;
cpw = 4190;
rhoa = 1.25;
rhow = 1027;
cd = 1.3;
a = 0.7859; b= 0.03477; c= 0.00412;
eps = 0.62197;
tau0 = 0.15;
K = 500;
k = Ly*K/H;
sbc = 5.670367 *10^(-8);
Ah = 10^(-3); % Take to be a constant but varies (Can be obtained from Table 2 of Smith 1988)
Pa = 10^5;
alpha = 1*10^-4;
beta = 8*10^-4;
T0 = 0;
S0 = 34;
rho0 = 1026;
g = 9.81;

% define the grid
dz = 5;
dy =25*1000;

z = -[hm:dz:2500];
y = [0:dy:2000*1000];

% forcing conditions
tau = tau0*sin(pi*y/Ly);
Tair = (12-1)*y.^2/Ly^2 +1;
Simp = (35-34)*y/Ly + 34;
uw = sqrt(tau*rhow/rhoa/cd);
% end

%% Get initial setup

% insert the initial profile
T_init = Tair;
S_init = Simp;
drho_init = -alpha*rho0*(T_init-T0) + beta*rho0*(S_init - S0);
Pinit = 0;
surf_flux_init = hflux(T_init);
Hinit = surf_flux_init.Hir + surf_flux_init.Hsens + surf_flux_init.Hlat;
EPinit=(surf_flux_init.E - Pinit);
Dinit = -alpha/cpw*Hinit + beta*rho0*S_init.*EPinit;
psi_res_ml_init = Dinit./(drho_init/dy);
b_ml_init = -drho_init*g/rho0;

% solve the characteristic equations (22 in MR 2003)
% dy/dl = 1; dz/dl = wc; where l is distance along chanracteristic
dl =dy;
l=[0:dl:Ly+2000];
chary = zeros(length(l),length(b_ml_init));
charz = zeros(length(l),length(b_ml_init));
for j=1:length(b_ml_init)
    chary(1,j) = y(j);
    charz(1,j) = z(1);
    for i =2:length(l)
        chary(i,j)=chary(i-1,j)+dl;
        idy = find(y<=chary(i,j),1,'last');
        charz(i,j)=charz(i-1,j) - dl*sqrt(-tau(idy)/f/k/rho0 - psi_res_ml_init(j)/k);
    end
end

% initial slope at the base
for j=1:length(b_ml_init)
    slope_init(j) = diff(charz(1:2,j))/diff(chary(1:2,j));
end

%% Make an iterative solver for mixed layer T,S equations (13, 14 in bd 2010)

dt= 10^3; % choose some time step
for k =1:1000
    if k==1
        T_old = T_init;
        S_old = S_init;
        H_old = Hinit;
        EP_old = EPinit;
        psi = psi_res_ml_init;
        slope = slope_init;
    end
    
    gradT = gradients(T_old,dy);
    gradS = gradients(S_old,dy);
    
    T_new = T_old + dt*(H_old/rho0/cpw/hm + k*slope.*gradT.d2 - psi./hm.*gradT.d1);
    S_new = S_old + dt*(EP_old.*S_old/hm + k*slope.*gradS.d2 - psi./hm.*gradS.d1);
    
    errT = mean((T_new-T_old).^2);
    errS = mean((S_new-S_old).^2);
%     
%     if errT<0.01 & errS<0.001
%         break
%     end
    
    T_old = T_new;
    S_old = S_new;
    
    drho_old = -alpha*rho0*(T_old-T0) + beta*rho0*(S_old - S0);
    P_old = 0;
    surf_flux_old = hflux(T_old);
    H_old = surf_flux_old.Hir + surf_flux_old.Hsens + surf_flux_old.Hlat;
    EP_old= (surf_flux_old.E - P_old);
    D_old = -alpha/cpw*H_old + beta*rho0*S_old.*EP_old;
    psi_res_ml_old = D_old./(drho_old/dy);
    b_ml_old = -drho_old*g/rho0;
    
    for j=1:length(b_ml_old)
        chary(1,j) = y(j);
        charz(1,j) = z(1);
        for i =2:length(l)
            chary(i,j)=chary(i-1,j)+dl;
            idy = find(y<=chary(i,j),1,'last');
            charz(i,j)=charz(i-1,j) - dl*sqrt(-tau(idy)/f/k/rho0 - psi_res_ml_old(j)/k);
        end
    end
    
    % initial slope at the base
    for j=1:length(b_ml_old)
        slope(j) = diff(charz(1:2,j))/diff(chary(1:2,j));
    end
    psi = psi_res_ml_old;
    
    figure(1)
    hold all
    plot(k, errT,'o')
end


%% Plot the characteristics
plot(chary,charz,'-o')
axis([min(y) max(y) -3000 0])
