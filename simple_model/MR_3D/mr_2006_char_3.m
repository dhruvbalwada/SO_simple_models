% Residual mean solution for the ACC and its assosciated overturning
% circulation
% The prognostic model in which the mixed layer buoyancy is relaxed to some
% atmospheric conditions
% see section 5 in Radko and Marshall 2006 for more details

clear all
close all
matlab_flag = 2014;
%% Find the zeroth order solution, characteristic lines first
% make some sort of a grid

%% define conditions
% Length scales and other physical parameters
Ly = 2000*1000;
Lx = 21000*1000;
Lz = 3000;
f = -10^(-4);
beta = 10^(-11);
ko = 10^6;           % eddy parameter
kx = 2*pi/Lx;         % wave number

dz = 5;
dy =100*1000;
dx = 100*1000;
hm = 200;          % mixed layer depth
z = -[hm:dz:Lz];
y = [0:dy:Ly];
x = [0:dx:Lx];
% surface wind stress
tau_o = 2*10^(-4);
tau0 = tau_o*(0.6 + sin(pi*y/Ly));
tau1 = 0.3*exp(0.5*pi*sqrt(-1))*tau0;
tau = zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        tau(i,j) = tau0(j) + real(tau1(j)*exp(sqrt(-1)*kx*x(i)));
    end
end
figure
contourf(x,y,tau')
title('Wind stress')

% target buoyancy at the surface
delbo = 0.015;
b0star = delbo*y/Ly;
b1star = -0.07*exp(-0.3*pi*sqrt(-1))*delbo*sin(pi*y/Ly);

bstar = zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        bstar(i,j) = b0star(j) + real(b1star(j)*exp(sqrt(-1)*kx*x(i)));
    end
end
figure
contourf(x,y,bstar')
title('Prescribed buoyancy contours')


%% Algorithm
% 1. Compute 2D zero order solution
% 2. Initial guess is made for b0m (start with linear)
% 3. compute the buoyancy flux B and integrate method of characteristics
% 4. compare solution at y=Ly(north) and the difference is used to adjust
% the surface value. repeat till northern value converges

% Implement the above algo

dl =75000;
l=[0:dl:3000*1000];
chary = zeros(length(l),length(y));
charz = zeros(length(l),length(y));
id_hm = find(z==-hm);
lambda = 9*10^(-6);
he = 1000;
for iter =1:100%:10
    if iter ==1
        % guess at a surface buoyancy
        b0m = b0star;
        B0 = -lambda*(b0m - b0star);
        psi_res_0(1) = B0(1)/((b0m(2)-b0m(1))/(y(2)-y(1)));
        for k =2:length(y)-1
            psi_res_0(k) = B0(k)/((b0m(k+1)-b0m(k-1))/(y(k+1)-y(k-1)));
        end
        psi_res_0(length(y)) = B0(length(y))/((b0m(length(y))-b0m(length(y)-1))/(y(end)-y(end-1)));
    else
        % The northern boundary condition
        idy = find(chary(:,1)<=2000*1000,1,'last');
        hthermo = charz(idy,1);
        A1 = delbo/(exp(-hm/he) - exp(hthermo/he));
        A2 = -A1*exp(hthermo/he);
        
        % update the surface value to the corresponding required value
        max_berr = 0;
        for k =1:length(y)
            idy = find(chary(:,k)<=2000*1000,1,'last');
            if ~isempty(idy)
                bN = A1*exp(charz(idy,k)/he)+A2;
                berr = b0m(k) -  bN;
                b0m(k) = b0m(k) - 0.01*(berr); % did this slow relaxation factor of 0.01 (don't know exactly why but worked)
                max_berr = max(max_berr, berr);
            end
        end
%         disp(max_berr);
        % relcalculate fluxes etc
        B0 = -lambda*(b0m - b0star);
        psi_res_0(1) = B0(1)/((b0m(2)-b0m(1))/(y(2)-y(1)));
        for k =2:length(y)-1
            psi_res_0(k) = B0(k)/((b0m(k+1)-b0m(k-1))/(y(k+1)-y(k-1)));
        end
        psi_res_0(length(y)) = B0(length(y))/((b0m(length(y))-b0m(length(y)-1))/(y(end)-y(end-1)));
    end
    % Solve for char curves for solution to first order of the equations
    for j=1:length(y)
        chary(1,j) = y(j);
        charz(1,j) = z(1);
        for i =2:length(l)
            chary(i,j)=chary(i-1,j)+1.0*dl;
            idy = find(y<=chary(i,j),1,'last');
            charz(i,j)=charz(i-1,j) - dl*sqrt(-tau0(idy)/f/ko - psi_res_0(j)/ko);
        end
    end
end
charz = real(charz);
b0m = real(b0m);
psi_res_0 = real(psi_res_0);
B0 = real(B0);
%% Plot the characteristics
close all
figure
plot(chary,charz,'-o','linewidth',1)
axis([min(y) max(y) -3000 0])
title('Characteristics')

figure
plot(y,B0,'o-')
title('Surface buoyancy flux')

figure
plot(y,b0m)
title('Surface buoyancy')
%%
% charz(2:end+1,:) = charz;
% charz(1,:) = 0;
% chary(2:end+1,:) = chary;
%%
b0grid = 0*charz;
for i =1:length(b0m);
    b0grid(:,i) = b0m(i);
end
if matlab_flag==2014
    b0grid_interp = scatteredInterpolant(chary(:),charz(:),b0grid(:));
else
    b0grid_interp = TriScatteredInterp(chary(:),charz(:),b0grid(:));
end

%%
[zi yi] = meshgrid([-3000:0], y);
vq = b0grid_interp(yi,zi);

figure
contourf(yi,zi,vq,[0:1:15]*10^-3)
axis([0 2*10^6 -3000 0])
%%
db0dy = nan*chary;
db0dz = nan*chary;
% calculate db0/dy
for j=1:length(b0m)
    for i=1:length(l)
        zlev= charz(i,j);
        ylev= chary(i,j);
        yi = [ylev-dy/2 ylev+dy/2];
        
        db0dy(i,j) = (b0grid_interp(yi(2),zlev)-b0grid_interp(yi(1),zlev))/(yi(2)-yi(1));
        
        zi = [zlev-dz/2 zlev+dz/2];
        db0dz(i,j) = (b0grid_interp(ylev, zi(2))-b0grid_interp(ylev, zi(1)))/(zi(2)-zi(1));
        
    end
end

% break

%% calculate the zonal velocity
if matlab_flag==2014
    db0dy_interp = scatteredInterpolant(chary(:),charz(:),db0dy(:));
else
    db0dy_interp = TriScatteredInterp(chary(:),charz(:),db0dy(:));
end

u0 = 0*chary;

for j=1:length(b0m)
    for i =1:length(l)
        ylev = chary(i,j) ;
        zlev = charz(i,j) ;
        idy = find(chary(:,1)<=ylev,1,'last');
        h = charz(idy,1); % find thermocline depth
        zint = [h:dz:zlev];
        yint = ylev*ones(1,length(zint));
        db0dyint = db0dy_interp(yint,zint);
        zint = [zint zlev+dz];
        
        u0(i,j) = -1/f*nansum(db0dyint.*diff(zint));
    end
end

%%
figure
% contourf(chary, real(charz), abs(u0),[0:0.4:5]*10^(-2))
contourf(chary, real(charz), abs(u0),[0:7]*10^(-2))
axis([min(y) max(y) -3000 0])
caxis([0 5]*10^-2)
title('Baroclinic flow')

%break

%% calculate slopes

s0 = -db0dy./db0dz;

%
figure
contourf(real(chary), real(charz), abs(ko*real(s0)),[0:100:2500])
axis([min(y) max(y) -3000 0])
caxis([0 1500])

% break
%% Calculate  order eps Psi_res at the base of the boundary layer
options = optimoptions('fsolve','Display','iter','Algorithm','levenberg-marquardt');
b1_ml = fsolve(@(b1m) solve_for_b1(b1m, b1star, chary, charz,psi_res_0,l, y,dy,dz, f, ko,...
    Lx, b0m,db0dy, db0dz, s0,u0, tau1, lambda, hm), b1star, options);

