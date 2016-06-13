% Numerical solutions to Marshall and Radko 2003
% Residual mean solution for the ACC and its assosciated overturning
% circulation
% most is for the diagnostic model and not the prognostic one

clear all
close all

matlab_flag =2014;
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
k = 2*pi/Lx;         % wave number

dz = 5;
dy =50*1000;
dx = 100*1000;
hm = 200;          % mixed layer depth
z = -[hm:dz:Lz];
y = [0:dy:Ly];
x = [0:dx:Lx];
% surface wind stress
tau_o = 1*10^(-4);
tau0 = tau_o*(0.6 + sin(pi*y/Ly));
tau1 = 0.3*exp(0.5*pi*sqrt(-1))*tau0;
tau = zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        tau(i,j) = tau0(j) + real(tau1(j)*exp(sqrt(-1)*k*x(i)));
    end
end
figure
contourf(x,y,tau')
title('Wind stress')

% buoyancy change across ACC
delbo = 0.015;
b0 = delbo*y/Ly;
b1 = -0.07*exp(-0.3*pi*sqrt(-1))*delbo*sin(pi*y/Ly);
b1m= b1;

bm = zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        bm(i,j) = b0(j) + real(b1(j)*exp(sqrt(-1)*k*x(i)));
    end
end
figure
contourf(x,y,bm')
title('Surface buoyancy contours')


% Surface byouancy forcing
Bo = 3*10^(-9);
B0 = Bo*sin(pi*y/Ly);
B1 = -0.5*(1-sqrt(-1))*B0;
for i=1:length(x)
    for j=1:length(y)
        B(i,j) = B0(j) + real(B1(j)*exp(sqrt(-1)*k*x(i)));
    end
end
figure
contourf(x,y,B')
title('Surface buoyancy flux')


%% find the conditions at the bottom of the homogeneous, steady mixed layer

id_hm = find(z==-hm);
% psi_res_0(:,id_hm) = B0/(delbo/Ly); % choose this because the gradient is a constant here
psi_res_0 = B0/(delbo/Ly);

figure 
plot(y,psi_res_0.*Lx/10^6)
% psi_res_1 = (B1 - psi_res_0*

%% Solve for char curves for solution to first order of the equations
dl =50000;
l=[0:dl:Ly];
chary = zeros(length(l),length(b0));
charz = zeros(length(l),length(b0));


for j=1:length(b0)
    chary(1,j) = y(j);
    charz(1,j) = z(1);
    for i =2:length(l)
        chary(i,j)=chary(i-1,j)+1.0*dl;
        idy = find(y<=chary(i,j),1,'last');
        charz(i,j)=charz(i-1,j) - dl*sqrt(-tau0(idy)/f/ko - psi_res_0(j)/ko);
    end
end

%%
charz= real(charz);

%% Plot the characteristics
plot(chary,charz,'-o')
axis([min(y) max(y) -3000 0])


%% Calculate the baroclinic zonal transport
% assumes that the flow below the thermocline is zero
% thermocline is defined as the characteristic that starts at y=0

b0grid = 0*charz;
for i =1:length(b0);
    b0grid(:,i) = b0(i);
end
if matlab_flag == 2014
    b0grid_interp = scatteredInterpolant(chary(:),charz(:),b0grid(:));
else
    b0grid_interp = TriScatteredInterp(chary(:),charz(:),b0grid(:));
end
% [zi yi] = meshgrid([-2000:-200], y);
% vq = b0grid_interp(yi,zi);

%%
db0dy = nan*chary;
db0dz = zeros*chary;
% calculate db0/dy
for j=1:length(b0)
    for i=1:length(l)
        zlev= charz(i,j);
        ylev= chary(i,j);
        if ylev == y(1)
                    yi = [ylev ylev+dy];
                elseif ylev == y(end)
                    yi = [ylev-dy ylev];
                else
                    yi = [ylev-dy/2 ylev+dy/2];
        end
        db0dy(i,j) = (b0grid_interp(yi(2),zlev)-b0grid_interp(yi(1),zlev))/(yi(2)-yi(1));
%         zi = [zlev-dz/2 zlev+dz/2];
             if zlev == z(1)
                    zi = [zlev zlev+dz];
                elseif zlev == z(end)
                    zi = [zlev-dz zlev];
                else
                    zi = [zlev-dz/2 zlev+dz/2];
                end            
    db0dz(i,j) = (b0grid_interp(ylev, zi(2))-b0grid_interp(ylev, zi(1)))/(zi(2)-zi(1));
        
    end
end

if matlab_flag==2012
    db0dy(isnan(db0dy))= 0;
    db0dz(isnan(db0dz))= 0;
end

%%
figure
contourf(real(chary), real(charz), db0dy)
axis([min(y) max(y) -3000 0])
%%
figure
contourf(real(chary), real(charz), db0dz)
axis([min(y) max(y) -3000 0])



%% calculate the zonal velocity
if matlab_flag == 2014
    db0dy_interp = scatteredInterpolant(chary(:),charz(:),db0dy(:));
else
    db0dy_interp = TriScatteredInterp(chary(:),charz(:),db0dy(:));
end

u0 = 0*chary;

for j=1:length(b0)
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
contourf(chary, real(charz), u0,[0:10]*10^(-2))
axis([min(y) max(y) -3000 0])

% break
%% Calculate  order eps Psi_res at the base of the boundary layer
b1= b1m;
clear db1mdy
db1mdy(1) = (b1(2)-b1(1))/dy;
for i=2:length(b1)-1
    db1mdy(i) = (b1(i+1)-b1(i-1))/2/dy;
end
db1mdy(length(b1)) = (b1(end)-b1(end-1))/dy;


psi_res_1_ml = (B1 - psi_res_0.*db1mdy - sqrt(-1)*k*squeeze(u0(1,:)).*b1m*hm)./db0dy(1,:);
if matlab_flag == 2012
    psi_res_1_ml(isnan(psi_res_1_ml)) = 0;
    psi_res_1_ml(isinf(psi_res_1_ml)) = 0;
end
%% calculate slopes

s0 = -db0dy./db0dz;

figure
contourf(real(chary), real(charz), abs(ko*(s0)),[0:100:2000],'showtext','on')
axis([min(y) max(y) -3000 0])
caxis([0 2500])
if matlab_flag==2012
    s0(isnan(s0)) = 0;
    s0(isinf(s0)) = 0;
end
% caxis([0 1500])
% break
%% Solve coupled ODEs of equation 21
% this is the iterative solution to the first order
b1 = 0*charz;
u1 = 0*charz;
psi_res_1 = 0*charz;
RHS1 = 0*charz;
RHS2 = 0*charz;
RHS1old = 0*charz;
RHS2old = 0*charz;
tol1 = 0;
tol2 = 0;

% psi_res_1_ml(:) = 0;
for iter = 1:10
    % calculate the RHS for equation 21 in iterative sense
    
    if iter ==1
        RHS1old(:,:) =1;
        RHS2old(:,:) =1;
        RHS1(:,:) = 0;
        RHS2(:,:) = 0;
    else
        RHS1old = RHS1;
        RHS2old = RHS2;
        
        % new
        if matlab_flag == 2014
            b1interp = scatteredInterpolant(chary(:),charz(:),b1(:));
        else
            b1interp = TriScatteredInterp(chary(:),charz(:),b1(:));
        end
        P1 =nan*b1;
        u1 = nan*b1;
        
        % calculate P1
        for j =1:length(b0)
            for i = 1:length(l)
                ylev = chary(i,j) ;
                zlev = charz(i,j) ;
                idy = find(chary(:,1)<=ylev,1,'last');
                h = charz(idy,1); % find thermocline depth
                zint = [h:dz:zlev];
                yint = ylev*ones(1,length(zint));
                b1int = b1interp(yint,zint);
                zint = [zint zlev+dz];
                P1(i,j) = nansum(b1int.*diff(zint));
            end
        end
        
        if matlab_flag == 2014
            P1interp = scatteredInterpolant(chary(:), charz(:), P1(:));
        else
            P1interp = TriScatteredInterp(chary(:), charz(:), P1(:));
        end
        % calculate u1
        for j=1:length(b0)
            for i=1:length(l)
                zlev= charz(i,j);
                ylev= chary(i,j);
                
                if ylev == y(1)
                    yi = [ylev ylev+dy];
                elseif ylev == y(end)
                    yi = [ylev-dy ylev];
                else
                    yi = [ylev-dy/2 ylev+dy/2];
                end
                
                u1(i,j) = -1/f*(P1interp(yi(2),zlev)-P1interp(yi(1),zlev))/(yi(2)-yi(1));
                
                %                 zlev= charz(i,j);
                %                 idy1 = find(charz(:,j-1)>=zlev,1,'last');
                %                 idy2 = find(charz(:,j+1)>=zlev,1,'last');
                %                 if (~isempty(idy1) && ~isempty(idy2))
                %                     u1(i,j) = -1/f*(P1(idy2,j+1)-P1(idy1,j-1))/(chary(idy2,j+1)-chary(idy1,j-1));
                %                 else
                %                     u1(i,j) = NaN;
                %                 end
            end
        end
        u1interp = scatteredInterpolant(chary(:),charz(:), u1(:));

        % calculate the surface b1 gradient
        for j =1:length(b0)
            zlev= charz(1,j);
            ylev= chary(1,j);
               if ylev == y(1)
                    yi = [ylev ylev+dy];
                elseif ylev == y(end)
                    yi = [ylev-dy ylev];
                else
                    yi = [ylev-dy/2 ylev+dy/2];
                end            
            b1y(j) = -1/f*(b1interp(yi(2),zlev)-b1interp(yi(1),zlev))/(yi(2)-yi(1));
        end
        
%         b1yinterp=scatteredInterpolant(y, 
        % calculate the new RHS
        for j = 1:length(b0)
            for i=1:length(l)
                ylev = chary(i,j);
                zlev = charz(i,j);
                
                %                 clear zint P1int u1int
                zint = []; P1int = []; u1int = []; b1int=[];
%                 zint(1) = charz(i,j);
                zint = [-200:-dz:zlev];
                yint = ylev*ones(1,length(zint));
                P1int = P1interp(yint, zint);
                u1int = u1interp(yint, zint);
                b1int = b1interp(yint, zint);
                zint = [zint(1)+dz zint];
                
                idsurf = find(y<=ylev,1,'last');
%                 n=2;
%                 for p=j+1:length(b0) % this loop runs opposite to other z loops because the integral is from surface
%                     % and not from the thermocline up.
%                     idy = find(chary(:,p)<=ylev,1,'last');
%                     if ~isempty(idy)
%                         zint(n) = charz(idy,p);
%                         %                     idy = find(chary(:,p)>=ylev,1);
%                         P1int(n-1) = P1(idy,p);
%                         u1int(n-1) = u1(idy,p);
%                         b1int(n-1) = b1(idy,p);
%                         n=n+1;
%                     end
%                 end
%                 if length(P1int)>0 
                    RHS1(i,j) = sqrt(-1)*k*db0dz(i,j)*(nansum(P1int.*diff(zint)) + P1int(1)*(-hm) + b1int(1)*(-(hm^2)/2) )/2/ko/s0(i,j)/f;
%                 else
%                     RHS1(i,j) = sqrt(-1)*k*db0dz(i,j)*(P1(i,j)*(-hm)+b1(i,j)*(-(hm^2)/2))/2/ko/s0(i,j)/f;
%                 end
                
%                 if abs(s0(i,j))<=10^-5
%                     RHS1(i,j) = 0;
%                 end
                
%                 if length(P1int)>0
                    RHS2(i,j) = sqrt(-1)*k*(nansum(u1int.*diff(zint)) + u1int(1)*(-hm) + b1y(idsurf)*(-(hm^2)/2));
%                 else
%                     RHS2(i,j) = sqrt(-1)*k*(u1(i,j)*(-hm) + b1y(j)*(-(hm^2)/2));
%                 end
%                 if length(P1int)>0
%                     RHS1(i,j) = sqrt(-1)*k*db0dz(i,j)*(nansum(P1int.*diff(-zint)))/2/ko/s0(i,j)/f;
%                     RHS2(i,j) = sqrt(-1)*k*(nansum(u1int.*diff(-zint)));
%                 else
%                     RHS1(i,j) = 0;
%                     RHS2(i,j) = 0;
%                 end
            end
        end
    end
    
    % calculate the residue
    idy1 = find(l<=200*1000,1,'last');
    idy2 = find(l<=600*1000,1,'last');
    
    %     idz1 = find(z<=200*1000,1,'last');
    %     idz2 = find(z<=1200*1000,1,'last');
    
    tol1(iter) = nanmean(nanmean(abs(RHS1(idy1:idy2,3:10)-RHS1old(idy1:idy2,3:10))))
    tol2(iter) = nanmean(nanmean(abs(RHS2(idy1:idy2,3:10)-RHS2old(idy1:idy2,3:10))));
    
    for char_index=1:length(b0)
        if char_index ==1
            dpsi_res_0db = (psi_res_0(char_index+1)-psi_res_0(char_index))/(b0(char_index+1)-b0(char_index));
        elseif char_index == length(b0)
            dpsi_res_0db = (psi_res_0(char_index)-psi_res_0(char_index-1))/(b0(char_index)-b0(char_index-1));
        else
            dpsi_res_0db = (psi_res_0(char_index+1)-psi_res_0(char_index-1))/(b0(char_index+1)-b0(char_index-1));
        end
        % first guess at solution
        a = warning('off','all');
        id = find(chary(:,char_index)<=2000*1000+dy,1,'last');
        
        [T, Y] = ode23s(@(t, ysol) ord_eps(t,ysol, char_index,chary, dpsi_res_0db, l , y, f, ko,k, db0dz, s0, u0, tau1, RHS1, RHS2) ...
            , l(1:id), [b1m(char_index) psi_res_1_ml(char_index)]);
        
        m = size(Y,1);
        b1(1:m,char_index) = Y(1:m,1);
        psi_res_1(1:m,char_index) = Y(1:m,2);
        disp(char_index)
        
        if length(find(isnan(Y(:,1))))>0
            disp('What ?')
        end
        
        
    end
    figure
    contourf(chary, -abs(charz), abs(b1),[0:30]*10^-4)
    caxis([-1 10]*10^-4)
    axis([0 2000*1000 -2000 0 ])
    drawnow
end


%%
    figure
    contourf(chary, -abs(charz), real(b1),[0:40]*10^-4, 'edgecolor','none')
    hold all
%     plot(chary, -abs(charz))
    caxis([0 50]*10^-4)
    axis([0 2000*1000 -3000 0 ])

%% Calculate the wres
w0(1) =(psi_res_0(2)-psi_res_0(1))/dy;
w1(1) =  -RHS2(1,1) + (psi_res_1_ml(2)-psi_res_1_ml(1))/dy;
for i =2:length(y)-1
    w1(i) = -RHS2(1,i)+ (psi_res_1_ml(i+1)-psi_res_1_ml(i-1))/2/dy;
    w0(i) = (psi_res_0(i+1)-psi_res_0(i-1))/2/dy;
end
w1(length(y)) =  -RHS2(1,length(y))+ (psi_res_1_ml(length(y))-psi_res_1_ml(length(y)-1))/dy;
w0(length(y)) = (psi_res_0(length(y))-psi_res_0(length(y)-1))/dy;

% w0(1) =(psi_res_0(2)-psi_res_0(1))/dy;
% w1(1) =  (psi_res_1_ml(2)-psi_res_1_ml(1))/dy;
% for i =2:length(y)-1
%     w1(i) = (psi_res_1_ml(i+1)-psi_res_1_ml(i-1))/2/dy;
%     w0(i) = (psi_res_0(i+1)-psi_res_0(i-1))/2/dy;
% end
% w1(length(y)) =  (psi_res_1_ml(length(y))-psi_res_1_ml(length(y)-1))/dy;
% w0(length(y)) = (psi_res_0(length(y))-psi_res_0(length(y)-1))/dy;

for i=1:length(x)
    for j=1:length(y)
        wres(i,j) = w0(j) + real(w1(j)*exp(sqrt(-1)*kx*x(i)));
    end
end

% wres = w1res +

wek = 0*wres;
for i=2:length(x)-1
    for j=2:length(y)-1
        wek(i,j) = -1/f*(tau(i,j+1)-tau(i,j-1))/2/dy;
    end
end
% w1res
%%
figure(20), clf
hold all 
contourf(x,y,wres',[-6:0.2:6]*10^-6,'edgecolor','none')
contour(x,y,wres',[0 0],'edgecolor','k')
caxis([-6 6]*10^-6)
% figure
% contourf(x,y,wek)
%
