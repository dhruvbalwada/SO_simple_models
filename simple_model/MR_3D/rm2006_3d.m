% Numerical solutions to Marshall and Radko 2003
% Residual mean solution for the ACC and its assosciated overturning
% circulation

clear all
close all

% Find the zeroth order solution first
dz = 5;
dy =50*1000;
hm = 200;          % mixed layer depth

CFL = dz/dy/10^(-4);  % 10^-4 is chosen based on manual entry, proper way would be to calculate CFL in the compute loop

if CFL >1
    disp('CFL too high')
    exit
end

z = -[hm:dz:2500];
y = [0:dy:2000*1000];

tau_o = 2*10^(-4); % surface wind stress
Ly = 2000*1000;     % Meridional scale
Lx = 21000*1000;

tau = tau_o*(0.3 + sin(pi*y/Ly));

delbo = 7*10^(-3); % buoyancy change across ACC
bo = delbo*y/Ly;

Bo = 2*10^(-9);    % Net byouancy forcing
B = -Bo*sin(2*pi*y/Ly);

f = -10^(-4);
beta = 10^(-11);
k = 10^6;           % eddy parameter

psi_bar  = nan(length(y),length(z));
psi_star = nan(length(y),length(z));
psi_res  = nan(length(y),length(z));

id_hm = find(z==-hm);
psi_res(:,id_hm) = B/(delbo/Ly); % choose this because the gradient is a constant here
% otherwise take the proper db/dy
% gradient;

psi_res_b = B/(delbo/Ly);

b=nan(length(y),length(z));
wc= b;
b(:,1) = bo; % specify the surface buoyancy
b(1,:) = bo(1); % southern wall has same bo as at the surface over there

for j=2:length(z)-1 % loop from surface to depth
    idb = find(bo<=b(1,j-1),1,'last'); % find the buoyancy at the
    if ~isempty(idb)
        wc(1,j-1) = -sqrt(-tau(1)/f/k - psi_res_b(idb)/k);
        b(1,j) = b(1,j-1) + (-dz)/wc(1,j-1)*(b(2,j-1)-b(1,j-1))/dy;
    else
        b(1,j)= 0;
    end
    for i =2:length(y)-1
        idb = find(bo<=b(i,j-1),1,'last'); % find the buoyancy at the
        if ~isempty(idb)
            wc(i,j-1) = -sqrt(-tau(i)/f/k - psi_res_b(idb)/k);
            b(i,j) = b(i,j-1) - (-dz)/wc(i,j-1)*(b(i+1,j-1)-b(i-1,j-1))/2/dy;
        else
            b(i,j)=0;
        end
    end
    idb = find(bo<=b(end,j-1),1,'last'); % find the buoyancy at the
    if ~isempty(idb)
        wc(end,j-1) = -sqrt(-tau(end)/f/k - psi_res_b(idb)/k);
        b(end,j) = b(end,j-1) - (-dz)/wc(end,j-1)*(b(end,j-1)-b(end-1,j-1))/dy;
    else
        b(end,j) = 0;
    end
end

%% Find the solutions for epsilon order terms at the bottom of the mixed layer
% make the distance axis
fac=1.2;
l = y*fac;
dl = dy*fac;
clear ly
for k=1:length(bo)
    for i=1:length(z)
        idy = find(abs(b(:,i)-bo(k))==min(abs(b(:,i)-bo(k))),1);
%          idy = find(b(:,i)<=bo(k),1,'last');
        ly(i,k) = y(idy);
    end
end


