function [north_err] = solve_for_b1(b1m, b1star, chary, charz,psi_res_0,l, y, dy,dz, f, ko, Lx,b0m, db0dy, db0dz, s0,u0, tau1, lambda, hm)

%% Calculate  order eps Psi_res at the base of the boundary layer
b1          = 0*charz;
u1          = 0*charz;
psi_res_1   = 0*charz;
kx          = 2*pi/Lx;         % wave number


db1mdy(1) = (b1m(2)-b1m(1))/dy;
for i=2:length(b1m)-1
    db1mdy(i) = (b1m(i+1)-b1m(i-1))/2/dy;
end
db1mdy(length(b1m)) = (b1m(end)-b1m(end-1))/dy;
B1 = -lambda*(b1m - b1star);
psi_res_1_ml = (B1 - psi_res_0.*db1mdy - sqrt(-1)*kx*u0(1,:).*b1m*hm)./db0dy(1,:);

% figure
% plot(y,b1m)
% drawnow

% Solve coupled ODEs of equation 21
RHS1    = 0*charz;
RHS2    = 0*charz;
RHS1old = 0*charz;
RHS2old = 0*charz;
tol1    = 0;
tol2    = 0;

for iter = 1:4
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
        b1interp = scatteredInterpolant(chary(:),charz(:),b1(:));
        P1 =nan*b1;
        u1 = nan*b1;
        
        % calculate P1
        for j =1:length(b0m)
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
        P1interp = scatteredInterpolant(chary(:), charz(:), P1(:));
        % calculate u1
        for j=1:length(b0m)
            for i=1:length(l)
                zlev= charz(i,j);
                ylev= chary(i,j);
                yi = [ylev-dy/2 ylev+dy/2];
                u1(i,j) = -1/f*(P1interp(yi(2),zlev)-P1interp(yi(1),zlev))/(yi(2)-yi(1));
            end
        end
        u1interp = scatteredInterpolant(chary(:),charz(:), u1(:));
        for j =1:length(b0m)
            zlev= charz(1,j);
            ylev= chary(1,j);
            yi = [ylev-dy/2 ylev+dy/2];
            
            b1y(j) = -1/f*(b1interp(yi(2),zlev)-b1interp(yi(1),zlev))/(yi(2)-yi(1));
        end
        % calculate the new RHS
        for j = 1:length(b0m)-1
            for i=1:length(l)
                ylev = chary(i,j);
                zlev = charz(i,j);
                %                 clear zint P1int u1int
                zint = []; P1int = []; u1int = [];
                zint(1) = charz(i,j);
                %                     zint = [-200+dz:-dz:zlev];
                %                     yint = ylev*ones(1,length(zint));
                %                     P1int = P1interp(yint, zint);
                %                     u1int = u1interp(yint, zint);
                %                     b1int = b1interp(yint, zint);
                %                     zint = [zint(1)-dz zint];
                n=2;
                for p=j+1:length(b0m) % this loop runs opposite to other z loops because the integral is from surface
                    %                         % and not from the thermocline up.
                    idy = find(chary(:,p)<=ylev,1,'last');
                    if ~isempty(idy)
                        zint(n) = charz(idy,p);
                        %                             %                     idy = find(chary(:,p)>=ylev,1);
                        P1int(n-1) = P1(idy,p);
                        u1int(n-1) = u1(idy,p);
                        b1int(n-1) = b1(idy,p);
                        
                        n=n+1;
                    end
                end
                if length(P1int)>0
                    RHS1(i,j) = sqrt(-1)*kx*db0dz(i,j)*(nansum(P1int.*diff(-zint))+P1int(end)*(-hm)+b1int(end)*(-(hm^2)/2))/2/ko/s0(i,j)/f;
                    RHS2(i,j) = sqrt(-1)*kx*(nansum(u1int.*diff(-zint)) + u1int(end)*(-hm) + b1y(j)*(-(hm^2)/2));
                else
                    RHS1(i,j) = sqrt(-1)*kx*db0dz(i,j)*(P1(i,j)*(-hm)+b1(i,j)*(-(hm^2)/2))/2/ko/s0(i,j)/f;
                    RHS2(i,j) = sqrt(-1)*kx*(u1(i,j)*(-hm) + b1y(j)*(-(hm^2)/2));
                end
                % RHS1(i,j) = sqrt(-1)*kx*db0dz(i,j)*nansum(P1int.*diff(zint))/2/ko/s0(i,j)/f;
                % RHS2(i,j) = sqrt(-1)*kx*nansum(u1int.*diff(zint));
            end
        end
    end
    
    % calculate the residue
    idy1 = find(l<=200*1000,1,'last');
    idy2 = find(l<=600*1000,1,'last');
    
    %     idz1 = find(z<=200*1000,1,'last');
    %     idz2 = find(z<=1200*1000,1,'last');
    
    %         tol1(iter) = nanmean(nanmean(abs(RHS1(idy1:idy2,3:end-25)-RHS1old(idy1:idy2,3:end-25))))
    %         tol2(iter) = nanmean(nanmean(abs(RHS2(idy1:idy2,3:end-25)-RHS2old(idy1:idy2,3:end-25))));
    
    for char_index=1:length(b0m)
        if char_index ==1
            dpsi_res_0db = (psi_res_0(char_index+1)-psi_res_0(char_index))/(b0m(char_index+1)-b0m(char_index));
        elseif char_index == length(b0m)
            dpsi_res_0db = (psi_res_0(char_index)-psi_res_0(char_index-1))/(b0m(char_index)-b0m(char_index-1));
        else
            dpsi_res_0db = (psi_res_0(char_index+1)-psi_res_0(char_index-1))/(b0m(char_index+1)-b0m(char_index-1));
        end
        % first guess at solution
        a = warning('off','all');
        id = find(chary(:,char_index)<=2400*1000,1,'last');
        
        if (id>2 & ~isempty(id))
            [T, Y] = ode23(@(t, ysol) ord_eps(t,ysol, char_index,chary, dpsi_res_0db, l , y, f, ko,kx, db0dz, s0, u0, tau1, RHS1, RHS2) ...
                , l(1:id), [b1m(char_index) psi_res_1_ml(char_index)]);
            
            m = size(Y,1);
            b1(1:m,char_index) = Y(1:m,1);
            psi_res_1(1:m,char_index) = Y(1:m,2);
%             disp(char_index)
        end
    end
    
end
figure
contourf(chary, charz, abs(b1),[0:30]*10^(-4))
caxis([1 30]*10^-4)
axis([0 2000*1000 -2000 0 ])
%      drawnow

b1N = 0;
for k =1:length(y)
    idy = find(chary(:,k) <=1800*1000,1,'last');
    if ~isempty(idy)
        north_err(k) = b1(idy,k) - b1N;
    else
        north_err(k) = 0;
    end
end
figure(2000), hold all 
plot(y,north_err) 
drawnow
