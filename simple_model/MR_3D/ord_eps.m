function [df] = ord_eps(t,ysol, char_index,chary,  dpsi_res_0db, l , y, f, ko, k, db0dz, s0, u0, tau1, RHS1, RHS2 )
df= zeros(2,1) ;

yinter = interp1(l, chary(:, char_index),t);
a = interp1(y, tau1, yinter);
if isnan(a) 
    a=0;
end
b = interp1(l, db0dz(:,char_index), t);
c = interp1(l, s0(:,char_index), t);
d = interp1(l, u0(:,char_index), t);
e = interp1(l, RHS1(:,char_index),t);
g = interp1(l, RHS2(:,char_index),t);


% if abs(c)>=10^-5
    df(1) = (ysol(2) + a/f)*b/(2*ko*c)+e ;% need to add the right hand side term here
    df(2) = dpsi_res_0db*df(1) - sqrt(-1)*k*d*ysol(1)/b+g;
% else
%     df(1) = 0 
%     df(2) = dpsi_res_0db*df(1) - sqrt(-1)*k*d*ysol(1)/b+g;
% end


    
if isnan(df(1))
    disp('why')
%     exit
end

end