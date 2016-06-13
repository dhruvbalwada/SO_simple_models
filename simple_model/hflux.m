function [surf_flux] = surf_flux(Tm)

    global f Lx Ly H hm cpa rhoa rhow cd a b c eps tau0 K k sbc Ah y z tau Tair uw Pa

    Hir = 4*sbc*(Tm - Tair).*(Tair+273.15).^3;
    Hsens = rhoa*Ah*uw*cpa.*(Tm-Tair);
    L = 2500800-2300*Tm;
    
    Hlat = 0.0015*L*rhoa.*uw*(rm_vap(Tm) - rm_vap(Tair)); 
    
    E = Hlat./L;
end

function [rm] = rm_vap(T)
    global a b c Pa eps
    
    ebar = 0.98*10^((a+b*T)./(1+c*T));
    rm = eps*ebar/Pa./(1-ebar/Pa*(1-eps));

end