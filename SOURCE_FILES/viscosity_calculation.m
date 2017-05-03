
%==============================================================
% Initial effective viscosity
eta_cur         =   etam;

% Power-law viscosity
if Powerlaw==1
    co=(Qm./(RK.*Tm));
    ind=co>150;
    co(ind)=150;
    
    %         eta_ddd     =   ((etam.*exp(co).*E2ndm).^(1./nm))./2./E2ndm;
%     eta_cur     =   2.^((1-nm)./nm).*3.^(-(1+nm)./(2.*nm)).*etam.^(1./nm).*E2ndm.^((1-nm)./nm).*exp(co./nm);    
    eta_cur     =   0.5.*etam.^(1./nm).*E2ndm.^((1-nm)./nm).*exp(co./nm);    
end

ind = eta_cur > eta_max;
eta_cur(ind) = eta_max;

% Brittle viscosity
Txxm_new        =   Txxm.*(1-(mum*Ddt)./(mum*Ddt+eta_cur))+2*eta_cur.*Exxm.*(mum*Ddt)./(mum*Ddt+eta_cur);
Txym_new        =   Txym.*(1-(mum*Ddt)./(mum*Ddt+eta_cur))+2*eta_cur.*Exym.*(mum*Ddt)./(mum*Ddt+eta_cur);

T2ndm_new       =   (Txxm_new.^2 + Txym_new.^2).^0.5;
T2ndm           =   (Txxm.^2 + Txym.^2).^0.5;

Yield       =   sind(phim).*(1 - lambdam) .* Pm + Cm .* cosd(phim);
ind = Yield < 0;
Yield(ind) = 0;

ind_plast = find(Yield < T2ndm_new & Im>0);
eta_cur(ind_plast)    =   0.5*Yield(ind_plast)./E2ndm(ind_plast);
Txxm(ind_plast)       =   Txxm(ind_plast).*Yield(ind_plast)./(T2ndm(ind_plast));
Txym(ind_plast)       =   Txym(ind_plast).*Yield(ind_plast)./(T2ndm(ind_plast));

ind = eta_cur > eta_max;
eta_cur(ind) = eta_max;
ind = eta_cur < eta_min;
eta_cur(ind) = eta_min;

eta_effm = eta_cur;
%==========================================================
