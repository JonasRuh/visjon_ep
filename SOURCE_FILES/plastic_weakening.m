% Calculating plastic strain
strainm(ind_plast) = strainm(ind_plast) + E2ndm(ind_plast)*dt;

% Strain weakeing
if strcmp('strain',weakening)==1
    ind = strainm > w_1;
    phim(ind)    = phi(Im(ind))-(strainm(ind)-w_1)/(w_2-w_1).*(phi(Im(ind))-phi_w(Im(ind)));
    Cm(ind)      = C(Im(ind))-(strainm(ind)-w_1)/(w_2-w_1).*(C(Im(ind))-C_w(Im(ind)));
    
    ind = strainm > w_2;
    phim(ind)    = phi_w(Im(ind));
    Cm(ind)      = C_w(Im(ind));
end


% Calculating plastic work
workm(ind_plast) = workm(ind_plast) + 2.*Txxm(ind_plast).*Exxm(ind_plast)*dt + 2.*Txym(ind_plast).*Exym(ind_plast)*dt;

% Work weakeing
if strcmp('work',weakening)==1
    ind = workm > w_1;
    phim(ind)    = phi(Im(ind))-(workm(ind)-w_1)/(w_2-w_1).*(phi(Im(ind))-phi_w(Im(ind)));
    Cm(ind)      = C(Im(ind))-(workm(ind)-w_1)/(w_2-w_1).*(C(Im(ind))-C_w(Im(ind)));
    
    ind = workm > w_2;
    phim(ind)    = phi_w(Im(ind));
    Cm(ind)      = C_w(Im(ind));
end



