% Surface process (diffusion)
if mod(timestep,Surface_dt)==0
    num1S=1;
    if SurfaceBC_left < 0
        iiS(num1S)   =   1;
        jjS(num1S)   =   1;
        ssS(num1S)   =   1;
        iiS(num1S+1) =   1;
        jjS(num1S+1) =   2;
        ssS(num1S+1) =  -1;
        RS(num1S)    =   0;
        num1S=3;
    else
        iiS(num1S)   =   1;
        jjS(num1S)   =   1;
        ssS(num1S)   =   1;
        RS(num1S)    =   SurfaceBC_left;
        num1S=2;
    end
    
    if SurfaceBC_right < 0
        iiS(num1S)   =   (nx-1)*surface_nodes+1;
        jjS(num1S)   =   (nx-1)*surface_nodes+1;
        ssS(num1S)   =   1;
        iiS(num1S+1) =   (nx-1)*surface_nodes+1;
        jjS(num1S+1) =   (nx-1)*surface_nodes+1-1;
        ssS(num1S+1) =  -1;
        RS((nx-1)*surface_nodes+1)       =   0;
        num1S=num1S+2;
    else
        iiS(num1S)   =   (nx-1)*surface_nodes+1;
        jjS(num1S)   =   (nx-1)*surface_nodes+1;
        ssS(num1S)   =   1;
        RS((nx-1)*surface_nodes+1)       =   SurfaceBC_right;
        num1S=num1S+1;
    end
    
    num2S = ((nx-1)*surface_nodes+1-3)+num1S;
    
    iiS(num1S:num2S) = 2:(nx-1)*surface_nodes+1-1;
    jjS(num1S:num2S) = 2:(nx-1)*surface_nodes+1-1;
    ssS(num1S:num2S) = 2*SurfaceCoeff/(dx/surface_nodes)^2 + 1/(time-time_old);
    
    num1S = num2S+1; num2S = ((nx-1)*surface_nodes+1-3)+num1S;
    
    iiS(num1S:num2S) = 2:(nx-1)*surface_nodes+1-1;
    jjS(num1S:num2S) = 3:(nx-1)*surface_nodes+1;
    ssS(num1S:num2S) = -SurfaceCoeff/(dx/surface_nodes)^2;
    
    num1S = num2S+1; num2S = ((nx-1)*surface_nodes+1-3)+num1S;
    
    iiS(num1S:num2S) = 2:(nx-1)*surface_nodes+1-1;
    jjS(num1S:num2S) = 1:(nx-1)*surface_nodes+1-2;
    ssS(num1S:num2S) = -SurfaceCoeff/(dx/surface_nodes)^2;
    RS(2:(nx-1)*surface_nodes+1-1) =   surface_y(2:(nx-1)*surface_nodes+1-1)/(time-time_old);
    
    LS = sparse(iiS,jjS,ssS,(nx-1)*surface_nodes+1,(nx-1)*surface_nodes+1);

    surface_y_new(1,:) = LS\RS;
    
    ind = find((surface_y_new >= WaterLevel & surface_y <= WaterLevel) | (surface_y_new <= WaterLevel & surface_y >= WaterLevel));
    surface_y(ind) = WaterLevel;

    % Sedimentation
    if SedimentationStyle == 1 || SedimentationStyle == 2 || SedimentationStyle == 4
        ind = find(surface_y_new > WaterLevel & surface_y_new < surface_y);
        surface_y(ind) = surface_y_new(ind);
    end
    
    % Erosion
    if SedimentationStyle == 1 || SedimentationStyle == 3 || SedimentationStyle == 4
        ind = find(surface_y_new < WaterLevel & surface_y_new > surface_y);
        surface_y(ind) = surface_y_new(ind);
    end
    
    nj  =   fix(xm/(dx/surface_nodes))+1;
    mid =   (xm-surface_x(nj))/(dx/surface_nodes);
    
    mym =   surface_y(nj)-mid.*(surface_y(nj)-surface_y(nj+1));
    
    %============================================================
    %============================================================
    
    if variable_lambda==2
        ind=ym>mym;
        lambdam(ind)=0.4+(ym(ind)-mym(ind))./(Ly-mym(ind))*(lambda_bottom-0.4);
        ind= lambdam>lambda_bottom;
        lambdam(ind)=lambda_bottom;
        ind = lambdam>0.8 & (xm<15e3 | xm>85e3);
        lambdam(ind)=0.8;
    elseif variable_lambda==3
        ind=ym>mym;
        lambdam(ind)=0.4+(ym(ind)-mym(ind))./(Ly-11e3)*(lambda_bottom-0.4);
        ind= lambdam>lambda_bottom;
        lambdam(ind)=lambda_bottom;
    end
    
    S_ind    =   find(ym > mym);
    I_ind    =   find(Im(S_ind) <= 2 & Im(S_ind) > 0);
    M_ind    =   S_ind(I_ind);
    
    if round(time/SecYear/SedChange)<time/SecYear/SedChange && length(SediMarker)==2
        ii=1;
    else
        ii=2;
    end
    
    Im(M_ind) = SediMarker(ii);
    % put values to minds
    etam(M_ind)    =   eta(Im(M_ind));
    mum(M_ind)     =   mu(Im(M_ind));
    rhom0(M_ind)   =   rho(Im(M_ind));
    Cm(M_ind)      =   C(Im(M_ind));
    phim(M_ind)    =   phi(Im(M_ind));
    lambdam(M_ind) =   lambda(Im(M_ind));
    if Temperature == 1
        if ra_heat == 1
            hrm(M_ind)     =   hr(Im(M_ind));
        end
        if adiab_heat == 1 && TPdep_dens == 0
            tam(M_ind)        =   ta(Im(M_ind));
        end
        if TPdep_dens == 1
            tam(M_ind)        =   ta(Im(M_ind));
            tbm(M_ind)        =   tb(Im(M_ind));
        end
        Tm(M_ind)      =   T(Im(M_ind));
        km(M_ind)      =   kk(Im(M_ind));
        cpm(M_ind)     =   cp(Im(M_ind));
        rhocpm0(M_ind) =   rho(Im(M_ind)).*cp(Im(M_ind));
    end
    if Powerlaw == 1
        nm(M_ind)      =   n(Im(M_ind));
        Qm(M_ind)      =   Q(Im(M_ind));
    end
    
    % Erosion
    S_ind    =   find(ym < mym);
    I_ind    =   find(Im(S_ind) > 2);
    M_ind    =   S_ind(I_ind);
    
    Im(M_ind) = ErosMarker;
    % put values to markers
    etam(M_ind)    =   eta(Im(M_ind));
    mum(M_ind)     =   mu(Im(M_ind));
    rhom0(M_ind)   =   rho(Im(M_ind));
    Cm(M_ind)      =   C(Im(M_ind));
    phim(M_ind)    =   phi(Im(M_ind));
    lambdam(M_ind) =   lambda(Im(M_ind));
    strainm(M_ind) =   0;
    workm(M_ind)   =   0;
    if Temperature == 1
        if ra_heat == 1
            hrm(M_ind)     =   hr(Im(M_ind));
        end
        if adiab_heat == 1 && TPdep_dens == 0
            tam(M_ind)        =   ta(Im(M_ind));
        end
        if TPdep_dens == 1
            tam(M_ind)        =   ta(Im(M_ind));
            tbm(M_ind)        =   tb(Im(M_ind));
        end
        Tm(M_ind)      =   T(Im(M_ind));
        km(M_ind)      =   kk(Im(M_ind));
        cpm(M_ind)     =   cp(Im(M_ind));
        rhocpm0(M_ind) =   rho(Im(M_ind)).*cp(Im(M_ind));
    end
    if Powerlaw == 1
        nm(M_ind)      =   n(Im(M_ind));
        Qm(M_ind)      =   Q(Im(M_ind));
    end

    clear nj surface_y_new;

    time_old = time;
end