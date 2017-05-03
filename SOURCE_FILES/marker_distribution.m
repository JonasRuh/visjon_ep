% Defining values of markers
marker = 0;
for im = 1:nym
    for jm = 1:nxm
        marker  =   marker+1;
        xm(marker)  =   Lx/nxm/2 + Lx/nxm*(jm-1) + (rand-0.5)*Lx/nxm;
        ym(marker)  =   Ly/nym/2 + Ly/nym*(im-1) + (rand-0.5)*Ly/nym;
        
        for iS = 1:size(SETUP,1)
            % Rock1
            if xm(marker) > SETUP(iS,2) && xm(marker) < SETUP(iS,3) && ym(marker) > SETUP(iS,4) && ym(marker) < SETUP(iS,5)
                rocktype = SETUP(iS,1);
            end
        end
        
        if ~isempty(find(pattern_marker==Phase(rocktype))==1)
            % Horizontal stripes
            if (round(ym(marker)/pattern_ydim)<ym(marker)/pattern_ydim) && (strcmp('horizontal',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                rocktype=rocktype+1;
            end
            % Vertical stripes            
            if (round(xm(marker)/pattern_xdim)<xm(marker)/pattern_xdim) && (strcmp('vertical',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                rocktype=rocktype+1;
            end
        end
        
        % put values to markers
        etam(marker)    =   eta(rocktype);
        mum(marker)     =   mu(rocktype);
        rhom0(marker)   =   rho(rocktype);
        Cm(marker)      =   C(rocktype);
        phim(marker)    =   phi(rocktype);
        lambdam(marker) =   lambda(rocktype);
        Im(marker)      =   Phase(rocktype);
        if Temperature == 1
            if ra_heat == 1
                hrm(marker)     =   hr(rocktype);
            end
            if adiab_heat == 1 && TPdep_dens == 0
                tam(marker)     =   ta(rocktype);
            end
            if TPdep_dens == 1
                tam(marker)     =   ta(rocktype);
                tbm(marker)     =   tb(rocktype);
            end
            Tm(marker)      =   T(rocktype);
            km(marker)      =   kk(rocktype);
            cpm(marker)     =   cp(rocktype);
            rhocpm0(marker) =   rho(rocktype)*cp(rocktype);
        end
        if Powerlaw == 1
            nm(marker)      =   n(rocktype);
            Qm(marker)      =   Q(rocktype);
        end
        if ym(marker)>11e3 && variable_lambda>=1
            lambdam(marker)=lambdam(marker)+(ym(marker)-11e3)/4e3*(lambda_bottom-0.4);
        end
        if lambdam(marker)>0.8 && (xm(marker)<15e3 || xm(marker)>85e3)
            lambdam(marker)=0.8;
        end
    end
end

% % Set initial temperature distribution
% if Temperature == 1
%     for marker=1:marknum
%         Tm(marker)=273;
%         % DEEPER MANTLE
%         if ym(marker)>=A_thick
%             Tm(marker)=273+L_A_temp+M_grad*(ym(marker)-A_thick-L_thick)/1e3;
%         end
%         % CRUST
%         if ym(marker)>=A_thick && ym(marker)<=A_thick+C_thick
%             Tm(marker)=273+Moho_temp/C_thick*(ym(marker)-A_thick);
%         end
%         
%         % PERTUBATION
%         value1 = (sin((xm(marker)-Pert_beg)/((Pert_end-Pert_beg)/2)*pi-pi/2)+1)/2;
%         
%         % LITHOSPHERIC MANTLE
%         if ym(marker)>=A_thick+C_thick && ym(marker)<=A_thick+L_thick && (xm(marker)<Pert_beg || xm(marker)>Pert_end)
%             Tm(marker)=273+Moho_temp+(L_A_temp-Moho_temp)/(L_thick-C_thick)*(ym(marker)-A_thick-C_thick);
%             % Lithosphere
%         elseif ym(marker)>=A_thick+C_thick && ym(marker)<=A_thick+L_thick-value1*Pert_add && (xm(marker)>=Pert_beg && xm(marker)<=Pert_end)
%             % y = (15*90/(15+2) - 90) / (0.5/(15+2)-1)
%             Tm(marker)=273+Moho_temp+(L_A_temp-(M_grad*Pert_add/1e3)*value1-Moho_temp)/(L_thick-Pert_add*value1-C_thick)*(ym(marker)-A_thick-C_thick);
%         end 
%     end 
% end

% Set initial temperature distribution
if Temperature == 1
    for marker=1:marknum
        Tm(marker)=273;
        % Background (mantle)
        if Im(marker)>2
            Tm(marker)=(L_thick/1e3*L_grad+273-M_grad*(L_thick+A_thick)/1e3)+ym(marker)/1e3*M_grad;
        end
        
%         % Lithosphere
%         if xm(marker)>=Pert_beg && xm(marker)<=Pert_end && Im(marker)>2 && Im(marker)<12
%             
%             value1 = (sin((xm(marker)-Pert_beg)/((Pert_end-Pert_beg)/2)*pi-pi/2)+1)/2;
%             
%             %      y = (15*90/(15+2) - 90) / (0.5/(15+2)-1)
%             if ym(marker)<=L_thick+A_thick - ((L_grad*L_thick/(L_grad+Pert_add))-L_thick) / (M_grad/(L_grad+Pert_add)-1)*value1;
%                 Tm(marker)=273+(ym(marker)-A_thick)*(L_grad+Pert_add*(value1))/1e3;  
%             end
%         elseif xm(marker)<Pert_beg || xm(marker)>Pert_end
%             if Im(marker)>2 && ym(marker)<=L_thick+A_thick
%                 Tm(marker)=273+(ym(marker)-A_thick)*L_grad/1e3;
%             end           
%         end
        
        % Lithosphere different
        if xm(marker)>=Pert_beg && xm(marker)<=Pert_end && Im(marker)>2
            
            value1 = (sin((xm(marker)-Pert_beg)/((Pert_end-Pert_beg)/2)*pi-pi/2)+1)/2;
            
            Tm(marker)=273+(ym(marker)-A_thick)*(L_grad+(value1*Pert_add))/1e3;
            
        elseif xm(marker)<Pert_beg || xm(marker)>Pert_end
            if Im(marker)>2 && ym(marker)<=L_thick+A_thick
                Tm(marker)=273+(ym(marker)-A_thick)*L_grad/1e3;
            end           
        end

    end
end    
value2 = (sin((x_Vy(2:end-1)-Pert_beg)/((Pert_end-Pert_beg)/2)*pi-pi/2)+1)/2;
ind=find(x_Vy(2:end-1)<Pert_beg);
value2(ind)=0;
ind=find(x_Vy(2:end-1)>Pert_end);
value2(ind)=0;
bottom_T    =   L_thick/1e3*(L_grad+(Pert_add.*value2))+273+(Ly-L_thick-A_thick)/1e3*M_grad;

% 
% value2 = (sin((x_Vy(2:end-1)-Pert_beg)/((Pert_end-Pert_beg)/2)*pi-pi/2)+1)/2;
% ind=find(x_Vy(2:end-1)<Pert_beg);
% value2(ind)=0;
% ind=find(x_Vy(2:end-1)>Pert_end);
% value2(ind)=0;
% bottom_T    =   bottom_T + value2.*Pert_add;

