
%----------------- ADD MARKERS TO EMPTY CELLS ---------------------
% Find empty cells RIGHT
right_vel = Vx((ny+1)*(nx-1)+2:1:(ny+1)*nx-1)';
right_income = right_income-right_vel.*dt;
ind=right_income<0;
right_income(ind)=0;

for i=1:ny-1
    if right_income(i)>dx/mx
        right_income(i)=right_income(i)-dx/mx;
        
        j = 1:my;
        ind_mark =  find(Im==-1);
        
        if length(j)>length(ind_mark)
            marker(1:length(ind_mark))              =    ind_mark;
            marker(length(ind_mark)+1:length(j))    =    marknum+[1:length(j)-length(ind_mark)];
            marknum    =   marker(end);
        else
            marker(1:length(j)) =   ind_mark(1:length(j));
        end
        
        xm(marker)  =   Lx - dx/2/mx + dx*(rand(size(j))-0.5)/mx - left_income(i);
        ym(marker)  =   y(i) + dy/my*(j-1) + dy/my/2 + dy*(rand(size(j))-0.5)/my;
        
        
        for iS = 1:size(SETUP,1)
            % Rock1
            ind=find(SETUP(iS,3) == Lx & ym(marker) > SETUP(iS,4) & ym(marker) < SETUP(iS,5));
            h(ind) = SETUP(iS,1);
        end
        
        for ii=1:length(pattern_marker)
            ind1=find((pattern_marker(ii)==h)==1);
            % Horizontal stripes
            if (strcmp('horizontal',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                ind2=find(round(ym(marker(ind1))./pattern_ydim)<ym(marker(ind1))./pattern_ydim);
                h(ind1(ind2))=h(ind1(ind2))+1;
            end
            % Vertical stripes
            if (strcmp('vertical',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                if (round(time/SecYear/(pattern_xdim/BC_right))>time/SecYear/(pattern_xdim/BC_right));
                    h(ind1)=h(ind1)+1;
                end
            end
        end
        
        etam(marker)       =   eta(h);
        mum(marker)        =   mu(h);
        rhom0(marker)      =   rho(h);
        Cm(marker)         =   C(h);
        phim(marker)       =   phi(h);
        Im(marker)         =   Phase(h);
        eta_effm(marker)   =   eta(h);
        lambdam(marker)    =   lambda(h);
        if Temperature == 1
            if ra_heat == 1
                hrm(marker)        =   hr(h);
            end
            if adiab_heat == 1 && TPdep_dens == 0
                tam(marker)        =   ta(h);
            end
            if TPdep_dens == 1
                tam(marker)        =   ta(h);
                tbm(marker)        =   tb(h);
            end
            %             Tm(marker)         =   T(h);
            km(marker)         =   kk(h);
            cpm(marker)        =   cp(h);
            rhocpm0(marker)    =   rho(h).*cp(h);
        end
        if Powerlaw == 1
            nm(marker)         =   n(h);
            Qm(marker)         =   Q(h);
        end
        
        %         Txxm(marker)       =   0;
        %         Txym(marker)       =   0;
        %         Exxm(marker)       =   0;
        %         Exym(marker)       =   0;
        %         E2ndm(marker)      =   0;
        %         Pm(marker)         =   0;
        %         wm(marker)         =   0;
        strainm(marker)     =   0;
        workm(marker)       =   0;
        
        if variable_lambda>=1
            find(ym(marker)>11e3);
            lambdam(marker)=lambdam(marker)+(ym(marker)-11e3)./4e3.*(lambda_bottom-0.4);
        end
        
        % Interpolate marker values from nodes
        jj = nx;
        ii = fix((ym(marker)+dy/2)/dy)+1;
        
        dxm =   xm(marker) - x_P(jj);
        dym =   ym(marker) - y_P(ii);
        
        ind_eq = (jj-1).*(ny+1)+ii;
        
        % for nodal points
        Exxm(marker)      =    Exx(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Exx(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Exx(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Exx(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        Txxm(marker)      =    Txx(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Txx(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Txx(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Txx(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        Pm(marker)        =    P(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +P(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +P(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
            +P(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        
        if Temperature == 1
            Tm(marker)        =    Temp(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
                +Temp(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
                +Temp(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
                +Temp(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        end
        
        jj = nx-1;
        ii = fix((ym(marker))/dy)+1;
        
        dxm =   xm(marker) - x_Vx(jj);
        dym =   ym(marker) - y_Vy(ii);
        
        ind_eq = (jj-1).*(ny+1)+ii;
        
        Exym(marker)      =    Exy(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Exy(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Exy(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Exy(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        E2ndm(marker)     =    E2nd_s(ind_eq)'      .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +E2nd_s(ind_eq+(ny+1))'                 .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +E2nd_s(ind_eq+1)'                      .* (1 - (dxm)/dx) .* (dym/dy) ...
            +E2nd_s(ind_eq+(ny+1)+1)'               .* (dxm/dx)       .* (dym/dy);
        
        Txym(marker)      =    Txy(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Txy(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Txy(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Txy(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        wm(marker)        =    w(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +w(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +w(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
            +w(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        
        clear marker
    end
end

%----------------- ADD MARKERS TO EMPTY CELLS ---------------------
% Find empty cells LEFT

left_vel = Vx(2:1:ny)';
left_income = left_income+left_vel.*dt;
ind=left_income<0;
left_income(ind)=0;

for i=1:ny-1
    if left_income(i)>dx/mx
        left_income(i)=left_income(i)-dx/mx;
        
        j = 1:my;
        ind_mark =  find(Im==-1);
        
        if length(j)>length(ind_mark)
            marker(1:length(ind_mark))              =    ind_mark;
            marker(length(ind_mark)+1:length(j))    =    marknum+[1:length(j)-length(ind_mark)];
            marknum    =   marker(end);
        else
            marker(1:length(j)) =   ind_mark(1:length(j));
        end
        
        xm(marker)  =   x(1) + dx/2/mx + dx*(rand(size(j))-0.5)/mx + left_income(i);
        ym(marker)  =   y(i) + dy/my*(j-1) + dy/my/2 + dy*(rand(size(j))-0.5)/my;
        
        for iS = 1:size(SETUP,1)
            % Rock1
            ind=find(SETUP(iS,2) == 0 & ym(marker) > SETUP(iS,4) & ym(marker) < SETUP(iS,5));
            h(ind) = SETUP(iS,1);
        end
        
        for ii=1:length(pattern_marker)
            ind1=find((pattern_marker(ii)==h)==1);
            % Horizontal stripes
            if (strcmp('horizontal',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                ind2=find(round(ym(marker(ind1))./pattern_ydim)<ym(marker(ind1))./pattern_ydim);
                h(ind1(ind2))=h(ind1(ind2))+1;
            end
            % Vertical stripes
            if (strcmp('vertical',pattern_type)==1 || strcmp('kaki',pattern_type)==1)
                if (round(time/SecYear/(pattern_xdim/BC_right))>time/SecYear/(pattern_xdim/BC_right));
                    h(ind1)=h(ind1)+1;
                end
            end
        end
        
        etam(marker)       =   eta(h);
        mum(marker)        =   mu(h);
        rhom0(marker)      =   rho(h);
        Cm(marker)         =   C(h);
        phim(marker)       =   phi(h);
        Im(marker)         =   Phase(h);
        eta_effm(marker)   =   eta(h);
        lambdam(marker)    =   lambda(h);
        if Temperature == 1
            if ra_heat == 1
                hrm(marker)        =   hr(h);
            end
            if adiab_heat == 1 && TPdep_dens == 0
                tam(marker)        =   ta(h);
            end
            if TPdep_dens == 1
                tam(marker)        =   ta(h);
                tbm(marker)        =   tb(h);
            end
            %             Tm(marker)         =   T(h);
            km(marker)         =   kk(h);
            cpm(marker)        =   cp(h);
            rhocpm0(marker)    =   rho(h).*cp(h);
        end
        if Powerlaw == 1
            nm(marker)         =   n(h);
            Qm(marker)         =   Q(h);
        end
        
        %         Txxm(marker)       =   0;
        %         Txym(marker)       =   0;
        %         Exxm(marker)       =   0;
        %         Exym(marker)       =   0;
        %         E2ndm(marker)      =   0;
        %         Pm(marker)         =   0;
        %         wm(marker)         =   0;
        strainm(marker)     =   0;
        workm(marker)       =   0;

        % Interpolate marker values from nodes
        jj = 1;
        ii = fix((ym(marker)+dy/2)/dy)+1;
        
        dxm =   xm(marker) - x_P(jj);
        dym =   ym(marker) - y_P(ii);
        
        ind_eq = (jj-1).*(ny+1)+ii;
        
        % for nodal points
        Exxm(marker)      =    Exx(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Exx(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Exx(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Exx(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        Txxm(marker)      =    Txx(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Txx(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Txx(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Txx(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        Pm(marker)        =    P(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +P(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +P(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
            +P(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        
        if Temperature == 1
            Tm(marker)        =    Temp(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
                +Temp(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
                +Temp(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
                +Temp(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        end
        
        jj = 1;
        ii = fix((ym(marker))/dy)+1;
        
        dxm =   xm(marker) - x_Vx(jj);
        dym =   ym(marker) - y_Vy(ii);
        
        ind_eq = (jj-1).*(ny+1)+ii;
        
        Exym(marker)      =    Exy(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Exy(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Exy(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Exy(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        E2ndm(marker)     =    E2nd_s(ind_eq)'      .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +E2nd_s(ind_eq+(ny+1))'                 .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +E2nd_s(ind_eq+1)'                      .* (1 - (dxm)/dx) .* (dym/dy) ...
            +E2nd_s(ind_eq+(ny+1)+1)'               .* (dxm/dx)       .* (dym/dy);
        
        Txym(marker)      =    Txy(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Txy(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Txy(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Txy(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        wm(marker)        =    w(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +w(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +w(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
            +w(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        
        clear marker
    end
end

%----------------- ADD MARKERS TO EMPTY CELLS ---------------------
% Find empty cells TOP
top_vel = Vy((ny+1)+1:(ny+1):(ny+1)*(nx-1)+1)';
top_income = top_income+top_vel.*dt;
ind=top_income<0;
top_income(ind)=0;

for i=1:nx-1
    if top_income(i)>dy/my
        top_income(i)=top_income(i)-dy/my;
        
        j = 1:mx;
        ind_mark =  find(Im==-1);
        
        if length(j)>length(ind_mark)
            marker(1:length(ind_mark))              =    ind_mark;
            marker(length(ind_mark)+1:length(j))    =    marknum+[1:length(j)-length(ind_mark)];
            marknum    =   marker(end);
        else
            marker(1:length(j)) =   ind_mark(1:length(j));
        end
        
        ym(marker)  =   y(1) + dy/2/my + dy*(rand(size(j))-0.5)/my + top_income(i);
        xm(marker)  =   x(i) + dx/mx*(j-1) + dx/mx/2 + dx*(rand(size(j))-0.5)/mx;

        h = mark_top;
        
        etam(marker)       =   eta(h);
        mum(marker)        =   mu(h);
        rhom0(marker)      =   rho(h);
        Cm(marker)         =   C(h);
        phim(marker)       =   phi(h);
        Im(marker)         =   Phase(h);
        eta_effm(marker)   =   eta(h);
        lambdam(marker)    =   lambda(h);
        if Temperature == 1
            if ra_heat == 1
                hrm(marker)        =   hr(h);
            end
            if adiab_heat == 1 && TPdep_dens == 0
                tam(marker)        =   ta(h);
            end
            if TPdep_dens == 1
                tam(marker)        =   ta(h);
                tbm(marker)        =   tb(h);
            end
            Tm(marker)         =   top_T;
            km(marker)         =   kk(h);
            cpm(marker)        =   cp(h);
            rhocpm0(marker)    =   rho(h)*cp(h);
        end
        if Powerlaw == 1
            nm(marker)         =   n(h);
            Qm(marker)         =   Q(h);
        end
        
        Txxm(marker)       =   0;
        Txym(marker)       =   0;
        Exxm(marker)       =   0;
        Exym(marker)       =   0;
        E2ndm(marker)      =   0;
        Pm(marker)         =   0;
        wm(marker)         =   0;
        strainm(marker)     =   0;
        workm(marker)       =   0;
        
        clear marker
    end
end


%----------------- ADD MARKERS TO EMPTY CELLS ---------------------
% Find empty cells BOTTOM if open lower boundary

bottom_vel = Vy((ny+1)+(ny):(ny+1):(ny+1)*(nx)-1)';
bottom_income = bottom_income-bottom_vel.*dt;
ind=bottom_income<0;
bottom_income(ind)=0;

for i=1:nx-1
    if bottom_income(i)>dy/my
        bottom_income(i)=bottom_income(i)-dy/my;
        
        j=1:mx;
        ind_mark =  find(Im==-1);
        
        if length(j)>length(ind_mark)
            marker(1:length(ind_mark))              =    ind_mark;
            marker(length(ind_mark)+1:length(j))    =    marknum+[1:length(j)-length(ind_mark)];
            marknum    =   marker(end);
        else
            marker(1:length(j)) =   ind_mark(1:length(j));
        end
        
        ym(marker)  =   Ly - dy/my/2 + dy*(rand(size(j))-0.5)/my - bottom_income(i);
        xm(marker)  =   x(i) + dx/mx*(j-1) + dx/mx/2 + dx*(rand(size(j))-0.5)/mx;
        
        h = mark_bottom;
        
        etam(marker)       =   eta(h);
        mum(marker)        =   mu(h);
        rhom0(marker)      =   rho(h);
        Cm(marker)         =   C(h);
        phim(marker)       =   phi(h);
        Im(marker)         =   Phase(h);
        eta_effm(marker)   =   eta(h);
        lambdam(marker)    =   lambda(h);
        if Temperature == 1
            if ra_heat == 1
                hrm(marker)        =   hr(h);
            end
            if adiab_heat == 1 && TPdep_dens == 0
                tam(marker)        =   ta(h);
            end
            if TPdep_dens == 1
                tam(marker)        =   ta(h);
                tbm(marker)        =   tb(h);
            end
            if Ext_T == 1
                Tm(marker)         =   bottom_T-M_grad.*(Ly-ym(marker))./1e3;
            end
            km(marker)         =   kk(h);
            cpm(marker)        =   cp(h);
            rhocpm0(marker)    =   rho(h)*cp(h);
        end
        if Powerlaw == 1
            nm(marker)         =   n(h);
            Qm(marker)         =   Q(h);
        end
        strainm(marker)     =   0;
        workm(marker)       =   0;
        
        dxm =   xm(marker) - x_P(i);
        dym =   ym(marker) - (Ly-dy/2);
        
        ind_eq = (i-1).*(ny+1)+ny;
        
        % for nodal points
        Exxm(marker)      =    Exx(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Exx(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Exx(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Exx(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        Txxm(marker)      =    Txx(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Txx(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Txx(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Txx(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        Pm(marker)        =    P(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +P(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +P(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
            +P(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        if Temperature == 1 && Ext_T < 1
            Tm(marker)        =    Temp(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
                +Temp(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
                +Temp(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
                +Temp(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        end
        
        
        dxm =   xm(marker) - x_Vx(i);
        dym =   ym(marker) - (Ly-dy);
        
        ind_eq = (i-1).*(ny+1)+ny-1;
        
        Exym(marker)      =    Exy(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Exy(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Exy(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Exy(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        E2ndm(marker)     =    E2nd_s(ind_eq)'      .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +E2nd_s(ind_eq+(ny+1))'                 .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +E2nd_s(ind_eq+1)'                      .* (1 - (dxm)/dx) .* (dym/dy) ...
            +E2nd_s(ind_eq+(ny+1)+1)'               .* (dxm/dx)       .* (dym/dy);
        
        Txym(marker)      =    Txy(ind_eq)'         .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +Txy(ind_eq+(ny+1))'                    .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +Txy(ind_eq+1)'                         .* (1 - (dxm)/dx) .* (dym/dy) ...
            +Txy(ind_eq+(ny+1)+1)'                  .* (dxm/dx)       .* (dym/dy);
        
        wm(marker)        =    w(ind_eq)'           .* (1 - (dxm/dx)) .* (1 - (dym/dy)) ...
            +w(ind_eq+(ny+1))'                      .* (dxm/dx)       .* (1 - (dym/dy)) ...
            +w(ind_eq+1)'                           .* (1 - (dxm)/dx) .* (dym/dy) ...
            +w(ind_eq+(ny+1)+1)'                    .* (dxm/dx)       .* (dym/dy);
        
        clear marker
    end
end

