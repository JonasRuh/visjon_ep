%% 3D Conservative Stokes / Marker in cell / plastic iterations
% VisoElastoPlasic approach
% by Jonas, March/April 2012
clear all; 

addpath('SOURCE_FILES')                        % FUNCTIONS
mkdir OUTPUT

create_output       = 50;
create_breakpoint   = 500;

brutus              = 1;
Temperature         = 1;
Powerlaw            = 1;
Surface_diffusion   = 1;

if brutus == 0
    clf;
    load Colormaps;
end

% Length of model
Lx   =  500e3;
Ly   =   35e3;

% number of nodes
nx      =    2001;
ny      =     141;
p_node  =    1999;   % between 3 and nx-2

% markers per node
mx  =   4;
my  =   4;

        % Phase    Visc       EModul    Dens     Phi     PhiW      Coh         CohW      lamb      Hr       Ta       Tb     cp      Temp       kk       n       Q       
ROCKS= [    1      1e17        1e11        1       0       1       1e20        1e20        0        0        0        0    3e6       273      200       1       0   ;   % AIR
            2      1e17        1e11      1e3       0       1       1e20        1e20        0        0        0        0    3e3       273      200       1       0   ;   % WATER

            3    2.0e17        1e11     2500      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     2.3   154e3   ;   % Quartzite1
            4    2.0e17        1e11     2500      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     2.3   154e3   ;   % Quartzite1
            5    2.0e17        1e11     2500      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     2.3   154e3   ;   % Quartzite1
            6    2.0e17        1e11     2500      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     2.3   154e3   ;   % Quartzite1
            
            7    2.0e17        1e11     2800      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     2.3   154e3   ;   % Quartzite2
            8    2.0e17        1e11     2800      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     2.3   154e3   ;   % Quartzite2
            9    4.8e22        1e11     3000      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     3.2   238e3   ;   % Quartzite2_lower_OC
           10    4.8e22        1e11     3000      30      15        2e7         2e5      0.0     2e-6     2e-5   45e-13    1e3       273      2.5     3.2   238e3   ;   % Quartzite2_lower_OC
            
           11    1.2e16        1e11     2800      30      15        2e7         2e5      0.0     2e-7     2e-5   45e-13    1e3       273      2.5     2.4   212e3   ;   % Quartz-diorite1
           12    1.2e16        1e11     2800      30      15        2e7         2e5      0.0     2e-7     2e-5   45e-13    1e3       273      2.5     2.4   212e3   ;   % Quartz-diorite1

           13    4.0e16        1e11     3300      30      15        2e7         2e7      0.0     2e-8     2e-5   45e-13    1e3       273      2.5     3.5   532e3   ;   % Mantle1
           14    4.0e16        1e11     3300      30      15        2e7         2e7      0.0     2e-8     2e-5   45e-13    1e3       273      2.5     3.5   532e3   ;  % Mantle2
           15    4.0e16        1e11     3300      30      15        2e7         2e7      0.0     2e-8     2e-5   45e-13    1e3       273      2.5     3.5   532e3   ;   % Mantle1
           16    4.0e16        1e11     3300      30      15        2e7         2e7      0.0     2e-8     2e-5   45e-13    1e3       273      2.5     3.5   532e3   ];  % Mantle2

       
Phase   =   ROCKS(:,1)';
eta     =   ROCKS(:,2)';
mu      =   ROCKS(:,3)';
rho     =   ROCKS(:,4)';
phi     =   ROCKS(:,5)';
phi_w   =   ROCKS(:,6)';
C       =   ROCKS(:,7)';
C_w     =   ROCKS(:,8)';
lambda  =   ROCKS(:,9)';
if Temperature == 1
    hr  =   ROCKS(:,10)';
    ta  =   ROCKS(:,11)';
    tb  =   ROCKS(:,12)';
    cp  =   ROCKS(:,13)';
    T   =   ROCKS(:,14)';
    kk  =   ROCKS(:,15)';
end
if Powerlaw == 1
    n   =   ROCKS(:,16);
    Q   =   ROCKS(:,17);
end

variable_lambda=0;      % 1: initial lambda; 2: equalized whole column; 3: equalized top 4 km
lambda_bottom=0.9;

% Geometry      Phase       x-min       x-max       y-min       y-max
SETUP   =   [   1           0           Lx          0           Ly      ;
               11         0e3           Lx     05.0e3       35.0e3      ];

% Plastic weakeing thresholds
weakening    = {'strain','velocity','work'};
weakening    = weakening(3);
w_1 = 0;
w_2 = 200e6;

% Parameters
gravity_y       =   9.81;
gravity_x       =   00;
SecYear         =   3600*24*365.25;

% Surface process
SedimentationStyle  =   2;      % 1 = sed and ero || 2 = only sed || 3 = only ero || 4 = linear below sealevel
SurfaceCoeff        =   1e-5;
SurfaceBC_left      =   5e3;       % For free slip: -1
SurfaceBC_right     =   5e3;       % For free slip: -1
Surface_dt          =   1;         % makes surface process every xx timestep
time_old            =   0;
SediMarker          =  [3 4];       % Marker type for sedimentation
ErosMarker          =   1;          % Marker type for erosion
SedChange           =   1e6;        % Marker change for sedimentation in years
WaterLevel          =   0;          % If WaterLevel == 0 it is switched off !!
surface_nodes       =   1;          % Times x nodes along surface
surface_init        =   5e3;
surface_x           =   [0:(Lx/(nx-1))/surface_nodes:Lx];
surface_y           =   surface_init*ones(size(surface_x));
surface_smoother    =   1;

% Cutoff viscosities
eta_max = 1e24;
eta_min = 1e17;

if Temperature==1
    % Temperature boundary conditions
    A_thick     =   05e3; % Thickness of sticky-air
    C_thick     =   30e3; % Thickness of Crust
    L_thick     =   30e3; % Thickness of Lithosphere
    Moho_temp   =   540;  % Temperature at Moho in °C
    L_A_temp    =   540; % Temperature at Lithosphere/Asthenosphere boundary in °C
    M_grad      =   0.5;  % Temperature gradient in Mantle (background)
    Pert_beg    =   125e3;
    Pert_end    =   375e3;
    Pert_add    =   3; % km difference for Lithosphere/Asthenosphere depth
    top_T       =   273;
    Ext_T       =   1;
    if Ext_T == 1
        bottom_T    =   L_A_temp+273+M_grad*(Ly-A_thick-L_thick)/1e3;
    else
        bottom_T    =   2;
    end
    shear_heat  =   0;
    ra_heat     =   0;
    adiab_heat  =   0;
    TPdep_dens  =   0;
    
L_grad      =   13;   % Temperature gradient in Lithosphere
M_grad      =   13;  % Temperature gradient in Mantle (background)

end

% Velocity boundary conditions

%   L T T T T T T T T T T T RR      T T T T T T T T T T T T T   
%   L                       RR      L                       R
%   L      x-velocity       RR      L      y-velocity       R
%   L                       RR      L                       R
%   L                       RR      B B B B B B B B B B B B B
%   L B B B B B B B B B B B RR      B B B B B B B B B B B B B


BC_left_init     = -0.005;
BC_right_init    =  0.005;
BC_top_init      =  0.01*35/500;
BC_bottom_init   =(-0.000*387/1000)/200;
bound            = {'freeslip','noslip','velocity','external','mixed'};
top_BC           = bound(1);
bottom_BC        = bound(1);
left_BC          = bound(1);
right_BC         = bound(1);
left_mix         = -ones(1,ny-2);
right_mix        = -ones(1,ny-2);
left_mix(1:110)  = 1;
right_mix(1:110) = 1;
top_mix          = ones(1,nx-2);
bottom_mix       = [ones(1,399)         -ones(1,1201)  ones(1,399); ...
                   -0.005*ones(1,399)    zeros(1,1201) 0.005*ones(1,399)];
Ext_BC           = 0.995;

% Inversion
inversion = [30e6 37.5e6 37.5e6 310e6];     % Time in Ma for velocity inversion

% Incoming marker type
mark_top    = 1;
mark_bottom = 16;

%=============================================
% Initialize matrices and coordinates
initialize_beg;
%=============================================

% Initial marker pattern
pattern_marker  =   [11];     % Types of markers with pattern (check marker phases !!!!)
pattern_type    =   {'horizontal','vertical','kaki'};   % Stripes or kaki
pattern_type    =   pattern_type(1);
pattern_xdim    =   2.5e3;    % thickness of vertical stripes
pattern_ydim    =   2.5e3;    % thickness of horizontal stripes

%=============================================
% Initialize phase and temperature distribution of marker 
marker_distribution;
%=============================================

% Initial time and iteration setup
time        = 0;
t_beg       = 1;
dt_value    = 0.25;          % threshold value for maximal timestep, 0.1 = 10% movement of dx or dy
Ddt         = 1000*SecYear;  % Initial timestep to pre-stress
short_dt    = 1000*SecYear;  % Short initial timestep
dt_max      = 5000*SecYear;  % Timestep
n_short_dt  = 25;            % Number of short initial timesteps
miniter     = 3;             % minimal number of iterations
maxiter     = 50;             % maximal number of iterations
error       = 1;             % 1 = Error 1 (average nodal velocity change); 2 = Error 2 (largest nodal velocity change)
vel_res     = 1e-14;         % sum(velocity) change fot iter break
strain_rate_smoother = 0.87;

% Load breakpoint file if necessary
if exist(char('Breakpoint.mat'))
    load Breakpoint.mat
    t_beg    = timestep+1;
end

%=========================== START TIME LOOP ==========================
%======================================================================

for timestep = t_beg:1e6
    tic

    dt=dt_max;
    if timestep <= n_short_dt
        dt = short_dt;  
    end
    
    velocity_inversion;
    
    for niter = 1:maxiter

        fprintf('Timestep: %d\n', timestep);
        fprintf('Iteration: %d\n', niter);
        
        %=============================================
        % Reload old values and initialize matrices
        initialize_iter;
        %=============================================
        
        %=============================================
        % Calculate power-law and brittle viscosities
        viscosity_calculation;
        %=============================================
        
        %=============================================
        % Fill nodal values from marker information
        marker_to_nodes; % dispatched loops, better for desktop
        toc, fprintf('for updating nodal values!\n'); 
        %=============================================   
        
        if Temperature == 1            
            % Boundary conditions to interpolated T
            % Upper BC
            Temp((ny+1)+1:(ny+1):(nx-1)*(ny+1)+1)   =   top_T*2 - Temp((ny+1)+2:(ny+1):(nx-1)*(ny+1)+2);
            % Lower BC
            if Ext_T == 1
                Temp(2*(ny+1):(ny+1):(nx)*(ny+1))       =   bottom_T'*2 - Temp(2*(ny)+1:(ny+1):(nx)*(ny+1)-1);
            else
                Temp(2*(ny+1):(ny+1):(nx)*(ny+1))       =   bottom_T + Temp(2*(ny)+1:(ny+1):(nx)*(ny+1)-1);
            end
            % Left BC
            Temp(1:ny+1)                            =   Temp(ny+2:2*(ny+1));
            % Right BC
            Temp(nx*(ny+1)+1:(nx+1)*(ny+1))         =   Temp((nx-1)*(ny+1)+1:(nx)*(ny+1));
        end
        
        %=====================================================
        
        %=============================================
        % Fill and solve matrix for stokes
        stokes_direct_solver;
        toc, fprintf('for solving the matrix!\n')
        %=============================================
                        
        % define optimal timestep
        Ddt = dt;
        
        Vx_max  =   max(max(abs(Vx2d_s)));
        Vy_max  =   max(max(abs(Vy2d_s)));
        
        if dt_value*dx/Vx_max < Ddt
            Ddt  =   dt_value*dx/Vx_max;
        end
        if dt_value*dy/Vy_max < Ddt
            Ddt  =   dt_value*dy/Vy_max;
        end
                
        %=============================================
        % Calculating strain rates and stresses
        strain_rate_stress_calculation;
        toc, fprintf('for strain/stress calculation and interpolation to markers!\n');
             fprintf('========= Ddt = %d years =============================\n',Ddt/SecYear);
        %=============================================
        
        fprintf('========= Min Pressure = %d MPa ===========================\n',min(P)/1e6);
        fprintf('========= Max Pressure = %d MPa ===========================\n',max(P)/1e6);
   
        %=============================================
        % Calculating strain rates and stresses
        vel_res_calculation;
        if timestep>1
            fprintf('========= VELOCITY ERROR %d = %d =========\n\n',error,iterations(error,niter+maxiter*(timestep-1)));
        end
        %=============================================

        % Plotting if running on desktop
        if brutus==0 && timestep>1
            E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
            eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
            strain_2d   =   reshape(strain,ny+1,nx+1);
            lambda_s_2d   =   reshape(lambda_s,ny+1,nx+1);
            
            figure(1), clf
            colormap jet
            
            subplot(321)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(E2nd_s_2d(1:end-1,1:end-1))), caxis([-18 -12])
            shading interp
            colorbar
            axis image, axis ij
            title('E2nd [1/s]')
            
            subplot(322)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(eta_s_2d(1:end-1,1:end-1))), caxis(log10([eta_min eta_max]))
            shading interp
            colorbar
            axis image, axis ij
            title(['\eta_{s} [Pa.s] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(323)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(lambda_s_2d(1:end-1,1:end-1)))
            shading interp
            colorbar
            axis image, axis ij
            title(['Lambda [-] ',    num2str(Ddt/SecYear)])
            drawnow
            
            subplot(324)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(T2nd_s_2d(1:end-1,1:end-1)))
            shading interp
            colorbar
            axis image, axis ij
            title(['\tau_{s} [Pa] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(325)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),Vx2d_s*SecYear)
            shading interp
            colorbar
            axis image, axis ij
            title(['V_{x} [m/s] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(326)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),Vy2d_s*SecYear)
            shading interp
            colorbar
            axis image, axis ij
            title(['V_{y} [m/s] ',    num2str(Ddt/SecYear)])
            drawnow

            
            figure(444)
            subplot(211), semilogy((niter-1)/maxiter+(timestep-1),iterations(1,niter+maxiter*(timestep-1)),'ko'), hold on
            title('Error 1: with sticky-air')
            subplot(212), semilogy((niter-1)/maxiter+(timestep-1),iterations(2,niter+maxiter*(timestep-1)),'ko'), hold on
            title('Error 2: with sticky-air')
        end
        
        % Exiting the iteration loop
        if (niter>=miniter && iterations(error,niter+maxiter*(timestep-1))<vel_res) || timestep==1
            break;
        end
    end

    %==========================================================================
    dt = Ddt;
    
    plastic_weakening;
    %==========================================================================

    %=============================================
    % Calculating temperature
    subgrid_diffusion_stresses;
    stress_rotation;
    %=============================================
    
    %=============================================
    % Calculating temperature
    if Temperature == 1
        temperature_direct_solver;
        toc, fprintf('for temperature calculation!\n');
    end
    %=============================================


    %=============================================
    % Move markers
    move_marker;
    move_surface;
    toc, fprintf('for moving makers!\n');
    %=============================================

    time = time + dt;

    %======================================================================
    % Visualization begin
    %======================================================================


    if brutus == 0
        figure(3), clf
        
        for i=1:20
            Phase1  =   find(Im == i);
            marksize = 5;
            hold on
            plot(xm(Phase1),ym(Phase1),'.','Color',Colormaps(i,:),'MarkerSize',marksize)
        end
        
        plot(surface_x,surface_y,'r'), hold on
        
        axis image, axis ij
        title(['Time: ',num2str((time)/SecYear/1e6),' Ma'])
%         plot(x2d_b,y2d_b,'k',x2d_b',y2d_b','k')
%         quiver(x2d_b,y2d_b,Vx2db,Vy2db,'r','LineWidth',1.5)
        E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
        eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
        T2nd_s_2d   =   reshape(T2nd_s,ny+1,nx+1);
                
        hold off
                
        figure(8)
        plot(surface_x,surface_y,'b'), hold off
        axis ij
        
        figure(4), clf
        colormap jet
        
        subplot(311)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(E2nd_s_2d(1:end-1,1:end-1))), caxis([-16 -12])
        shading interp
        colorbar
        axis image, axis ij
        title('E2nd ')
        
        subplot(312)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(eta_s_2d(1:end-1,1:end-1))), caxis(log10([eta_min eta_max]))
        shading interp
        colorbar
        axis image, axis ij
        title(['\eta_{s}   ',    num2str(dt/SecYear)])
        
        subplot(313)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(T2nd_s_2d(1:end-1,1:end-1)))
        shading interp
        colorbar
        axis image, axis ij
        title(['T2nd'])
        drawnow
    end

    %==========================================================================
    % Visualization end
    %==========================================================================

    %=============================================
    % Delete markers out of grid
    outgrid_marker;
    %=============================================

    %=============================================
    % New markers from sides
    incoming_marker;
    toc, fprintf('for outgoing/incoming markers!\n');
    %=============================================
    
    %=============================================
    % Surface process
    if WaterLevel > 0
        ind=find(Im==1 & ym>=WaterLevel);
        Im(ind)=2;
        ind=find(Im==2 & ym<WaterLevel);
        Im(ind)=1;

        if SedimentationStyle == 4
            ind=find(Im==2 & ym>WaterLevel);        
            if round(time/SecYear/SedChange)<time/SecYear/SedChange && length(SediMarker)==2
                ii=1;
            else
                ii=2;
            end
            Im(ind) = SediMarker(ii);
            
            ind=surface_y>WaterLevel;
            surface_y(ind)=WaterLevel;
        end
    end
    if  Surface_diffusion == 1
        surface_calculation;
        toc, fprintf('for surface process!\n');
    end
    %=============================================


    fprintf('dt   = %d years\n', dt/SecYear)
    fprintf('Time = %d years\n\n', time/SecYear)

    toc, fprintf('for complete timestep\n===============\nMarknum: %d\n===============\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n', marknum);    

    %======================== SAVE OUTPUT FILE ============================
    if mod(timestep,create_output)==0 || timestep==2
        
        fname = ['Compression_',num2str(1e6+timestep),'.mat'];
        cd OUTPUT
        save([char(fname)]   ,'Vx2d_s','Vy2d_s','E2nd_s','T2nd_s','xm','Txx','Txy',...
            'Temp','Tm',...
            'ym','Im','eta_s','Pm','P2d','strainm','work','workm','Lx','Ly','rho_s','strain','lambda_s',...
            'nx','ny','time','dt','SecYear','x','y','surface_y','iterations')
        cd ..
    end

    % Save breakpoint file
    if mod(timestep,create_breakpoint)==0
        fname = ['Breakpoint.mat'];
        save(char(fname))
%         break
    end
end
