% length of delta x and delta y
dx  =   Lx/(nx-1);
x   =   0:dx:Lx;
dy  =   Ly/(ny-1);
y   =   0:dy:Ly;

% initiate coordinates
x_Vy    =   -dx/2:dx:Lx+dx/2;
y_Vy    =   0:dy:Ly+dy;

x_Vx    =   0:dx:Lx+dx;
y_Vx    =   -dy/2:dy:Ly+dy/2;

x_P     =   -dx/2:dx:Lx+dx/2;
y_P     =   -dy/2:dy:Ly+dy/2;

% % initiate 2D grids
% % Basic points
% [x2d_b,y2d_b]   =   meshgrid(x_Vx(1:end-1),y_Vy(1:end-1));
% % Pressure points
% [x2d_p,y2d_p]   =   meshgrid(x_P,y_P);
% % Velocity points (without ghost points)
% [x2d_vx,y2d_vx]   =   meshgrid(x_Vx(1:end-1),y_Vx);
% [x2d_vy,y2d_vy]   =   meshgrid(x_Vy,y_Vy(1:end-1));

% Number of Markers
nxm     =   mx*(nx-1);
nym     =   my*(ny-1);
marknum =   nxm*nym;

left_income   = zeros(1,ny-1);
right_income  = zeros(1,ny-1);
top_income    = zeros(1,nx-1);
bottom_income = zeros(1,nx-1);

xm      =   zeros(1,marknum);
ym      =   zeros(1,marknum);
etam    =   zeros(1,marknum);
mum     =   zeros(1,marknum);
wm      =   zeros(1,marknum);
rhom0   =   zeros(1,marknum);
phim    =   zeros(1,marknum);
Cm      =   zeros(1,marknum);
Txxm    =   zeros(1,marknum);
Txym    =   zeros(1,marknum);
txxm    =   zeros(1,marknum);
txym    =   zeros(1,marknum);
dtxxm    =   zeros(1,marknum);
dtxym    =   zeros(1,marknum);
Exxm    =   zeros(1,marknum);
Exym    =   zeros(1,marknum);
Exxn    =   zeros(1,marknum);
Exyn    =   zeros(1,marknum);
E2ndm   =   zeros(1,marknum);
T2ndm   =   zeros(1,marknum);
E2ndn   =   zeros(1,marknum);
Pm      =   zeros(1,marknum);
Im      =   zeros(1,marknum);
lambdam =   zeros(1,marknum);
strainm =   zeros(1,marknum);
workm   =   zeros(1,marknum);
alpha   =   zeros(1,marknum);
if Temperature == 1
    if ra_heat == 1
        hrm     =   zeros(1,marknum);
    end
    if adiab_heat == 1 && TPdep_dens == 0
        tam     =   zeros(1,marknum);
    end
    if TPdep_dens == 1
        tam     =   zeros(1,marknum);
        tbm     =   zeros(1,marknum);
    end 
    Tm      =   zeros(1,marknum);
    km      =   zeros(1,marknum);
    cpm     =   zeros(1,marknum);
    rhocpm0 =   zeros(1,marknum);
end
if Powerlaw == 1
    RK = 8.3145;
    nm      =   zeros(1,marknum);
    Qm      =   zeros(1,marknum);
end 
    
eta_s   =   zeros((ny+1)*(nx+1),1);
Toxy    =   zeros((ny+1)*(nx+1),1);
mu_s    =   zeros((ny+1)*(nx+1),1);
w       =   zeros((ny+1)*(nx+1),1);
rho_s   =   zeros((ny+1)*(nx+1),1);
eta_p   =   zeros((ny+1)*(nx+1),1);
Toxx    =   zeros((ny+1)*(nx+1),1);
mu_p    =   zeros((ny+1)*(nx+1),1);
rho_p   =   zeros((ny+1)*(nx+1),1);
rho_vy  =   zeros((ny+1)*(nx+1),1);
rho_vx  =   zeros((ny+1)*(nx+1),1);
strain  =   zeros((ny+1)*(nx+1),1);
work    =   zeros((ny+1)*(nx+1),1);
lambda_s=   zeros((ny+1)*(nx+1),1);
if Temperature == 1
    if ra_heat == 1
        Hr      =   zeros((ny+1)*(nx+1),1);
    end
    if adiab_heat == 1 || TPdep_dens == 1
        Ta      =   zeros((ny+1)*(nx+1),1);
    end 
    k_vx    =   zeros((ny+1)*(nx+1),1);
    k_vy    =   zeros((ny+1)*(nx+1),1);
    Temp    =   zeros((ny+1)*(nx+1),1);
    cp_p    =   zeros((ny+1)*(nx+1),1);
    rhocp_p =   zeros((ny+1)*(nx+1),1);
end    
    
% Matrix formation
R   =   zeros((nx+1)*(ny+1)*3,1);
LT  =   spalloc((nx+1)*(ny+1),(nx+1)*(ny+1),((nx+1)*(ny+1))*11);    % number of temperature coefficients
RT  =   zeros((nx+1)*(ny+1),1);
LS  =   spalloc((nx-1)*surface_nodes+1,(nx-1)*surface_nodes+1,((nx-1)*surface_nodes+1)*5);
RS  =   zeros((nx-1)*surface_nodes+1,1);
% Right hand side
S   =   zeros((nx+1)*(ny+1)*3,1);

iterations=zeros(2,1);

