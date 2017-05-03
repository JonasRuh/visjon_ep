  
%=====================================================
% Calculate visco-elastic viscosity and stresses
% Pressure scaling variable
ind     =   eta_s ~= 0;
etasmin =   min(eta_s(ind));
pscale  =   etasmin/dx;
kcont   =   2*etasmin/(dx+dy);
kbond   =   4*etasmin/(dx+dy)^2;

% Find equation positions for x-stokes
if timestep == 1
    matrix_x_eq                         = rho_vx;
    matrix_x_eq(1:ny+1)                 = 0;
    matrix_x_eq(end-(2*(ny+1))+1:end)   = 0;
    
    % Fill up matrix
    kx_eq       =   find(matrix_x_eq ~= 0);
    kx_leftBC   =   1:(ny+1);
    kx_rightBC  =   (ny+1)*(nx-1)+1:(ny+1)*(nx+1);
    kx_topBC    =   (ny+1)+1:(ny+1):(ny+1)*(nx-2)+1;
    kx_bottomBC =   2*(ny+1):(ny+1):(ny+1)*(nx-1);
end

num1 = 1;
num2 = length(kx_eq);

% x-Stokes
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq-(ny+1);
ss(num1:num2)   =   (dt*eta_p(kx_eq)        .*mu_p(kx_eq))       ./(dx^2*(dt*mu_p(kx_eq)          +eta_p(kx_eq)));                % Vx1
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq-1;
ss(num1:num2)   =   (dt*eta_s(kx_eq-1)      .*mu_s(kx_eq-1))     ./(dy^2*(dt*mu_s(kx_eq-1)        +eta_s(kx_eq-1)));        % Vx2
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq;
ss(num1:num2)   =  -(dt*eta_p(kx_eq+(ny+1)) .*mu_p(kx_eq+(ny+1)))./(dx^2*(dt*mu_p(kx_eq+(ny+1))   +eta_p(kx_eq+(ny+1))))...
    -(dt*eta_p(kx_eq)        .*mu_p(kx_eq))       ./(dx^2*(dt*mu_p(kx_eq)          +eta_p(kx_eq)))...
    -(dt*eta_s(kx_eq)        .*mu_s(kx_eq))       ./(dy^2*(dt*mu_s(kx_eq)          +eta_s(kx_eq)))...
    -(dt*eta_s(kx_eq-1)      .*mu_s(kx_eq-1))     ./(dy^2*(dt*mu_s(kx_eq-1)        +eta_s(kx_eq-1)));        % Vx3
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+1;
ss(num1:num2)   =   (dt*eta_s(kx_eq)        .*mu_s(kx_eq))       ./(dy^2*(dt*mu_s(kx_eq)          +eta_s(kx_eq)));                % Vx4
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+(ny+1);
ss(num1:num2)   =   (dt*eta_p(kx_eq+(ny+1)) .*mu_p(kx_eq+(ny+1)))./(dx^2*(dt*mu_p(kx_eq+(ny+1))   +eta_p(kx_eq+(ny+1))));        % Vx5
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+(ny+1)*(nx+1)-1;
ss(num1:num2)   =   (dt*eta_s(kx_eq-1)      .*mu_s(kx_eq-1))     ./(dy*dx*(dt*mu_s(kx_eq-1)       +eta_s(kx_eq-1)))...
    -(dt*eta_p(kx_eq)        .*mu_p(kx_eq))       ./(dy*dx*(dt*mu_p(kx_eq)         +eta_p(kx_eq)));     % Vy1
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+(ny+1)*(nx+1);
ss(num1:num2)   =  -(dt*eta_s(kx_eq)        .*mu_s(kx_eq))       ./(dy*dx*(dt*mu_s(kx_eq)         +eta_s(kx_eq)))...
    +(dt*eta_p(kx_eq)        .*mu_p(kx_eq))       ./(dy*dx*(dt*mu_p(kx_eq)         +eta_p(kx_eq)));             % Vy2
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+(ny+1)*(nx+1)+(ny+1)-1;
ss(num1:num2)   =  -(dt*eta_s(kx_eq-1)      .*mu_s(kx_eq-1))     ./(dy*dx*(dt*mu_s(kx_eq-1)       +eta_s(kx_eq-1)))...
    +(dt*eta_p(kx_eq+(ny+1)) .*mu_p(kx_eq+(ny+1)))./(dx*dy*(dt*mu_p(kx_eq+(ny+1))  +eta_p(kx_eq+(ny+1))));    % Vy3
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+(ny+1)*(nx+1)+(ny+1);
ss(num1:num2)   =   (dt*eta_s(kx_eq)        .*mu_s(kx_eq))       ./(dy*dx*(dt*mu_s(kx_eq)         +eta_s(kx_eq)))...
    -(dt*eta_p(kx_eq+(ny+1)) .*mu_p(kx_eq+(ny+1)))./(dx*dy*(dt*mu_p(kx_eq+(ny+1))  +eta_p(kx_eq+(ny+1))));        % Vy4
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+2*(ny+1)*(nx+1);
ss(num1:num2)   =   kcont/dx;   % P1
num1 = num2+1; num2 = num2+length(kx_eq);
ii(num1:num2)   =   kx_eq;
jj(num1:num2)   =   kx_eq+2*(ny+1)*(nx+1)+(ny+1);
ss(num1:num2)   =  -kcont/dx;  % P2
num1 = num2+1; num2 = num2+length(kx_leftBC);
R(kx_eq)        =  -((eta_s(kx_eq)          .*Toxy(kx_eq))          ./(dt*mu_s(kx_eq)       +eta_s(kx_eq))          -(eta_s(kx_eq-1).*Toxy(kx_eq-1))./(dt*mu_s(kx_eq-1) +eta_s(kx_eq-1)))   /dy...
    -((eta_p(kx_eq+(ny+1))   .*Toxx(kx_eq+(ny+1)))   ./(dt*mu_p(kx_eq+(ny+1))+eta_p(kx_eq+(ny+1)))   -(eta_p(kx_eq)  .*Toxx(kx_eq))  ./(dt*mu_p(kx_eq)   +eta_p(kx_eq)))     /dx...
    -rho_vy(kx_eq)*gravity_x;
% x-stokes BC
% Left boundary
ii(num1:num2)   =   kx_leftBC;
jj(num1:num2)   =   kx_leftBC;
ss(num1:num2)   =   1*kbond;
R(kx_leftBC)    =   BC_left./SecYear*kbond;
num1 = num2+1; num2 = num2+length(kx_rightBC);
% Right boundary
ii(num1:num2)   =   kx_rightBC;
jj(num1:num2)   =   kx_rightBC;
ss(num1:num2)   =   1*kbond;
R(kx_rightBC)   =   BC_right./SecYear*kbond;
num1 = num2+1; num2 = num2+length(kx_topBC);
% Upper boundary: free slip
ii(num1:num2)   =   kx_topBC;
jj(num1:num2)   =   kx_topBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(kx_topBC);
ii(num1:num2)   =   kx_topBC;
jj(num1:num2)   =   kx_topBC+1;
if strcmp('freeslip',top_BC)==1
    ss(num1:num2)   =  -1*kbond;
    R(kx_topBC)     =   0;
elseif strcmp('noslip',top_BC)==1
    ss(num1:num2)   =   1*kbond;
    R(kx_topBC)     =   0;
elseif strcmp('mixed',top_BC)==1
    ss(num1:num2)   =   top_mix.*kbond;
    R(kx_topBC)     =   0;
end
num1 = num2+1; num2 = num2+length(kx_bottomBC);
% Lower boundary: velocity
ii(num1:num2)   =   kx_bottomBC;
jj(num1:num2)   =   kx_bottomBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(kx_bottomBC);
ii(num1:num2)   =   kx_bottomBC;
jj(num1:num2)   =   kx_bottomBC-1;
if strcmp('freeslip',bottom_BC)==1 || strcmp('external',bottom_BC)==1
    ss(num1:num2)   =  -1*kbond;
    R(kx_bottomBC)  =   0;
elseif strcmp('velocity',bottom_BC)==1
    ss(num1:num2)   =   1*kbond;
    R(kx_bottomBC)  =   2*BC_bottom/SecYear*kbond;
elseif strcmp('noslip',bottom_BC)==1
    ss(num1:num2)   =   1*kbond;
    R(kx_bottomBC)  =   0;
elseif strcmp('mixed',bottom_BC)==1
    ss(num1:num2)   =   bottom_mix(1,:).*kbond;
    R(kx_bottomBC)  =   0;
    if size(bottom_mix,1)>1
        R(kx_bottomBC)  =   2.*bottom_mix(2,:)./SecYear*kbond;
    end        
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Find equation positions for y-stokes
if timestep == 1
    matrix_y_eq                         = rho_vy;
    matrix_y_eq(1:ny+1:(ny+1)*nx+1)     = 0;
    matrix_y_eq(ny:(ny+1):(ny+1)*nx)    = 0;
    
    % Fill up matrix
    ky_eq       =   find(matrix_y_eq ~= 0);
    ky_ind      =   ky_eq+(ny+1)*(nx+1);
    ky_leftBC   =   [2:ny-1]+(ny+1)*(nx+1);
    ky_rightBC  =   [(ny+1)*(nx)+2:(ny+1)*(nx+1)-2]+(ny+1)*(nx+1);
    ky_topBC    =   [1:(ny+1):(ny+1)*(nx)+1]+(ny+1)*(nx+1);
    ky_bottomBC =   [(ny):(ny+1):(ny+1)*(nx+1)-1 (ny+1):(ny+1):(ny+1)*(nx+1)]+(ny+1)*(nx+1);
end

num1 = num2+1; num2 = num2+length(ky_eq);
% y-Stokes
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind-(ny+1);
ss(num1:num2)   =   (dt*eta_s(ky_eq-(ny+1)) .*mu_s(ky_eq-(ny+1)))./(dx^2*(dt*mu_s(ky_eq-(ny+1))   +eta_s(ky_eq-(ny+1))));       % Vy1
num1 = num2+1; num2 = num2+length(ky_eq);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind-1;
ss(num1:num2)   =   (dt*eta_p(ky_eq)        .*mu_p(ky_eq))       ./(dy^2*(dt*mu_p(ky_eq)          +eta_p(ky_eq)));               % Vy2
num1 = num2+1; num2 = num2+length(ky_eq);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind;
ss(num1:num2)   =  -(dt*eta_p(ky_eq+1)      .*mu_p(ky_eq+1))     ./(dy^2*(dt*mu_p(ky_eq+1)        +eta_p(ky_eq+1)))...
    -(dt*eta_p(ky_eq)        .*mu_p(ky_eq))       ./(dy^2*(dt*mu_p(ky_eq)          +eta_p(ky_eq)))...
    -(dt*eta_s(ky_eq-(ny+1)) .*mu_s(ky_eq-(ny+1)))./(dx^2*(dt*mu_s(ky_eq-(ny+1))   +eta_s(ky_eq-(ny+1))))...
    -(dt*eta_s(ky_eq)        .*mu_s(ky_eq))       ./(dx^2*(dt*mu_s(ky_eq)          +eta_s(ky_eq)));
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind+1;
ss(num1:num2)   =   (dt*eta_p(ky_eq+1)      .*mu_p(ky_eq+1))     ./(dy^2*(dt*mu_p(ky_eq+1)        +eta_p(ky_eq+1)));       % Vy4
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind+(ny+1);
ss(num1:num2)   =   (dt*eta_s(ky_eq)        .*mu_s(ky_eq))       ./(dx^2*(dt*mu_s(ky_eq)          +eta_s(ky_eq)));               % Vy5
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind-(ny+1)*(nx+1)-(ny+1);
ss(num1:num2)   =   (dt*eta_s(ky_eq-(ny+1)) .*mu_s(ky_eq-(ny+1)))./(dy*dx*(dt*mu_s(ky_eq-(ny+1))+eta_s(ky_eq-(ny+1))))...
    -(dt*eta_p(ky_eq)        .*mu_p(ky_eq))       ./(dy*dx*(dt*mu_p(ky_eq)  +eta_p(ky_eq)));%-gravity/4*((rho_s(i,j)-rho_s(i,j-1))/dx)*(dt/2);      % Vx1
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind-(ny+1)*(nx+1)-(ny+1)+1;
ss(num1:num2)   =  -(dt*eta_s(ky_eq-(ny+1)) .*mu_s(ky_eq-(ny+1)))./(dy*dx*(dt*mu_s(ky_eq-(ny+1))+eta_s(ky_eq-(ny+1))))...
    +(dt*eta_p(ky_eq+1)      .*mu_p(ky_eq+1))     ./(dx*dy*(dt*mu_p(ky_eq+1)+eta_p(ky_eq+1)));%-gravity/4*((rho_s(i,j)-rho_s(i,j-1))/dx)*(dt/2);              % Vx2
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind-(ny+1)*(nx+1);
ss(num1:num2)   =  -(dt*eta_s(ky_eq)        .*mu_s(ky_eq))       ./(dy*dx*(dt*mu_s(ky_eq)  +eta_s(ky_eq)))...
    +(dt*eta_p(ky_eq)        .*mu_p(ky_eq))       ./(dy*dx*(dt*mu_p(ky_eq)  +eta_p(ky_eq)));%-gravity/4*((rho_s(i,j)-rho_s(i,j-1))/dx)*(dt/2);      % Vx3
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind-(ny+1)*(nx+1)+1;
ss(num1:num2)   =   (dt*eta_s(ky_eq)        .*mu_s(ky_eq))       ./(dy*dx*(dt*mu_s(ky_eq)  +eta_s(ky_eq)))...
    -(dt*eta_p(ky_eq+1)      .*mu_p(ky_eq+1))     ./(dx*dy*(dt*mu_p(ky_eq+1)+eta_p(ky_eq+1)));%-gravity/4*((rho_s(i,j)-rho_s(i,j-1))/dx)*(dt/2);              % Vx4
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind+(ny+1)*(nx+1);
ss(num1:num2)   =   kcont/dy;   % P1
num1 = num2+1; num2 = num2+length(ky_ind);
ii(num1:num2)   =   ky_ind;
jj(num1:num2)   =   ky_ind+(ny+1)*(nx+1)+1;
ss(num1:num2)   =  -kcont/dy;  % P2
num1 = num2+1; num2 = num2+length(ky_leftBC);
R(ky_ind)       =  -((eta_s(ky_eq)  .*Toxy(ky_eq))  ./(dt*mu_s(ky_eq)  +eta_s(ky_eq))  -(eta_s(ky_eq-(ny+1)).*Toxy(ky_eq-(ny+1)))./(dt*mu_s(ky_eq-(ny+1))+eta_s(ky_eq-(ny+1))))/dx...
    -((eta_p(ky_eq+1).*(-Toxx(ky_eq+1)))./(dt*mu_p(ky_eq+1)+eta_p(ky_eq+1))-(eta_p(ky_eq)  .*(-Toxx(ky_eq)))  ./(dt*mu_p(ky_eq)  +eta_p(ky_eq)))  /dy...
    -rho_vy(ky_eq)*gravity_y;

% y-Stokes BC
% Left boundary
ii(num1:num2)   =   ky_leftBC;
jj(num1:num2)   =   ky_leftBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(ky_leftBC);
ii(num1:num2)   =   ky_leftBC;
jj(num1:num2)   =   ky_leftBC+(ny+1);
if strcmp('freeslip',left_BC)==1
    ss(num1:num2)   =  -1*kbond;
elseif strcmp('noslip',left_BC)==1
    ss(num1:num2)   =   1*kbond;
elseif strcmp('mixed',left_BC)==1
    ss(num1:num2)   =   left_mix.*kbond;
end
num1 = num2+1; num2 = num2+length(ky_rightBC);
R(ky_leftBC)    =   0;
% Right boundary
ii(num1:num2)   =   ky_rightBC;
jj(num1:num2)   =   ky_rightBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(ky_rightBC);
ii(num1:num2)   =   ky_rightBC;
jj(num1:num2)   =   ky_rightBC-(ny+1);
if strcmp('freeslip',right_BC)==1
    ss(num1:num2)   =  -1*kbond;
elseif strcmp('noslip',right_BC)==1
    ss(num1:num2)   =   1*kbond;
elseif strcmp('mixed',right_BC)==1
    ss(num1:num2)   =   right_mix.*kbond;
end
num1 = num2+1; num2 = num2+length(ky_topBC);
R(ky_rightBC)   =   0;
% Upper boundary: velocity
ii(num1:num2)   =   ky_topBC;
jj(num1:num2)   =   ky_topBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(ky_bottomBC);
R(ky_topBC)     =   BC_top./SecYear*kbond;
% Lower boundary: zero velocity
ii(num1:num2)   =   ky_bottomBC;
jj(num1:num2)   =   ky_bottomBC;
ss(num1:num2)   =   1*kbond;
if strcmp('freeslip',bottom_BC)==1
    R(ky_bottomBC)  =   BC_bottom./SecYear*kbond;
elseif strcmp('velocity',bottom_BC)==1
    R(ky_bottomBC)  =   0;
elseif strcmp('external',bottom_BC)==1
    num1 = num2+1; num2 = num2+length(ky_bottomBC);
    ii(num1:num2)   =   ky_bottomBC;
    jj(num1:num2)   =   ky_bottomBC-1;
    ss(num1:num2)   =  -Ext_BC*kbond;
    R(ky_bottomBC)  =   BC_bottom./SecYear*kbond;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Find equation positions for continuity equation
if timestep == 1
    matrix_p_eq                     = eta_p;
    matrix_p_eq(ny+1+2)             = 0;
    matrix_p_eq(ny+1+ny)            = 0;
    matrix_p_eq((ny+1)*(nx-1)+2)    = 0;
    matrix_p_eq((ny+1)*(nx-1)+ny)   = 0;
    matrix_p_eq((ny+1)*(p_node-1)+2)= 0;
    
    % Fill up matrix
    kp_ind      =   find(matrix_p_eq ~= 0)+2*(ny+1)*(nx+1);
    kp_sideBC   =   find(eta_p == 0)+2*(ny+1)*(nx+1);
    kp_leftBC   =   [(ny+1+2) (ny+1+ny)]+2*(ny+1)*(nx+1);
    kp_rightBC  =   [(ny+1)*(nx-1)+2 (ny+1)*(nx-1)+ny]+2*(ny+1)*(nx+1);
    kp_pBC      =   [(ny+1)*(p_node-1)+2]+2*(ny+1)*(nx+1);
end

num1 = num2+1; num2 = num2+length(kp_ind);
% Continuity
ii(num1:num2)   =   kp_ind;
jj(num1:num2)   =   kp_ind-2*(ny+1)*(nx+1);
ss(num1:num2)   =   1*kcont/dx;
num1 = num2+1; num2 = num2+length(kp_ind);
ii(num1:num2)   =   kp_ind;
jj(num1:num2)   =   kp_ind-2*(ny+1)*(nx+1)-(ny+1);
ss(num1:num2)   =  -1*kcont/dx;
num1 = num2+1; num2 = num2+length(kp_ind);
ii(num1:num2)   =   kp_ind;
jj(num1:num2)   =   kp_ind-(ny+1)*(nx+1);
ss(num1:num2)   =   1*kcont/dy;
num1 = num2+1; num2 = num2+length(kp_ind);
ii(num1:num2)   =   kp_ind;
jj(num1:num2)   =   kp_ind-(ny+1)*(nx+1)-1;
ss(num1:num2)   =  -1*kcont/dy;
num1 = num2+1; num2 = num2+length(kp_sideBC);
R(kp_ind)       =   0;  %-(reallog(rho_new(i,j))-reallog(rho_old(i,j))/dt);

% continuity BC
ii(num1:num2)   =   kp_sideBC;
jj(num1:num2)   =   kp_sideBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(kp_leftBC);
R(kp_sideBC)    =   0;
% Upper and lower left corners dP/dx=0 => P(i,j)-P(i,j+1)=0
ii(num1:num2)   =   kp_leftBC;
jj(num1:num2)   =   kp_leftBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(kp_leftBC);
ii(num1:num2)   =   kp_leftBC;
jj(num1:num2)   =   kp_leftBC+(ny+1);
ss(num1:num2)   =  -1*kbond;
num1 = num2+1; num2 = num2+length(kp_rightBC);
R(kp_leftBC)    =   0;
% Upper and lower right corners dP/dx=0 => P(i,j)-P(i,j-1)=0
ii(num1:num2)   =   kp_rightBC;
jj(num1:num2)   =   kp_rightBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(kp_rightBC);
ii(num1:num2)   =   kp_rightBC;
jj(num1:num2)   =   kp_rightBC-(ny+1);
ss(num1:num2)   =  -1*kbond;
num1 = num2+1; num2 = num2+length(kp_pBC);
R(kp_rightBC)       =   0;
% Pressure point
ii(num1:num2)   =   kp_pBC;
jj(num1:num2)   =   kp_pBC;
ss(num1:num2)   =   1*kbond;
num1 = num2+1; num2 = num2+length(kp_pBC);
R(kp_pBC)       =   0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ind = find(ii<=(nx+1)*(ny+1)*2 & jj<=(nx+1)*(ny+1)*2);
% L1 = sparse(ii(ind),jj(ind),ss(ind),(nx+1)*(ny+1)*2,(nx+1)*(ny+1)*2);
% ind = find(ii>(nx+1)*(ny+1)*2 & jj<=(nx+1)*(ny+1)*2);
% L2 = sparse(ii(ind)-(nx+1)*(ny+1)*2,jj(ind),ss(ind),(nx+1)*(ny+1),(nx+1)*(ny+1)*2);
% ind = find(ii>(nx+1)*(ny+1)*2 & jj>(nx+1)*(ny+1)*2);
% L3 = sparse(ii(ind)-(nx+1)*(ny+1)*2,jj(ind)-(nx+1)*(ny+1)*2,ss(ind),(nx+1)*(ny+1),(nx+1)*(ny+1));
% 
% L = L1+L2'*(L3^-1)*L2;

L = sparse(ii,jj,ss,(nx+1)*(ny+1)*3,(nx+1)*(ny+1)*3);

toc, fprintf('for filling the matrix!\n')

S = L\R;    % using MATLAB back slash solver

Vx    =   S(1:(ny+1)*(nx+1));
Vy    =   S(1+(ny+1)*(nx+1):2*(ny+1)*(nx+1));
P     =   S(2*(ny+1)*(nx+1)+1:end).*kcont;

Vx2d    =   reshape(Vx,ny+1,nx+1);
Vy2d    =   reshape(Vy,ny+1,nx+1);
P2d     =   reshape(P,ny+1,nx+1);

%===================================================================
% redistributing S

%   1--6--11--16
%   |  |  |    |
%   2--7--12--17
%   |  |  |    |
%   3--8--13--18
%   |  |  |    |
%   4--9--14--19
%   |  |  |    |
%   5-10--15--20

Vx2d_s   =   (Vx2d(1:end-1,1:end-1)+Vx2d(2:end,1:end-1))./2;
Vy2d_s   =   (Vy2d(1:end-1,1:end-1)+Vy2d(1:end-1,2:end))./2;
P2d_p    =   P2d;
