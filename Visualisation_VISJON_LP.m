%% Visualisation VISJON-LP

clear all, clf
for num = 1200;
    load Colormaps.mat
    
    cd OUTPUT
    load (['Compression_',num2str(1e6+num),'.mat']);
    cd ..
    
    figure(33), subplot(211)
    
        for i=1:20
            Phase1  =   find(Im == i);
            marksize = 1;
            hold on
            plot(xm(Phase1),ym(Phase1),'.','Color',Colormaps(i,:),'MarkerSize',marksize)
        end
    shading interp
    axis image, axis ij
    title(['Time: ',num2str((time)/SecYear/1e6),' Ma'])
%     plot(x2d_b,y2d_b,'k',x2d_b',y2d_b','k')
%     quiver(x2d_b,y2d_b,Vx2db,Vy2db,'r','LineWidth',1.5)
    hold on
    
    E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
    eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
    T2nd_s_2d   =   reshape(T2nd_s,ny+1,nx+1);
    rho_s_2d    =   reshape(rho_s,ny+1,nx+1);
    Temp_2d     =   reshape(Temp,ny+1,nx+1);
    
    contour(x,y,(Temp_2d(1:end-1,1:end-1)-273),[100:200:1300],'w')
    
    hold off
    
    figure(33), subplot(212)
    
    ind = find(rho_s_2d < 200);
    E2nd_s_2d(ind)  = 1e-14;   
    T2nd_s_2d(ind)  = 1e20;   
    eta_s_2d(ind)   = 1e21;   
    
%     subplot(311)
    pcolor(x,y,log10(eta_s_2d(1:end-1,1:end-1))), caxis([17 25])
    shading interp
%     colorbar
    axis image, axis ij
    title('E2nd ')
    hold on
    fill([0 x Lx 0],[0 surface_y 0 0],'w','EdgeColor','w','FaceColor','w','LineWidth',0.001)
    
%         
%      comp='norm_08'
%      print('-r900', '-djpeg', comp)
% 
%     subplot(312)
%     pcolor(x,y,log10(eta_s_2d(1:end-1,1:end-1))), caxis([17 24])
%     shading interp
%     colorbar
%     axis image, axis ij
%     title(['\eta_{s}   ',    num2str(dt/SecYear)])
%     
%     %         subplot(223)
%     %         pcolor(x_Vx(1:end-1),y_Vy(1:end-1),Txy(1:end-1,1:end-1))% caxis([-15 -12])
%     %         shading interp
%     %         colorbar
%     %         axis image, axis ij
%     %         title(['Txy'])
%     
%     subplot(313)
%     pcolor(x,y,T2nd_s_2d(1:end-1,1:end-1)), caxis([0 2e8])
%     shading interp
%     colorbar
%     axis image, axis ij
%     title(['T2nd'])
    
    % plot phases (not markers)
%     multi = 2;
%     
%     I_dx  =   Lx/(nx-1)/multi;
%     I_x   =   0:I_dx:Lx+(I_dx);
%     I_dy  =   Ly/(ny-1)/multi;
%     I_y   =   0:I_dy:Ly+(I_dy);
%     
%     II  =   zeros(length(I_y),length(I_x));
%     wts =   zeros(length(I_y),length(I_x));
%     
%     marknum = length(xm);
% 
%     for marker = 1:marknum
%         %             if xm(marker) >= 0 && xm(marker) <= Lx && ym(marker) >= 0 && ym(marker) <=Ly
%         % define indexes
%         I_j = fix(xm(marker)/I_dx)+1;
%         I_i = fix(ym(marker)/I_dy)+1;
%         if I_j<1
%             I_j = 1;
%         end
%         if I_j>nx*multi
%             I_j = nx*multi;
%         end
%         if I_i<1
%             I_i = 1;
%         end
%         if I_i>ny*multi
%             I_i = ny*multi;
%         end
%         
%         % updating eta_s (i,j)
%         if xm(marker)-I_x(I_j) <= I_dx/2 && ym(marker)-I_y(I_i) <= I_dy/2
%             wtijm               =   (1-(xm(marker)-I_x(I_j))/I_dx)*(1-(ym(marker)-I_y(I_i))/I_dy);
%             II(I_i,I_j)         =   II(I_i,I_j)+Im(marker)*wtijm;
%             wts(I_i,I_j)        =   wts(I_i,I_j) + wtijm;
%         end
%         % (i,j+1)
%         if xm(marker)-I_x(I_j) >= I_dx/2 && ym(marker)-I_y(I_i) <= I_dy/2
%             wtij1m              =   ((xm(marker)-I_x(I_j))/I_dx)*(1-(ym(marker)-I_y(I_i))/I_dy);
%             II(I_i,I_j+1)       =   II(I_i,I_j+1)+Im(marker)*wtij1m;
%             wts(I_i,I_j+1)      =   wts(I_i,I_j+1) + wtij1m;
%         end
%         % (i+1,j)
%         if xm(marker)-I_x(I_j) <= I_dx/2 && ym(marker)-I_y(I_i) >= I_dy/2
%             wti1jm          =   (1-(xm(marker)-I_x(I_j))/I_dx)*((ym(marker)-I_y(I_i))/I_dy);
%             II(I_i+1,I_j)       =   II(I_i+1,I_j)+Im(marker)*wti1jm;
%             wts(I_i+1,I_j)      =   wts(I_i+1,I_j) + wti1jm;
%         end
%         % (i+1,j+1)
%         if xm(marker)-I_x(I_j) >= I_dx/2 && ym(marker)-I_y(I_i) >= I_dy/2
%             wti1j1m             =   ((xm(marker)-I_x(I_j))/I_dx)*((ym(marker)-I_y(I_i))/I_dy);
%             II(I_i+1,I_j+1)     =   II(I_i+1,I_j+1)+Im(marker)*wti1j1m;
%             wts(I_i+1,I_j+1)    =   wts(I_i+1,I_j+1) + wti1j1m;
%         end
%     end
%     
%     % Divide by weights (i,j)
%     for I_i = 1:(ny-1)*multi+1
%         for I_j = 1:(nx-1)*multi+1
%             if wts(I_i,I_j) > 0
%                 II(I_i,I_j)         = II(I_i,I_j)/wts(I_i,I_j);
%             else
%                 II(I_i,I_j)         = 10;
%             end
%         end
%     end
%     
%     figure(33), subplot(212)
%     
%     pcolor(I_x(1:end-1),I_y(1:end-1),II(1:end-1,1:end-1))% caxis([18 24])
%     shading interp
%     colormap('Colormaps')
%     axis image, axis ij
%     caxis([1 20])

%     cd ..
end