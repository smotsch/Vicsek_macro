#! /usr/bin/octave -qf

%%
%%-- Script (executable) pour afficher le résultat d'une simulation de MacroVic_2D
%%

%% Parameters visualisation
lengthArrow = 5;		# the arrows
headSize    = .1;

shouldPlotEnd = 1;
jumpTime      = 1;
jumpSpace     = 1;
saveVideo     = true;
% parameters color
z_min = 0;
z_max = 2;


%%------------------ 0.1) Read parameters ------------------%%
%%----------------------------------------------------------%%
fid = fopen('../bin/PARAMETER_2D.txt');
C   = textscan(fid, '%s','delimiter', '\n');
Lx   = str2num(C{1}{11});
Ly   = str2num(C{1}{12});
dx   = str2num(C{1}{13});
dy   = str2num(C{1}{14});
Time = str2num(C{1}{16});
dt   = str2num(C{1}{17});
fclose(fid);

nCellx = floor(Lx/dx+.5);
nCelly = floor(Ly/dy+.5);
dxy = min(dx,dy);



%%------------------ 0.2) Initialisation ------------------%%
%%---------------------------------------------------------%%
%% load data
function matrixData=loadBinary2D(nameFile,numberRow)
%% Read the binary data in the file 'nameFile' in 'l'.
%% We also need to precise the number of rows (numberRow).
    fid = fopen(nameFile, 'r');
    fseek(fid, 4, 'bof');
    matrixData = fread(fid,[numberRow,Inf],'double');
    fclose(fid);
endfunction
%% init
nTime = floor(Time/dt+.5);
[X,Y] = meshgrid(dx*(0:nCellx),dy*(0:nCelly));
[X_conv,Y_conv] = meshgrid(dx*((1:nCellx)-.5),dy*((1:nCelly)-.5));
lengthArrow = .8*max(dx,dy);
%% colormap hot
l_red   = [ones(1,38)  , (26:-1:1)/26];
l_green = [ones(1,13)   , (26:-1:1)/26 , zeros(1,25)];
l_blue  = [(12:-1:1)/12 , zeros(1,52)];
A_hot   = [l_red',l_green',l_blue'];

%% on ajuste
if (jumpTime==-1) jumpTime = nTime; end


%%---------------------------------------------------------------%%
%%------------------------  la boucle  --------------------------%%
%%---------------------------------------------------------------%%

for iTime=0:jumpTime:nTime

    %---------------    A) load data     ---------------%
    %---------------------------------------------------%
    rho   = loadBinary2D(['../data/rho_',   num2str(iTime,'%09d'), ".dat"],nCellx);
    theta = loadBinary2D(['../data/theta_', num2str(iTime,'%09d'), ".dat"],nCellx);

    %---------------    B) plot data     ---------------%
    %---------------------------------------------------%

    hold on
    %------ B.1)  Density
    % on evite la saturation
    rho = rho - (rho>z_max).*(rho-z_max);
    imagesc(dx*(1:nCellx)-dx/2,dy*(1:nCelly)-dy/2,rho');
    %% deco      
    set(gca,'YDir','normal')
    axis([0 Lx 0 Ly],"equal")
    colormap(A_hot);
    caxis([z_min z_max])
    c = colorbar;
    %set(c,'yticklabel',"");
    %------ B.2)  Velocity
    % little trick
    Kernel2 = .25*[1 1;1 1];
    f_u = rho .* cos(theta);
    f_v = rho .* sin(theta);
    rho_conv = conv2(rho,Kernel2);
    u_conv = conv2(f_u,Kernel2)./rho_conv;
    v_conv = conv2(f_v,Kernel2)./rho_conv;
    % here it is...
    if (jumpSpace>1)
        quiver(X_conv',Y_conv',lengthArrow*u_conv(2:(nCellxJump+1),2:(nCellyJump+1)),...
               lengthArrow*v_conv(2:(nCellxJump+1),2:(nCellyJump+1)))
    else
        quiver(X_conv',Y_conv',lengthArrow*u_conv(2:(nCellx+1),2:(nCelly+1)),...
               lengthArrow*v_conv(2:(nCellx+1),2:(nCelly+1)))
    end
    hold off
    % deco
    title(["Density and velocity at  t = ",num2str(iTime*dt,"%10.2f")],'fontsize',14)
    xlabel("x",'fontsize',20)
    ylabel("y",'fontsize',20)
    %legend('off')
    %axis off;

    if (saveVideo)
        l = ['images/rho_uv_' , num2str(iTime,'%09d') , '.png'];
        print(sprintf(l),'-F:Helvetica:15','-tight');
    else
        pause(.01)
    end
    
    % small break
    pause(.01)

    % total mass
    %printf("Total mass: %2.6f\n",totalMass)
    
end

%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%

break

if (saveVideo)
    %% we make a movie
    timeNow   = clock;
    extension = [num2str(timeNow(1),"%04d"),...
    num2str(timeNow(2),"%02d"),...
    num2str(timeNow(3),"%02d"),"_",...
    num2str(timeNow(4),"%02d"),"h",...
    num2str(timeNow(5),"%02d")];
    name = [" videos/MacroVic_2D_",extension,".avi"];
    %% we make a movie :
    system(["mencoder ''mf://../images/rho_uv_*.png'' -mf fps=",...
            num2str( 20 )," -o ",name," -ovc  x264"]);
    %% num2str( 10 )," -o ",name," -ovc lavc -lavcopts vcodec=mpeg4"]);
    %% mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc x264 -x264encopts bitrate=1600:pass=1:subq=5
end

%% system(["mencoder ''mf://images/dens_velo*.png'' -mf fps=",num2str( 10 )," -o ",name," -ovc lavc -lavcopts vcodec=mpeg4"]);
