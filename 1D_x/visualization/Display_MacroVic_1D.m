#! /usr/bin/octave -qf

%%
%%-- Script (executable) pour afficher le résultat d'une simulation de MacroVic_1D
%%

%%-- Parameters display
choicePlot = 1;			% 1: ρ,u, 2: ρ,θ, 3: ρ,u,v
jumpTime   = -1;
shouldSave = 0;
%% parameters color
yMin = -.5;
yMax = 2.5;


%%------------------ 0.1) Read parameters ------------------%%
%%----------------------------------------------------------%%
fid = fopen("../bin/PARAMETER_1D.txt");
for i=1:10, temp = fgetl(fid); endfor
Lx = str2num(fgetl(fid));
for i=1:2, temp = fgetl(fid); endfor
Time = str2num(fgetl(fid));
dt   = str2num(fgetl(fid));
fclose(fid);


%%------------------ 0.2) Initialisation ------------------%%
%%---------------------------------------------------------%%
## load data
function l=loadBinary(nameFile)
  ## read the binary data in the file 'nameFile' in 'l'
  fid = fopen(nameFile, 'r');
  fseek(fid, 4, 'bof');
  l = fread(fid, 'double')';
  fclose(fid);
endfunction
## init
nTime = floor(Time/dt+.5);
rho   = loadBinary(["../data/rho_",   num2str(0,"%09d"), ".dat"]);
theta = loadBinary(["../data/theta_", num2str(0,"%09d"), ".dat"]);
n_x = length(rho)-1;
dx  = Lx/n_x;
xx  = dx*(0:n_x);
%% graphic
if (shouldSave==1)
  figure("visible","off")
  system("rm -rf ../images/*")
else
  figure("visible","on")
endif
clf

%% ajustement de dernière minute...
if (jumpTime==-1) jumpTime = nTime; endif
if (choicePlot==3) l = 1:jumpTheta:(n_x+1);  endif


%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%
%%------------------------  la boucle  --------------------------%%
%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%

for iTime=0:jumpTime:nTime

  ##---------------    A) load data     ---------------##
  ##---------------------------------------------------##
  %%----  is there a file at this time ?
  Bool = file_in_path("../data",["rho_", num2str(iTime,"%09d"), ".dat"]);
  if (length(Bool)==0) break endif
  %% load file
  rho   = loadBinary(["../data/rho_",   num2str(iTime,"%09d"), ".dat"]);
  theta = loadBinary(["../data/theta_", num2str(iTime,"%09d"), ".dat"]);

  ##---------------    B) plot data     ---------------##
  ##---------------------------------------------------##
  switch (choicePlot)
    case(1)
      ## la vitesse u ou cos(theta)
      plot(xx,rho,'b','linewidth',3,xx,cos(theta),'g','linewidth',3);
    case(2)
      ## l'angle      
      plot(xx,rho,'b','linewidth',3,xx,theta,'g','linewidth',3);
    case(3)
      ## le vecteur (cos(theta),sin(theta))
      clf
      hold on
      plot(xx,rho,'b',0,0,'k')
      quiver(xx(l),zeros(1,length(l)),1*cos(theta(l)),1*sin(theta(l)),0,'k')
      hold off
  endswitch
  %% deco
  title(["t = ",num2str(iTime*dt,"%10.2f")],'fontsize',14)
  xlabel("x",'fontsize',14)
  axis([0 Lx yMin yMax],"normal");
  switch (choicePlot)
    case(1)
      legend("rho","cos(theta)")
    case(2)
      legend("rho","theta")
    case(3)
      legend("rho","vec(theta)")
  endswitch

  %% save in a png file
  if shouldSave==1
    l = ["../images/rhoTheta_",num2str(iTime,"%09d"),".png"];
    print(sprintf(l));
  endif

  %% small break
  pause(.01)

endfor

%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%

if (shouldSave==1)
  %% we make a movie :
  timeNow   = clock;
  extension = [num2str(timeNow(1),"%04d"),...
	       num2str(timeNow(2),"%02d"),...
	       num2str(timeNow(3),"%02d"),"_",...
	       num2str(timeNow(4),"%02d"),"h",...
	       num2str(timeNow(5),"%02d")];
  name = ["../videos/MacroVic_1D_",extension,".avi"];
  ## the command  
  system(["mencoder ''mf://../images/rhoTheta_*.png'' -mf fps=",...
	  num2str( floor(1/dt/jumpTime) )," -o ",name," -ovc lavc -lavcopts vcodec=mpeg4"]);
else
  pause
endif


break

## save
intXmv = xx;
rhoMV_T4   = rho;
thetaMV_T4 = theta;

save solDC_T4_macro.mat intXmv rhoMV_T4 thetaMV_T4
