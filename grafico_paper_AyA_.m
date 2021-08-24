function [bx_nube,by_nube,bz_nube,date_,theta,phi,R,bx_nube_ext,by_nube_ext,bz_nube_ext,date_ext,eval1,eval2,eval3,evec1,evec2,evec3,B_,B,cociente,cocientefrente,cocienteback]=grafico_paper_AyA_(DIR_WORK,bx_gse_ext_,by_gse_ext_,bz_gse_ext_,date_ext_,Vsw,initial_date,end_date,graf,gap,B__,...
kind_mv,initial_date_graf,end_date_graf,valordelafechaenfrente,valordelafechaenback,jj,flag_sat)

% Orientation of main axes of a cylindrical magnetic cloud
% from b_gse components (using a minimum variance method) and rotate B accordingly.
% [bx_nube,by_nube,bz_nube,date_,theta,phi,R,...
% bx_nube_ext,by_nube_ext,bz_nube_ext,date_ext,eval1,eval2,eval3,evec1,evec2,evec3,B_]...
% =min_var(bx_gse_ext_,by_gse_ext_,bz_gse_ext_,date_ext_,Vsw,initial_date,end_date,graf,gap,B)
%
% Input
% bx_gse,by_gse,bz_gse: Extended time series of B (GSE, in nT) to process with MV method
% date_ext: Time series of the times associated with the field (it has to be in the date format of matlab)
% Vsw: Solar Wind velocity in X_GSE (positive scalar), in Km/sec
% initial_date: Date and Time of the beginning of the cloud (it has to be in the date format of matlab)
% end_date: Date and Time of the end of the cloud (it has to be in the date format of matlab)
% graf: flag to decide if produce or not plots of rotated field (0:no, 1:yes).
% gap: flag value to indicate data gap in b
% B__: absolute value of B series [array].
% kind_mv: key to decide kind of method (1 for MCs, 2 for shocks)
% initial_date_graf: initial date to plot (only useful if graf=1)
% end_date_graf:     end     date to plot (only useful if graf=1)
%
% Output
% bx_nube, by_nube, bz_nube: Magnetic field (in nT) in the coordinates system oriented as the cloud.
%                             (z_n parallel to the clous axis, x_n as r outbound)
% theta: angle (deg) between z_cloud and ecliptic (-90,+90), theta=90 corresponds to north
% phi: angle (deg) between the proj of z_n on ecliptic and x_gse, counterclock (0,360)
% R: Radius of the cloud from time lapse, Vsw, and orientation (in AU).
% bx_nube_ext, by_nube_ext, bz_nube_ext: Idem bi_nube, but to the extended range of time.
% date_ext: extended range of time.
% eval1,eval2,eval3: eigen values ordered from min to max (from normalized MV)
% evec1,evec2,evec3: eigen vectors (in GSE, same order as eigenvalues from normalized MV)
% B_: absolute value of B series, inside the flux rope [array].
%
% v1.2, by S. Dasso, Dec 20, 2007.
% Instituto de Astronom� y F�ica del Espacio (IAFE, CONICET-UBA)
% Departamento de F�ica, FCEN (UBA).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Definitions:
AU=1.49597870691e8 %AU=1.49598e8; %1AU=1.5e8km
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove data gaps for bx_gse_ext
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon=1e-4;
if gap>0,
   cond_good=find(bx_gse_ext_<gap-epsilon & by_gse_ext_<gap-epsilon & bz_gse_ext_<gap-epsilon);
elseif gap<0,
   cond_good=find(bx_gse_ext_>gap+epsilon & by_gse_ext_>gap+epsilon & bz_gse_ext_>gap+epsilon);
else
   'ERROR: gap flag cannot be zero'
end
bx_gse_ext=bx_gse_ext_(cond_good);
by_gse_ext=by_gse_ext_(cond_good);
bz_gse_ext=bz_gse_ext_(cond_good);
date_ext=date_ext_(cond_good);
B=B__(cond_good);
clear cond_good
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove data gaps and produce for bx_gse inside the cloud
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon=1e-4;
if gap>0,
   cond_good=find(bx_gse_ext_<gap-epsilon & by_gse_ext_<gap-epsilon & ...
   bz_gse_ext_<gap-epsilon & date_ext_>=initial_date & date_ext_<=end_date);
elseif gap<0,
   cond_good=find(bx_gse_ext_>gap+epsilon & by_gse_ext_>gap+epsilon & ...
   bz_gse_ext_>gap+epsilon & date_ext_>=initial_date & date_ext_<=end_date);
else
   'ERROR: gap flag cannot be zero'
end
% Comment: for simplicity in notation define bi as b_gse_i (i=x,y,z)
bx=bx_gse_ext_(cond_good);
by=by_gse_ext_(cond_good);
bz=bz_gse_ext_(cond_good);
date_=date_ext_(cond_good);
B_=B__(cond_good);
clear cond cond_good
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MV matrix (M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma^2_n=M*n*n
M(1,1)=mean(bx.^2)-mean(bx)^2;
M(2,2)=mean(by.^2)-mean(by)^2;
M(3,3)=mean(bz.^2)-mean(bz)^2;
M(1,2)=mean(bx.*by)-mean(bx)*mean(by);
M(1,3)=mean(bx.*bz)-mean(bx)*mean(bz);
M(2,3)=mean(by.*bz)-mean(by)*mean(bz);
M(2,1)=M(1,2);
M(3,1)=M(1,3);
M(3,2)=M(2,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigenvalues and eigenvectors
[V,D] =eig(M);
% Eigenvalues nor ordered at this stage.
lambda(1)=D(1,1);
lambda(2)=D(2,2);
lambda(3)=D(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of three unit vectors to get orientation of the structure:
n(1,1)=V(1,1); n(1,2)=V(2,1); n(1,3)=V(3,1);
n(2,1)=V(1,2); n(2,2)=V(2,2); n(2,3)=V(3,2);
n(3,1)=V(1,3); n(3,2)=V(2,3); n(3,3)=V(3,3);
clear V
%% This convention satisfies: M*n(1,:)=lambda(1)*n(1,:)
%% (idem to 2 y 3).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here the eigenvalues and unit vector are ordered
[lambda_ordenado,indices_lambda]=sort(lambda);
lambda_ordenado;
n_min=n(indices_lambda(1),:);
n_int=n(indices_lambda(2),:);
n_max=n(indices_lambda(3),:);
eval1=lambda_ordenado(1);
eval2=lambda_ordenado(2);
eval3=lambda_ordenado(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of x_versor_nube, y_versor_nube, z_versor_nube
%%%%% x_versor_nube=(+/-)n_min [+/-r_nube in the case p=0].
%%%%% y_versor_nube=(+/-)n_max [+/-r_phi_nube in the case p=0].
%%%%% z_versor_nube=(+/-)n_int [z_nube, main cylinder axis].
%%%%% Guess def due to possible changes of signs:
%%%%% Bz(0)>0, x_cloud=r(out_bound), right hand
x_versor_nube=n_min;
y_versor_nube=n_max;
z_versor_nube=n_int;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation of B (transformation to the cloud coordinates).
% ejes x_veror_nube, y_veror_nube, z_veror_nube.
bx_nube=bx*x_versor_nube(1)+by*x_versor_nube(2)+bz*x_versor_nube(3);
by_nube=bx*y_versor_nube(1)+by*y_versor_nube(2)+bz*y_versor_nube(3);
bz_nube=bx*z_versor_nube(1)+by*z_versor_nube(2)+bz*z_versor_nube(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check z_versor_nube sign: Bz(0)>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=fix(size(bz_nube,2)/2);
if bz_nube(NN-1)<0 & bz_nube(NN)<0 & bz_nube(NN+1)<0
   z_versor_nube=-n_int; % z_nube
   bx_nube=bx*x_versor_nube(1)+by*x_versor_nube(2)+bz*x_versor_nube(3);
   by_nube=bx*y_versor_nube(1)+by*y_versor_nube(2)+bz*y_versor_nube(3);
   bz_nube=bx*z_versor_nube(1)+by*z_versor_nube(2)+bz*z_versor_nube(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check x_versor_nube sign: x_cloud=r(out_bound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The convention is x_cloud is in agreement
% with r_versor for the out-bound branch,
% This is checked computung alfa (angle between x_versor_nube and x_GSE)
% |alfa| needs to be lower than 90 deg.
% |alfa| = arc cos (x_versor_nube  x_GSE)
% If |alfa| < 90 is ok, if not sign if definition of x_versor_nube need to be changed
mod_alfa=acos(x_versor_nube(1))*180/pi;
if mod_alfa < 85
 'Def. of x_versor_nube is Ok.';
elseif mod_alfa > 95
 'Def. of x_versor_nube is OPOSITE. INVERSE SIGN';
  x_versor_nube=-n_min; % r
  bx_nube=bx*x_versor_nube(1)+by*x_versor_nube(2)+bz*x_versor_nube(3);
  by_nube=bx*y_versor_nube(1)+by*y_versor_nube(2)+bz*y_versor_nube(3);
  bz_nube=bx*z_versor_nube(1)+by*z_versor_nube(2)+bz*z_versor_nube(3);
  mod_alfa=acos(x_versor_nube(1))*180/pi
  if mod_alfa >= 85
     'ERROR: CAMBIAR sgn x_versor no arregla condicion x_cloud=r(out_bound)'
     stop
  end
else
 'ERROR alfa GIVES near 90 grados !'
 'THE S/C PASS TROUGHT THE CLOUD AXIS (THE LEG)'
 'PATHOLOGIC CASE !'
 stop
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check y_versor_nube sign, such that the 3 new base be right-hand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if {x_versor_nube, y_versor_nube, z_versor_nube} is right-handed
MM_provi(1,:)=x_versor_nube;
MM_provi(2,:)=y_versor_nube;
MM_provi(3,:)=z_versor_nube;
MM_provi_det=det(MM_provi);
if abs(MM_provi_det-1) < 1e-10, % MM_provi_det == +1
 'Def. de y_versor_nube Ok.';
elseif abs(MM_provi_det+1) < 1e-10, % MM_provi_det == -1
 'Def. de y_versor_nube AL REVES. INVERTIR SIGNO.';
  y_versor_nube=-n_max; % phi
  bx_nube=bx*x_versor_nube(1)+by*x_versor_nube(2)+bz*x_versor_nube(3);
  by_nube=bx*y_versor_nube(1)+by*y_versor_nube(2)+bz*y_versor_nube(3);
  bz_nube=bx*z_versor_nube(1)+by*z_versor_nube(2)+bz*z_versor_nube(3);
  MM_provi(1,:)=x_versor_nube;
  MM_provi(2,:)=y_versor_nube;
  MM_provi(3,:)=z_versor_nube;
  MM_provi_det=det(MM_provi)
  if abs(MM_provi_det-1) > 1e-10,
     'ERROR: NO SE SOLUCIONA AL CAMBIO A TERNA DERECHA'
     stop
  end
else
 'ERROR, EL DET NO DA NI 1 NI -1 !!!!'
 stop
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the rotation matrix from range of B with exact boundaries
% here the extended rotated B is generated
bx_nube_ext=bx_gse_ext*x_versor_nube(1)+by_gse_ext*x_versor_nube(2)+ ...
bz_gse_ext*x_versor_nube(3);
by_nube_ext=bx_gse_ext*y_versor_nube(1)+by_gse_ext*y_versor_nube(2)+...
bz_gse_ext*y_versor_nube(3);
bz_nube_ext=bx_gse_ext*z_versor_nube(1)+by_gse_ext*z_versor_nube(2)+...
bz_gse_ext*z_versor_nube(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evec1=x_versor_nube;
evec2=z_versor_nube;
evec3=y_versor_nube;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if graf is diff from zero, plot B_cloud extended components.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graf ~= 0,
   cond_graf=find(date_ext>initial_date_graf & date_ext<end_date_graf);
   bmax=max([max(bx_nube_ext(cond_graf).*B(cond_graf)),...
   max(by_nube_ext(cond_graf).*B(cond_graf)),max(bz_nube_ext(cond_graf).*B(cond_graf))]);
   bmin=min([min(bx_nube_ext(cond_graf).*B(cond_graf)),...
   min(by_nube_ext(cond_graf).*B(cond_graf)),min(bz_nube_ext(cond_graf).*B(cond_graf))]);

   
   % B
   real_mod_B=sqrt((bx_nube_ext(cond_graf).*B(cond_graf)).^2+...
   (by_nube_ext(cond_graf).*B(cond_graf)).^2+...
   (bz_nube_ext(cond_graf).*B(cond_graf)).^2);


   % by and Fy
   clear fechaenh funciondelflujo fechaenhy funciondey  maximo indicemax fechaenhmax valordelafecha

   Fy=compute_Fy(by_nube,date_,Vsw);
   Fy_ext=compute_Fy(by_nube_ext,date_ext,Vsw);
   funciondelflujo=Fy/(max(abs(Fy)))*mean(B(cond_graf));
   funciondelflujo_ext=Fy_ext/(max(abs(Fy)))*mean(B(cond_graf));
   funciondey= by_nube_ext(cond_graf).*B(cond_graf);
   fechaenh=(date_-date_ext(cond_graf(1)))*24;
   fechaenhy=(date_ext(cond_graf)-date_ext(cond_graf(1)))*24;
 
  maximo=max(abs(funciondelflujo));
  maximosinmodulo=max(funciondelflujo);
  minimo=min(abs(funciondelflujo));
  minimosinmodulo=min(funciondelflujo);
  fin=funciondelflujo(end);
  indicefin=find(funciondelflujo==fin);
  indicemax=find(abs(funciondelflujo)==maximo);
  indicemaxsinmodulo=find(funciondelflujo==maximosinmodulo);
  indicemin=find(abs(funciondelflujo)==minimo);
  indiceminsinmodulo=find(funciondelflujo==minimosinmodulo);
  fechaenhmax=(fechaenh(indicemax))
  fechaenhmin=(fechaenh(indicemin))
  fechaenhmaxsinmodulo=(fechaenh(indicemaxsinmodulo))
  fechaenhminsinmodulo=(fechaenh(indiceminsinmodulo))
  fechaenfin=fechaenh(indicefin)
  indicedelfinenhy=find(fechaenfin==fechaenhy);
 
% clear fechaebhmax
 %fechaenhmax=27.5
 epsilon=0.01; 
%fechaenfrente=21.29
%fechaenback=28.99
% valordelafechaenfrente=datestr(date_ext(cond_graf(1))+(fechaenfrente/24))
 %valordelafechaenback=datestr(date_ext(cond_graf(1))+(fechaenback/24))
   
   fechaenback=(datenum(valordelafechaenback)-date_ext(cond_graf(1)))*24
   fechaenfrente=(datenum(valordelafechaenfrente)-date_ext(cond_graf(1)))*24
   indicefechaenback=min(find(fechaenhy>=fechaenback-epsilon&fechaenhy<=fechaenback+epsilon))
   indicefechaenfrente=min(find(fechaenhy>=fechaenfrente-epsilon&fechaenhy<=fechaenfrente+epsilon))
   
   
   
  valordelafecha=datestr(date_ext(cond_graf(1))+(fechaenhmax/24))
  funciondelflujo_ext1=funciondelflujo_ext(cond_graf);
  

    
    if funciondelflujo(end)>=funciondelflujo_ext1(indicedelfinenhy)
        delta=0;
         delta=funciondelflujo(end)-funciondelflujo_ext1(indicedelfinenhy);
        funciondelflujo_ext2=funciondelflujo_ext1+delta;
       
    end
    if funciondelflujo(end)<funciondelflujo_ext1(indicedelfinenhy)
         delta=0;
         delta=funciondelflujo_ext1(indicedelfinenhy)-funciondelflujo(end);
        funciondelflujo_ext2=funciondelflujo_ext1-delta;
    end
    
%   if funciondelflujo(indicemaxsinmodulo)==funciondelflujo(indicemax)
%      cociente=funciondelflujo(end)/funciondelflujo(indicemaxsinmodulo)
%      cocientefrente=funciondelflujo_ext2(indicefechaenfrente)/funciondelflujo(indicemaxsinmodulo)
%      cocienteback=funciondelflujo_ext2(indicefechaenback)/funciondelflujo(indicemaxsinmodulo)
%    end
%   if funciondelflujo(indiceminsinmodulo)==funciondelflujo(indicemax)
%      cociente=funciondelflujo(end)/funciondelflujo(indiceminsinmodulo)  
%      cocientefrente=funciondelflujo_ext2(indicefechaenfrente)/funciondelflujo(indiceminsinmodulo)
%      cocienteback=funciondelflujo_ext2(indicefechaenback)/funciondelflujo(indiceminsinmodulo)
%    end

    funcionenfrente=funciondelflujo_ext2(indicefechaenfrente);
    funcionenback=funciondelflujo_ext2(indicefechaenback);
    funciondelflujo_extnueva=funciondelflujo_ext2-funcionenfrente;
    funciondelflujo_nueva=funciondelflujo-funcionenfrente;

  if funciondelflujo(indicemaxsinmodulo)==funciondelflujo(indicemax)
     cociente=funciondelflujo_nueva(end)/funciondelflujo_nueva(indicemaxsinmodulo)
     cocientefrente=funciondelflujo_extnueva(indicefechaenfrente)/funciondelflujo_nueva(indicemaxsinmodulo)
     cocienteback=funciondelflujo_extnueva(indicefechaenback)/funciondelflujo_nueva(indicemaxsinmodulo)
   end
   if funciondelflujo(indiceminsinmodulo)==funciondelflujo(indicemax)
      cociente=funciondelflujo_nueva(end)/funciondelflujo_nueva(indiceminsinmodulo)  
      cocientefrente=funciondelflujo_extnueva(indicefechaenfrente)/funciondelflujo_nueva(indiceminsinmodulo)
      cocienteback=funciondelflujo_extnueva(indicefechaenback)/funciondelflujo_nueva(indiceminsinmodulo)
    end
 
fechaenhy_=fechaenhy-(initial_date-date_ext(cond_graf(1))).*24;
fechaenh_=fechaenhy-(initial_date-date_ext(cond_graf(1))).*24;
clear  fechaenhy  fechaenh
fechaenhy= fechaenhy_;
fechaenh= fechaenh_;
fechaenhmaxsinmodulo_=fechaenhmaxsinmodulo-(initial_date-date_ext(cond_graf(1))).*24;
fechaenhminsinmodulo_=fechaenhminsinmodulo-(initial_date-date_ext(cond_graf(1))).*24;
clear fechaenhmaxsinmodulo fechaenhminsinmodulo
fechaenhmaxsinmodulo=fechaenhmaxsinmodulo_;
fechaenhminsinmodulo=fechaenhminsinmodulo_;
fechaend=(end_date-initial_date).*24;
funciondez= bz_nube_ext(cond_graf).*B(cond_graf);
funciondex= bx_nube_ext(cond_graf).*B(cond_graf);

 clear comando;
if flag_sat==1
 'PRESS ANY KEY TO CONTINUE'
comando=strcat('save',DIR_WORK,'campoBy_H1_',num2str(jj),' fechaenhy initial_date end_date fechaenh bmin bmax fechaenhmaxsinmodulo fechaenhminsinmodulo fechaend funciondey funciondelflujo_extnueva fechaenhy   funciondelflujo indicemaxsinmodulo indiceminsinmodulo funciondez funciondex');
eval(comando); clear comando;
end
 clear comando;
if flag_sat==2
 'PRESS ANY KEY TO CONTINUE'
comando=strcat('save',DIR_WORK,'campoBy_H2_',num2str(jj),' fechaenhy initial_date end_date fechaenh bmin bmax fechaenhmaxsinmodulo fechaenhminsinmodulo fechaend funciondey funciondelflujo_extnueva fechaenhy   funciondelflujo indicemaxsinmodulo indiceminsinmodulo funciondez funciondex');
eval(comando); clear comando;
end
if flag_sat==3
comando=strcat('save',DIR_WORK,'campoBy_U_',num2str(jj),' fechaenhy initial_date end_date fechaenh bmin bmax fechaenhmaxsinmodulo fechaenhminsinmodulo fechaend funciondey funciondelflujo_extnueva fechaenhy   funciondelflujo indicemaxsinmodulo indiceminsinmodulo funciondez funciondex');
eval(comando); clear comando;
 figure
%By
 
   plot(fechaenhy,funciondey,'.','MarkerSize',graf)

   % plot(fechaenhminsinmodulo,funciondelflujo(indiceminsinmodulo),'d','MarkerFaceColor','r','MarkerSize',10)
   hold on

 plot(fechaenhy,funciondelflujo_extnueva,'.k','MarkerSize',graf)
 hold on
 hold on
 
  plot(fechaenhmaxsinmodulo,funciondelflujo(indicemaxsinmodulo),'d','MarkerFaceColor','r','MarkerSize',10)
   hold on
% plot(fechaenfrente,funciondelflujo_extnueva(indicefechaenfrente),'h','MarkerFaceColor','k','MarkerSize',10)
 hold on
% plot(fechaenback,funciondelflujo_extnueva(indicefechaenback),'h','MarkerFaceColor','r','MarkerSize',10)
 hold on
  %  axis([0 (date_ext(cond_graf(end))-date_ext(cond_graf(1)))*24 bmin bmax])
  % line([(initial_date-date_ext(cond_graf(1)))*24,(initial_date-date_ext(cond_graf(1)))*24],...
  % [bmin,bmax],'LineStyle','--','Color','k')
 line([0,0],[bmin,bmax],'LineStyle','--','Color','k')
 line([fechaend,fechaend],[bmin,bmax],'LineStyle','--','Color','k')

 %axis([-20 40 bmin bmax])
%%%%%%  line([fechaenfrente,fechaenfrente],[bmin,bmax],'LineStyle','--','Color','r')
%%%%%%%%   line([fechaenback,fechaenback],[bmin,bmax],'LineStyle','--','Color','r')
 
  % line([(end_date-date_ext(cond_graf(1)))*24,(end_date-date_ext(cond_graf(1)))*24],...
%   [bmin,bmax],'LineStyle','--','Color','k')
 %  line([0,(date_ext(cond_graf(end))-date_ext(cond_graf(1)))*24],[0,0],'LineStyle',':','Color','k')
 line([fechaenhy(1),fechaenhy(end)],[0,0],'LineStyle','-','Color','k')

   ylabel('B_{y,MV} (nT)')
  % xlabel(strcat('Time (Hours), after:',datestr(date_ext(cond_graf(1)))))
     xlabel(strcat('Time (Hours), after:',datestr(initial_date)))

%  print('-depsc',strcat(DIR_WORK(2:end),'figure1b'))

%pause


  
   close all
    clear fechaenh funciondelflujo fechaenhy funciondey  maximo indicemax fechaenhmax valordelafecha fechaenhmin fechaenhminsinmodulo
    clear fechaenhmaxsinmodulo fechaenfin indicemin indicemaxsinmodulo indiceminsinmodulo indicefin
    

  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if kind_mv==1,
   main_versor=z_versor_nube;
elseif kind_mv==2,
   main_versor=x_versor_nube;
else
   'kind_mv needs to be 1 or 2'
   stop
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of theta, phi, and gamma from the computed rotation matrix
% theta: latitud, on ecliptic theta=0, range:[-90,+90]
% phi: angle between projection of z_n on ecpliptic, measured from x_gse (phi=0) to +y_gse (phi=90)
% and beyond, range:[0,+360].
% Fixme, improve comment on gamma
% gamma: angulo que rota el plano perp al eje de la nube.
% Notar que NI theta NI phi dependen de como se elijen x_n e y_n
% (es decir que an tengo libertad para rotar el plano perp al eje del tubo)
% Esta libertad es usada (MV la usa) al elegir gamma.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta:
% the scalar product between z_nube and z_gse is: n_z_nube(3)
% and thus z_versor_nube(3)=cos(theta_monio)
% theta_monio_min_var is inside [0,180]
theta_monio_min_var=acos(main_versor(3))*180/pi;
% theta_monio: angle between z_gse and z_cloud, so:
theta_min_var=90-theta_monio_min_var;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod=sqrt(main_versor(1)^2+main_versor(2)^2);
p_zn_eclip_versor_x_gse=main_versor(1)/mod;
p_zn_eclip_versor_y_gse=main_versor(2)/mod;
cos_de_phi_min_var=p_zn_eclip_versor_x_gse;
sin_de_phi_min_var=p_zn_eclip_versor_y_gse;
'\phi (entre "proy. del eje nube en plano ecliptica" , "eje x_GSE"';
if sin_de_phi_min_var>0 & cos_de_phi_min_var>0,
 cuadrante=1;
 'Primer Cuadrante (0,90)';
 phi_min_var=acos(cos_de_phi_min_var)*180/pi;
elseif sin_de_phi_min_var>0 & cos_de_phi_min_var<0,
 cuadrante=2;
 'Segundo Cuadrante (90,180)';
 phi_min_var=acos(cos_de_phi_min_var)*180/pi;
elseif sin_de_phi_min_var<0 & cos_de_phi_min_var<0,
 cuadrante=3;
 'Tercer Cuadrante (180,270)';
 phi_min_var=360-acos(cos_de_phi_min_var)*180/pi;
elseif sin_de_phi_min_var<0 & cos_de_phi_min_var>0,
 cuadrante=4;
 'Cuarto Cuadrante (270,360)';
 phi_min_var=360-acos(cos_de_phi_min_var)*180/pi;
else
  'Error'
  stop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radio del tubo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance.
Delta_t=abs((initial_date-end_date))*24*3600;
% Delta_t is the time that the crafts pass trought the cloud.
Dist=Vsw*Delta_t/AU;
% Dist: Distance walked by the spacecraft inside cloud (AU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RADIO del tubo de flujo.
R_min_var=sqrt((sin(phi_min_var*pi/180)*Dist/2)^2 + ...
(cos(phi_min_var*pi/180)*sin(theta_min_var*pi/180)*Dist/2)^2);
% sqrt(x_n^2+y_n^2)
'Radio del tubo (suponiendo p=0):';
%Este calculo es equivalente a:
xn=(Dist/2)*x_versor_nube(1);
yn=(Dist/2)*y_versor_nube(1);
R_MV_verif=sqrt(xn^2+yn^2)
		gamma=calculo_gamma((phi_min_var*pi/180),(theta_min_var*pi/180));
		t=(date_-date_(1))*24*3600;
		[rho]=calculo_rho(Vsw/AU,t,gamma,(theta_min_var*pi/180),(phi_min_var*pi/180),0);
R_otro=(abs(rho(1))+abs(rho(end)))/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta=theta_min_var;
phi=phi_min_var;
R=R_min_var;
 

end