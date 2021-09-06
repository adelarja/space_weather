function[delta_t,deltat,MCsize,MCsizepartelineal,sigma,travel_time_,travel_timepartelineal,Vfit_center,coef_l,coef_lpartelineal,dV_fit_nube,dV_fit_partelineal,mean_distance_to_sun,Pears_datos_,Vc_centropartelineal,a,b,vantes,vdespues]=compute_Vlinearpart_mejorada(Vr,d,date_,date_ini,date_end,AU,flag_graf,flag_sat,nMC,distance_to_sun,DIR_WORK,date_ini_new,date_end_new,quality,overtak)

% Function to compute expansion coef from radial component of V
% [travel_time,Vfit_center,coef_l,dV_fit]=compute_expansion(Vr,d,date_,date_ini,date_end,AU,flag_graf)
%
% Input
% Vr: radial velocity [array]
% d:  distance to the Sun [array], in AU
% date_: date [array], absolute matlab format
% date_ini: initial date for the cloud, absolute matlab format
% date_end: end     date for the cloud, absolute matlab format
% AU:       astronomical unit (in km)
% flag_graf: 1 for plot / 0 for no plot
% flag_sat: 1 for Helios1 / 2 for Helios2
% nMC: number of MC analized
% distance_to_sun: distance to the Sun in AUs [array]
% DIR_WORK:
% date_ini_new:
% date_end_new:
%
% Output
% MCsize:        (t_end-t_start)*Vfit_center (in AU)
% travel_time_: transit time, assuming constant velocity (Vfit_center), in days
% Vfit_center: V fitted evaluated at middle time (center)
% coef_l:      m*travel_time/Vfit_center (m is fitted slope)
% dV_fit:      Vfit_in-Vfit_out
% mean_distance_to_sun: mean value of distance to the Sun in AUs

% cloud itself
 index=find(date_>=date_ini & date_<=date_end);
 datestr(date_ini)
 index_cloud=find(date_>=date_ini_new & date_<=date_end_new);%parte lineal
 mean_distance_to_sun=mean(distance_to_sun(index));

% change of time reference and units for time
 v=Vr(index_cloud);% parte lineal
D_=d(index_cloud)*AU; % D_ in km complete 


time=(date_(index)-date_(index(1)))*24*3600; % time in seconds from date_ini
delta_t=time(end);%duration


tiempo=(date_(index_cloud)-date_(index_cloud(1)))*24*3600; % time in seconds from date_ini_new
deltat=tiempo(end);% duration of linear part
t_center=tiempo(end)/2;  % time at middle of linear part (proxy of linear center) desde date_ini_new
t_=(date_(index))*24*3600;%original
tiemp=(date_(index_cloud))*24*3600;
t=t_-tiemp(1);% tiempo medido a partir de date_ini_new
tiempoenelcentro=t_center-(1);%center if linear part
%tenelcentromenoserror=t_center-errorcentro_enhoras*3600-tiemp(1);
%tenelcentromaserror=t_center+errorcentro_enhoras*3600-tiemp(1);
t_centro=t(1)+delta_t/2;% time at middle of cloud (proxy of center) desde date_ini_new


clear y m d h mi s
[y,m,d,h,mi,s]=datevec(date_ini);
antes=datenum(y,m,d,h-delta_t/3600,mi,s);
clear y m d h mi s
[y,m,d,h,mi,s]=datevec(date_end);
despues=datenum(y,m,d,h+delta_t/3600,mi,s);
clear y m d h mi s


indexantes=find(date_>=antes & date_<=date_ini);
indexdespues=find(date_>=date_end & date_<=despues);
vantes=mean(Vr(indexantes));
vdespues=mean(Vr(indexdespues));

%pause
% Fitting
[Pol,S]=polyfit(tiempo,v,1);

Vout_linearpart=Pol(1)*tiempo(end)+Pol(2);% linear part
Vin_linearpart=Pol(1)*tiempo(1)+Pol(2);% linear part
Vout=Pol(1)*t(end)+Pol(2);% cloud
Vin=Pol(1)*t(1)+Pol(2);% cloud
%Voutorigen=Pol(1)*(t(end)+tiemp(1))+Pol(2);
%Vinorigen=Pol(1)*(t(1)+tiemp(1))+Pol(2);
V_ajuste=Pol(1)*t+Pol(2);% cloud
%V_ajusteorigen=Pol(1)*(t+tiemp(1))+Pol(2);

[Vfit,delta]=polyval(Pol,tiempo,S);
ene=length(v);
norma=((v-Vfit));




sigma=sqrt(sum((v-Vfit).^2)./(ene-2));

deltaa=sqrt(sum(tiempo.^2))*sigma/sqrt(ene.*sum(tiempo.^2)-(sum(tiempo))^2);% acording to numerical recipes pag 773 and sigs.
deltab=deltaa.*sqrt(sum(tiempo.^2)/ene);
a=(ene*sum(tiempo.*v)-sum(tiempo)*sum(v))/(ene*sum(tiempo.^2)-(sum(tiempo))^2);%v=a*tiempo+b
b=(sum(v)-a*sum(tiempo))/ene;


Vfit_center=polyval(Pol,t_centro,S);% en el medio de la nube
Vc_centropartelineal=polyval(Pol,tiempoenelcentro,S); %en el medio de la parte lineal
%Vc_centrofrommvmenoserror=polyval(Pol,tenelcentromenoserror,S);
%Vc_centrofrommvmaserror=polyval(Pol,tenelcentromaserror,S);
%deltaVc=abs(Vc_centrofrommvmaserror-Vc_centrofrommvmenoserror)/2;
Pears_datos=corrcoef(Vfit,v);
Pears_datos_=Pears_datos(1,2);
  





figure(10188)
    hold on
  plot((tiempo-t(1))/3600,norma./Vc_centropartelineal,'k*'); 
    axis tight
    xlabel('Time (in hours) from beginning of MC')
    ylabel('(V-Vfit)/Vc lineal part')
if flag_sat==1
    title(strcat('Cloud ',{' '},datestr(date_ini,0)))
     print('-depsc',strcat(DIR_WORK(2:end),'Vrlineal_',num2str(nMC)))
    % print('-dtiff',strcat(DIR_WORK(2:end),'tiffVr_',num2str(nMC)))

    

    end
  
% Compute coefficient
travel_time=mean(D_)/Vfit_center;% en segundos
travel_timepartelineal=(mean(D_)/Vc_centropartelineal)/(3600*24);%en dias
%traveltimeerrorendias=sqrt((mean(D_)^2/Vc_centrofrommv^4)*deltaVc^2)/(3600*24);
coef_l=Pol(1)*travel_time/Vfit_center;
coef_lpartelineal=a*mean(D_)/Vc_centropartelineal^2;
%coef_lerror=abs(sqrt((mean(D_)^2/Vc_centrofrommv^4)*deltaa^2+(2*mean(D_)*a/Vc_centrofrommv^3)*deltaVc^2));
%errorterminodelapendiente=(mean(D_)^2/Vc_centrofrommv^4)*deltaa^2;
%errorterminodelavelocidadc=(2*mean(D_)*a/Vc_centrofrommv^3)*deltaVc^2;
dV_fit_nube=Vin-Vout;
dV_fit_partelineal=Vin_linearpart-Vout_linearpart;
MCsize=delta_t*Vfit_center/AU;
MCsizepartelineal=deltat*Vc_centropartelineal/AU;
%MCsizeerror=abs(sqrt((delta_t^2/AU^2)*deltaVc^2));
travel_time_=travel_time/(3600*24);
travel_timepartelineal_=travel_timepartelineal/(3600*24);
%save provi

   initial_date_graf=date_ini-(delta_t/(3600*12));
   end_date_graf=date_end+(delta_t/(3600*12));
   index_graf=find(date_>=initial_date_graf & date_<=end_date_graf);


    figure(10101)
    hold on
   % plot((date_-date_(index(1)))*24,Vr,'*k','MarkerSize',1)
 plot((date_(index_graf)-date_(index(1)))*24,Vr(index_graf),'*k','MarkerSize',1)
    plot((tiempo-t(1))/3600,v,'g');
    plot((tiempo-t(1))/3600,Vfit,'r')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % plot(t_center/3600,Vfit_center,'sr','MarkerSize',10)
    plot((t_centro-t(1))/3600,Vfit_center,'sr','MarkerSize',10) 
     % plot((tiempoenelcentro-t(1))/3600,Vc_centropartelineal,'dm','MarkerFaceColor','m','MarkerSize',10) 
   % %  plot((tenelcentromenoserror-t(1))/3600,Vc_centrofrommvmenoserror,'dm','MarkerFaceColor','m','MarkerSize',10) 
  % %  plot((tenelcentromaserror-t(1))/3600,Vc_centrofrommvmaserror,'dm','MarkerFaceColor','m','MarkerSize',10) 
     plot((t-t(1))/3600,V_ajuste,'b')
     plot((t(end)-t(1))/3600,Vout,'ob','MarkerSize',10) 
     plot((t(1)-t(1))/3600,Vin,'ob','MarkerSize',10)
      plot((tiempo(1)-t(1))/3600,Vin_linearpart,'rh','MarkerSize',10)
       plot((tiempo(end)-t(1))/3600,Vout_linearpart,'rh','MarkerSize',10)
     axis tight
    xlabel('Time (in hours) from beginning of MC')
    ylabel('V (km/sec)')
    where=axis;
    x1=(where(2)-where(1))*(1/40)+where(1);
    y1=where(4)-(where(4)-where(3))/20;
    y2=where(4)-2*(where(4)-where(3))/20;
    y3=where(3)+abs(where(4)-where(3))/100;
    y4=y3+4*abs(where(4)-where(3))/100;
    y5=vantes+vantes.*0.2;


 %   text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center,3),' MCsize=',{' '},num2str(MCsize,3),' T=',{' '},num2str(travel_time_,3),' Dv=',{' '},num2str(dV_fit,3)))
   %text(x1,y1,strcat('coef l=',{' '},num2str(coef_l),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Dv=',{' '},num2str(dV_fit)))

 %  text(x1,y2,strcat('coef correl=',{' '},num2str(Pears_datos(1,2)),' distsun=',{' '},num2str(mean_distance_to_sun),{' '},'q=',num2str(quality),{' '},'ot=',num2str(overtak),'slope in km/sec= ',num2str(Pol(1))))
    %text(x1,y1,strcat('coef l=',{' '},num2str(coef_l),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Dv=',{' '},num2str(dV_fit)))
%     text(x1,y1,strcat('coef_lmv=',{' '},num2str(coef_lfrommv),'coef_lerr=',{' '},num2str(coef_lerror),' Vcmv=',{' '},num2str(Vc_centrofrommv),' Vcerr=',{' '},num2str(deltaVc),' distsun=',{' '},num2str(mean_distance_to_sun),{' '},'q=',num2str(quality),{' '},'ot=',num2str(overtak)))

% text(x1,y2,strcat(' MCsizemv=',{' '},num2str(MCsizefrommv),' MCsizeerr=',{' '},num2str(MCsizeerror),'slope in km/sec= ',num2str(Pol(1)),' T=',{' '},num2str(travel_timefrommv),' Terr=',{' '},num2str(traveltimeerrorendias)))
 
    
switch nMC
    case 1
    text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 2
    text(x1,665,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))

    case 3
  text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))

     case 4
    text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
   
    case 5
      text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),'Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))

    case 6

text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
   
    case 7
    text(x1,400,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
    case 8
    text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
  case 9
    text(x1,530,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
  case 10
    text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
    case 11
    text(x1,510,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
    case 12
    text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Nopert %=',{' '},num2str((deltat/delta_t)*100,3)))
 
    case 13
     %text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
     text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str(100)))
case 14
      text(x1,510,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' M3size=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
    case 15
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 16
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 17
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 18
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 19
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 20
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 21
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 22
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 23
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
case 24
         text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' pert %=',{' '},num2str((deltat/delta_t)*100,3)))
end
  text(x1,y2,strcat('coef correl=',{' '},num2str(Pears_datos(1,2)),' distsun=',{' '},num2str(mean_distance_to_sun),{' '},'q=',num2str(quality),{' '},'ot=',num2str(overtak),'slope in km/sec= ',num2str(Pol(1))))
 %text(x1,y1,strcat('coef l=',{' '},num2str(coef_l,3),' Vc=',{' '},num2str(Vfit_center),' MCsize=',{' '},num2str(MCsize),' T=',{' '},num2str(travel_time_),' Dv=',{' '},num2str(dV_fit_partelineal,3)))

 

 %text(x1,y3,strcat('coef_lerror=',{' '},num2str(coef_lerror),' Vcerror=',{' '},num2str(deltaVc)))
 %text(x1,y4,strcat(' MCsizeerror=',{' '},num2str(MCsizeerror),' Terror=',{' '},num2str(traveltimeerrorendias)))

    line([(date_ini-date_(index(1)))*24,(date_ini-date_(index(1)))*24],...
    [min(Vr),max(Vr)],'LineStyle','--','Color','k')
    line([(date_end-date_(index(1)))*24,(date_end-date_(index(1)))*24],...
    [min(Vr),max(Vr)],'LineStyle','--','Color','k')


    line([(date_(indexantes(1))-date_(index(1)))*24,(date_(indexantes(end))-date_(index(1)))*24],...
    [vantes,vantes],'LineStyle','-','Color','r','LineWidth',3)
    line([(date_(indexdespues(1))-date_(index(1)))*24,(date_(indexdespues(end))-date_(index(1)))*24],...
    [vdespues,vdespues],'LineStyle','-','Color','r','LineWidth',3)

        
 %text(x1,y5,strcat('vantes= ',num2str(vantes),'vdespues=',num2str(vdespues)))

 %   line([(date_ini_new-date_(index(1)))*24,(date_ini_new-date_(index(1)))*24],...
 %   [min(Vr),max(Vr)],'LineStyle','--','Color','g')
 %   line([(date_end_new-date_(index(1)))*24,(date_end_new-date_(index(1)))*24],...
 %   [min(Vr),max(Vr)],'LineStyle','--','Color','g')

%pause

   
    line([(date_(indexantes(1))-date_(index(1)))*24,(date_(indexantes(end))-date_(index(1)))*24],...
    [vantes,vantes],'LineStyle','-','Color','r','LineWidth',3)
    line([(date_(indexdespues(1))-date_(index(1)))*24,(date_(indexdespues(end))-date_(index(1)))*24],...
    [vdespues,vdespues],'LineStyle','-','Color','r','LineWidth',3)
   title(strcat('Cloud ',{' '},datestr(date_ini,0)))
     print('-depsc',strcat(DIR_WORK(2:end),'NVr_',num2str(nMC)))
    print('-dtiff',strcat(DIR_WORK(2:end),'tiffVr_',num2str(nMC)))


    pause

%save provi
