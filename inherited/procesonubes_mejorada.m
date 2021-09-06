%DIR_WORK=' C:\nubes_Luciano\datos\ '
DIR_WORK=' C:\Users\Gulisano\Desktop\gemma\nube_13feb2009\ '
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AU_=1.49597870691e8; % 1AU = 149,597,870,691 ± 30 metros from Wikipedia
for jj=1:1
%comando=strcat('load ',DIR_WORK,'omni',num2str(jj));
comando=strcat('load ',DIR_WORK,'file_mat_plasmaWind');
eval(comando)
clear comando

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir=DIR_WORK;
kind_mv=1;
% following is only valid if is MC
normalized=1; % 0:no-normalized, 1:normalized
% flag value for data gap
gap=1e10;

graf=6; % graf=x (diff from zero) means plot and size of points (graf=0 not plot)

switch jj
    
    case 1
       fecha_omni=fecha;   
    clear y m d h mi s
      fecha_inicio(jj)=datenum(2009,02,14,00,0,0);
      %fecha_inicio(jj)=datenum(2009,02,14,00,0,0);
   [y,m,d,h,mi,s]=datevec(fecha_inicio(jj));
 initial_date_mv=datenum(y,m,d,h,mi,s);
   clear y m d h mi s     
        fecha_fin(jj)=datenum(2009,02,14,16,0,0);
        %fecha_fin(jj)=datenum(2009,02,14,16,0,0);
        [y,m,d,h,mi,s]=datevec(fecha_fin(jj));
     end_date_mv=datenum(y,m,d,h,mi,s); 
      initial_date=fecha_inicio(jj);
      end_date=fecha_fin(jj); 
      clear y m d h mi s    
      [y,m,d,h,mi,s]=datevec(fecha_fin(jj));    
   valordelafechaenfrente=datestr(initial_date);
   valordelafechaenback=datestr(end_date_mv);
    clear y m d h mi s
   date_ini_new=initial_date;
   date_end_new=end_date;
 clear y m d h mi s
   [y,m,d,h,mi,s]=datevec(initial_date-1);
   initial_date_graf=datenum(y,m,d,0,0,0);
   clear y m d h mi s
   [y,m,d,h,mi,s]=datevec(end_date+2);
   end_date_graf=datenum(y,m,d);
   clear y m d h mi s
   quality=1
  overtak=1
   
end


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 index_omni=find(fecha_omni>=initial_date & fecha_omni<=end_date);
   index_graf=find(fecha_omni>=initial_date_graf & fecha_omni<=end_date_graf);
   
   Vsw=mean(Vx_gse(index_omni));
 file_mfi='datos_mag_nube_provisorio';
 BX_GSE=Bx_gse;
 BY_GSE=By_gse;
 BZ_GSE=Bz_gse;
 x=X;
clear Bx_gse By_gse Bz_gse X

Vr=(Vx_gse(index_graf));
Vsw=mean(Vx_gse(index_omni));
x=x(index_graf);
fecha_omni_=fecha_omni(index_graf);
 
 
  B_= sqrt(BX_GSE.^2+ BY_GSE.^2+ BZ_GSE.^2);
  indicedebuenos=(find( abs(B_) < 9999));
  B=B_(indicedebuenos);
 Bx_gse=BX_GSE(indicedebuenos);
 By_gse=BY_GSE(indicedebuenos);
 Bz_gse=BZ_GSE(indicedebuenos);
 fecha_=fecha(indicedebuenos);
 
  cond_graf=find(fecha_>=initial_date_graf & fecha_<=end_date_graf);
   comando=strcat('save',dir,file_mfi,' fecha_ Bx_gse By_gse Bz_gse jj');
   eval(comando), clear comando
%comando=strcat('load ',DIR_WORK,'WF',num2str(jj));
%eval(comando)
%clear comando
buenos=(find(abs(Vr)<9999));

date__=fecha_omni_(buenos);
%date__=fecha_omni(buenos);
%d=abs(x); %in AUs 
mejores=(find(abs(x)<9999));
d=abs(x(mejores).*(637800)/AU_); %to convert Re in AUs 

nMC=jj;

Vr=abs(Vr(buenos));
flag_graf=1;
%distance_to_sun=abs(X);
%distance_to_sun=abs(X.*6378)/AU;
distance_to_sun=abs(x(mejores).*637800)./(AU_);

 %[MCsize,travel_time_,Vfit_center,coef_l,dV_fit,mean_distance_to_sun,Pears_datos_]=compute_expansion(Vr,d,date__,initial_date ,end_date ,AU_,flag_graf,nMC,distance_to_sun,DIR_WORK,date_ini_new,date_end_new)
[delta_t,deltat,MCsize,MCsizepartelineal,sigma,travel_time_,travel_timepartelineal,Vfit_center,coef_l,coef_lpartelineal,dV_fit_nube,dV_fit_partelineal,mean_distance_to_sun,Pears_datos_,Vc_centropartelineal,a,b,vantes,vdespues]=compute_Vlinearpart_mejorada(Vr,d,date__,initial_date ,end_date, AU_,1,1,jj,distance_to_sun,DIR_WORK,date_ini_new,date_end_new,quality,overtak)
 
%[delta_t,deltat,MCsize,MCsizepartelineal,sigma,travel_time_,travel_timepartelineal,Vfit_center,coef_l,coef_lpartelineal,dV_fit_nube,dV_fit_partelineal,mean_distance_to_sun,Pears_datos_,Vc_centropartelineal,a,b,vantes,vdespues]=compute_Vlinearpart(Vr,d,date__,initial_date ,end_date, AU_,1,1,jj,distance_to_sun,DIR_WORK,date_ini_new,date_end_new,quality,overtak)
  
 %[Bx_nube,By_nube,Bz_nube,date,theta,phi,R,Bx_nube_ext,By_nube_ext,Bz_nube_ext,date_ext,eval1,eval2,eval3,evec1,evec2,evec3]=main(DIR_WORK,dir,Vsw,initial_date_mv,end_date_mv,graf,gap,normalized,kind_mv,initial_date_graf,end_date_graf,file_mfi)

  %[Bx_nube,By_nube,Bz_nube,date,theta,phi,R,Bx_nube_ext,By_nube_ext,Bz_nube_ext,date_ext,eval1,eval2,eval3,evec1,evec2,evec3]=main(DIR_WORK,dir,Vsw,initial_date_mv,end_date_mv,graf,gap,normalized,kind_mv,initial_date_graf,end_date_graf,file_mfi)
   [Bx_nube,By_nube,Bz_nube,date,theta,phi,R,Bx_nube_ext,By_nube_ext,Bz_nube_ext,date_ext,eval1,eval2,eval3,evec1,evec2,evec3,cociente,cocientefrente,cocienteback]=main(DIR_WORK,dir,Vsw,initial_date_mv,end_date_mv,graf,gap,normalized,kind_mv,initial_date_graf,end_date_graf,file_mfi,valordelafechaenfrente,valordelafechaenback)

%comando=strcat('save nube_',num2str(jj));
  % eval(comando), clear comando
 comando=strcat('save ',DIR_WORK, 'nube_expansion_',num2str(jj));
 eval(comando), clear comando

 
delta_t
deltat
mean_distance_to_sun
pause
clear Bx_nube By_nube Bz_nube date theta phi R Bx_nube_ext By_nube_ext Bz_nube_ext date_ext eval1 eval2 eval3 evec1 evec2 evec3 cociente cocientefrent cocienteback DIR_ Vs winitial_date_mv end_date_mv graf initial_date_graf end_date_graf file_mfi valordelafechaenfrente valordelafechaenback 

clear MCsize travel_time_ Vfit_center coef_l dV_fit mean_distance_to_sun Pears_datos_ Vr d date__ initial_date end_date nMC distance_to_sun date_ini_new date_end_new

clear delta_t deltat MCsize MCsizepartelineal sigma travel_time_ travel_timepartelineal Vfit_center coef_l coef_lpartelineal dV_fit_nube  dV_fit_partelineal mean_distance_to_sun Pears_datos_ Vc_centropartelineal a b vantes vdespues
close all

end  



%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




