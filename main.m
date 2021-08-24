function [Bx_nube,By_nube,Bz_nube,date,theta,phi,R,Bx_nube_ext,By_nube_ext,Bz_nube_ext,date_ext,eval1,eval2,eval3,evec1,evec2,evec3,cociente,cocientefrente,cocienteback]=main(DIR_WORK,dir_,Vsw,initial_date,end_date,graf,gap,normalized,kind_mv,initial_date_graf,end_date_graf,file_mfi_,valordelafechaenfrente,valordelafechaenback)
%
% Main sequence to get orientation of IP structures (as MCs/CS/Shocks) from MV method.
% v2.1, by S. Dasso, March 2008.
% Instituto de Astronom?a y F?sica del Espacio (IAFE, CONICET-UBA), Argentina.
% Departamento de F?sica, FCEN (UBA), Argentina.
%
% Input
%
% Output
%

comando=strcat('load',dir_,file_mfi_);
eval(comando)
clear comando

if normalized==0,
  B=Bx_gse-Bx_gse+1; % array of ones
elseif normalized==1 & kind_mv==1,
  % Normalization of vector series to catch rotation of vector
  B=sqrt(Bx_gse.^2+By_gse.^2+Bz_gse.^2);
else
  'Normalized needs to be 0 or 1, valid only for MCs'
  stop
end
bx_gse=Bx_gse./B;
by_gse=By_gse./B;
bz_gse=Bz_gse./B;

[bx_nube,by_nube,bz_nube,date,theta,phi,R,bx_nube_ext,by_nube_ext,bz_nube_ext,date_ext,eval1,eval2,eval3,evec1,evec2,evec3,B_,B,cociente,cocientefrente,cocienteback]=min_var_oct09(DIR_WORK,bx_gse,by_gse,bz_gse,fecha_,Vsw,initial_date,end_date,graf,gap,B,kind_mv,initial_date_graf,end_date_graf,valordelafechaenfrente,valordelafechaenback,jj);

Bx_nube_ext=bx_nube_ext.*B;
By_nube_ext=by_nube_ext.*B;
Bz_nube_ext=bz_nube_ext.*B;
Bx_nube=bx_nube.*B_;
By_nube=by_nube.*B_;
Bz_nube=bz_nube.*B_;
clear cond

