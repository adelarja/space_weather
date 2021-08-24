function [Fy]=compute_Fy(by_nube,date_b,Vsw)
%
% function to compute acumulated Flux Fy, assuming constant velocity
%
% Input
% by_nube:
% date_b:
% Vsw:
%
% Output
% Fy: Acumulated flux Fy per unit length (as in Dasso et al., 2006, A&A) [nT]

t_hours=(date_b-date_b(1))*24;
Integrando=by_nube;
for gg=2:size(t_hours,2),
   signed_flux_phi_(gg)=Vsw*trapz(t_hours(1:gg),Integrando(1:gg));
end

Fy=signed_flux_phi_;






