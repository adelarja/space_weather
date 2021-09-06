function [rho]=calculo_rho(Vsw_,t,gama_rad,theta_rad,phi_rad,y0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo rho(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OJO AQUI HAY QUE DECIDIR SI SE DEFINE EPSILON POSITIVO 
% O NEGATIVO. PARA EVITAR QUE CAMBIE DE LADO Y APAREZCA
% SALTO EN rho(t). CON ESTE EVENTO CORRESPONDE negativo.
signo_epsilon=-1;
epsilon=signo_epsilon*(t(fix(end/2)+1)-t(fix(end/2)))/10;

%epsilon=0;
M_AU=-Vsw_*t(end)/2+Vsw_*t; % [M_AU]=AU
AAA=real(M_AU*(cos(gama_rad)*sin(theta_rad)*cos(phi_rad)-sin(gama_rad)*sin(phi_rad)));
BBB=real(y0+M_AU*(-sin(gama_rad)*sin(theta_rad)*cos(phi_rad)-cos(gama_rad)*sin(phi_rad)));
%rho=real(sqrt(AAA.^2+BBB.^2));

rho_old=(real(sqrt(AAA.^2+BBB.^2))).*sign(epsilon+t-t(fix(end/2)));

rho_nuevo=(real(sqrt(AAA.^2+BBB.^2)));
indexmin=min(find(rho_nuevo==(min(rho_nuevo))));
rho=[rho_nuevo(1:indexmin).*(-1),rho_nuevo(indexmin+1:end)];
%rho=rho_old;
%plot (abs(rho),'bd')
%pause
%plot(rho,'*')
%stop


%y0
%rho=real(sqrt(AAA.^2+BBB.^2).*sign(epsilon+t-t(fix(end/2))));
if (abs(rho)>=y0),
 %'ok'
else
 'error hay al menos un valor de rho que es menor que y0'
 clc
 rho
 y0
 stop
end

%plot(rho,'*')
%stop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



