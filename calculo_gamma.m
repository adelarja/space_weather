function [gamma]=calculo_gamma(phi_rad,theta_rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to compute gamma function [gamma]=calculo_gamma(phi_rad,theta_rad)
%input
%phi_rad,theta_rad
%output
%gamma

% Criterio para elegir gama: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%que x_n coincida con la direccion de la trayectoria.
epsilon=1e-15;
gama_rad=real(atan(-tan(phi_rad)/sin(theta_rad+epsilon)));
xn_dot_xGSE=cos(gama_rad)*sin(theta_rad)*cos(phi_rad)-sin(gama_rad)*sin(phi_rad);
if xn_dot_xGSE<0,
 gama_rad=real(pi+atan(-tan(phi_rad)/sin(theta_rad+epsilon)));
 xn_dot_xGSE=cos(gama_rad)*sin(theta_rad)*cos(phi_rad)-sin(gama_rad)*sin(phi_rad);
 if xn_dot_xGSE<0,
  ' ERROR CON EL CUADRANTE DE GAMMA !!'
 end
end
gamma=gama_rad;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


