function [dummy]=parametros()
% Definicion de parametros.

% ------------------------------------------------------------
%file_swe='wi_h1s_swe_19970108000145_19970113235734.cdf';
file_swe='omni_hros_1min_20090205000000_20090220000000.cdf';
% carefull, keep this blanck space at the beginning!
%dir_=' /home/agulisano/codigoscdf/codigos_para_adri_dic7_2007/';
dir_=' \Users\Gulisano\Desktop\gemma\nube_13feb2009\ ';
spacecraft='Wind';
% ------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_mat_plasma='file_mat_plasma';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save parametros dir_ file_swe file_mat_plasma spacecraft
