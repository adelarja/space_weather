function [dummy]=proceso_cdf_wind_swe_omni()
% Para procesar archivos .cdf que se bajan de pags de NASA.

load parametros
file_cdf=strcat(dir_(2:end),file_swe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OMNI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comando=strcat('info = cdfinfo(''',file_cdf,''')');
eval(comando)
%info = cdfinfo(file_cdf)
vars = info.Variables

% Fecha
clear comando
comando=strcat('data = cdfread(''',file_cdf,''',''','variable''',',','''','Epoch',''');')
eval(comando)

for ii=1:size(data,1),
	 fecha_(ii)=data{ii};
         fecha(ii)=todatenum(fecha_(ii));
end

% Componentes de Bgse
clear comando
comando=strcat('Bxgse = cdfread(''',file_cdf,''',''','variable''',',','''','BX_GSE',''');')
eval(comando)
clear comando
comando=strcat('Bygse = cdfread(''',file_cdf,''',''','variable''',',','''','BY_GSE',''');')
eval(comando)
clear comando
comando=strcat('Bzgse = cdfread(''',file_cdf,''',''','variable''',',','''','BZ_GSE',''');')
eval(comando)

for ii=1:size(Bxgse,1),
  clear Bgse_
  Bxgse_ = Bxgse{ii};
  Bygse_ = Bygse{ii};
  Bzgse_ = Bzgse{ii};
  Bx_gse(ii)=double(Bxgse_);
  By_gse(ii)=double(Bygse_);
  Bz_gse(ii)=double(Bzgse_);
end
% Componentes de Vgse
clear comando
comando=strcat('Vxgse = cdfread(''',file_cdf,''',''','variable''',',','''','Vx',''');')
eval(comando)
clear comando
comando=strcat('Vygse = cdfread(''',file_cdf,''',''','variable''',',','''','Vy',''');')
eval(comando)
clear comando
comando=strcat('Vzgse = cdfread(''',file_cdf,''',''','variable''',',','''','Vz',''');')
eval(comando)

for ii=1:size(Vxgse,1),
  clear Vgse_
  Vxgse_ = Vxgse{ii};
  Vygse_ = Vygse{ii};
  Vzgse_ = Vzgse{ii};
  Vx_gse(ii)=double(Vxgse_);
  Vy_gse(ii)=double(Vygse_);
  Vz_gse(ii)=double(Vzgse_);
end



%densidad de protones


clear comando
comando=strcat('Np_ = cdfread(''',file_cdf,''',''','variable''',',','''','proton_density',''');')
eval(comando)
for ii=1:size(Np_,1),
  Np = Np_{ii};
  np(ii)=double(Np);
end

% temperatura
clear comando
comando=strcat('temperature_ = cdfread(''',file_cdf,''',''','variable''',',','''','T',''');')
eval(comando)
for ii=1:size(temperature_,1),
  temperature = temperature_{ii};
  Temperature(ii)=double(temperature);
end

%Beta 
clear comando
comando=strcat('beta_ = cdfread(''',file_cdf,''',''','variable''',',','''','Beta',''');')
eval(comando)
for ii=1:size(beta_,1),
  beta = beta_{ii};
  Beta(ii)=double(beta);
end

% posicion Space Craft 
clear comando
comando=strcat('x_ = cdfread(''',file_cdf,''',''','variable''',',','''','x',''');')
eval(comando)
for ii=1:size(x_,1),
  x = x_{ii};
  X(ii)=double(x);
end

% AE_INDEX 
clear comando
comando=strcat('AE_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','AE_INDEX',''');')
eval(comando)
for ii=1:size(AE_INDEX_,1),
  AE_INDEX = AE_INDEX_{ii};
 AE(ii)=double(AE_INDEX);
end

% AL_INDEX 
clear comando
comando=strcat('AL_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','AL_INDEX',''');')
eval(comando)
for ii=1:size(AL_INDEX_,1),
  AL_INDEX = AL_INDEX_{ii};
 AL(ii)=double(AL_INDEX);
end

% AU_INDEX 
clear comando
comando=strcat('AU_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','AU_INDEX',''');')
eval(comando)
for ii=1:size(AU_INDEX_,1),
  AU_INDEX = AU_INDEX_{ii};
 AU(ii)=double(AU_INDEX);
end

% SYM_D_INDEX 
clear comando
comando=strcat('SYM_D_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','SYM_D',''');')
eval(comando)
for ii=1:size(SYM_D_INDEX_,1),
  SYM_D_INDEX = SYM_D_INDEX_{ii};
 SYM_D(ii)=double(SYM_D_INDEX);
end

% SYM_H_INDEX 
clear comando
comando=strcat('SYM_H_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','SYM_H',''');')
eval(comando)
for ii=1:size(SYM_H_INDEX_,1),
  SYM_H_INDEX = SYM_H_INDEX_{ii};
 SYM_H(ii)=double(SYM_H_INDEX);
end

%ASY_D_INDEX 
clear comando
comando=strcat('ASY_D_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','ASY_D',''');')
eval(comando)
for ii=1:size(ASY_D_INDEX_,1),
  ASY_D_INDEX = ASY_D_INDEX_{ii};
 ASY_D(ii)=double(ASY_D_INDEX);
end

%ASY_H_INDEX 
clear comando
comando=strcat('ASY_H_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','ASY_H',''');')
eval(comando)
for ii=1:size(ASY_H_INDEX_,1),
  ASY_H_INDEX = ASY_H_INDEX_{ii};
 ASY_H(ii)=double(ASY_H_INDEX);
end


%PC_N_INDEX
clear comando
comando=strcat('PC_N_INDEX_ = cdfread(''',file_cdf,''',''','variable''',',','''','PC_N_INDEX',''');')
eval(comando)
for ii=1:size(PC_N_INDEX_,1),
  PC_N_INDEX = PC_N_INDEX_{ii};
 PC_N(ii)=double(PC_N_INDEX);
end

% NRL, page 30: vTp=9.79*10^5*sqrt(Tp)   (cm/sec)
%eV=1.14e4;
%tp=(V_th/9.79).^2*eV; % tpr is in Kelvin

comando=strcat('save',dir_,file_mat_plasma,spacecraft,' fecha Vx_gse Vy_gse Vz_gse Bx_gse By_gse Bz_gse np Temperature Beta X AE AL AU SYM_D SYM_H ASY_D ASY_H PC_N');
eval(comando)
clear comando
close all, clear, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
