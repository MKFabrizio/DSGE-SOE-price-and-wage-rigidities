clc; clear all;
addpath('C:\dynare\5.0\matlab')


if ~isdir('Latex/Efficient')
    mkdir('Latex/Efficient')
end

%==========================================================================
% Main Parameters
%==========================================================================
COEFS      = 12;

BBETA      = 0.995; % 0.995
RSTRSTSS   = 1/BBETA;
GGAMMAC    = 1; %(Above 1 - 1.4 with current paramterization explodes in NK)
GGAMMA     = 0.85;  %Between 0 and  3
CCHI       = 1.5; %1.5;
GGAMMAD    = GGAMMAC;
EPS        = 1.5; % ELasticity of substitution over varieties.
GGAMMAST   = GGAMMA;
MMBAR      = 1;
OOMEGA     = 500;
SSIGMA     = 1;
EPSH       = 1; % (1-3)
EPSF       = 1; % (1-3)
MMU        = 1.1;
TAUH       = (EPS*GGAMMA - EPS + 1)/(EPS*GGAMMA);
TTHETAH     = 0.75;
TTHETAF     = 0;
RRHOI      = 0;
RHOA       = 0.74;
RHOZ       = 0.6;
RHOBCS     = 0;
RHONST     = 0;
OPT        = 0;
FXIR       = 0;
STD_PSI    = 0.01; %SWITCH
STD_A      = 0.0064; 
STD_Z      = 0.0289;
B0=0;
B1=0;
B2=0;
B3=0;
PHI_PIH    = 2.78;
PHI_G   = 0;
TTHETAH= 0.5;
FXIR=0;
EPSW= 2; % ELasticity of substitution over varieties LABOR
PPSI=1/3;
TTHETAHW = 0.75;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA  SSIGMA RHOA RHONST STD_PSI STD_A TAUH FXIR GGAMMAC TTHETAH; 

Ramsey      = 0;  %do not do Ramsey, not implemented (Needs to write trend depreciation model)
osr         = 1;  %do all optimal rules
logutility  = 1;
dixit       = 0;
OSR_OC_FLEX = zeros(4,1);
OSR_CES_FLEX = zeros(4,1); 
OSR_DS_FLEX = zeros(4,1);
OSR_DSCES_FLEX = zeros(4,1);
OSR_OC_STICKY = zeros(4,1);
OSR_CES_STICKY = zeros(4,1); 
OSR_DS_STICKY = zeros(4,1);
OSR_DSCES_STICKY = zeros(4,1);
results_opt_rule_F = zeros(3,1);
results_opt_rule_S = zeros(3,1);
results_OC_TAYLOR = zeros(12,1);

%% Flexible Prices

%A. Baseline Taylor CPI

PHI_PIH    = 2.78;

save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr.mod -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0
pause(1);
save OC_TAYLOR0 oo_ M_;
results_OC_TAYLOR(1:9,1)= 100*values_all;
results_OC_TAYLOR(10,1)= -FXIR;
results_OC_TAYLOR(11,1)= PHI_PIH;
results_OC_TAYLOR(12,1)= 0;
%%
%B. Optimal Taylor PPI. OSR in Taylor rule PPI
PHI_PIH    = 2.78;

save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0 
save OC_TAYLOR1 oo_ M_;
results_OC_TAYLOR(1:9,2)= 100*values_all;
results_OC_TAYLOR(10,2)= -FXIR;
results_OC_TAYLOR(11,2)= PHI_PIH;
results_OC_TAYLOR(12,2)= 0;
%%
%C. Optimal Taylor CPI. OSR in CPI
PHI_PIH    = 2.78;

save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1 
save OC_TAYLOR2 oo_ M_;
results_OC_TAYLOR(1:9,3)= 100*values_all;
results_OC_TAYLOR(10,3)= -FXIR;
results_OC_TAYLOR(11,3)= PHI_PIH;
results_OC_TAYLOR(12,3)= 0;
%%
%D. Optimal Taylor CPI. OSR in CPI
PHI_PIH = 0;
PHI_G = 1.2;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1 
pause(1);
PHI_G = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1 
save OC_TAYLOR3 oo_ M_;
results_OC_TAYLOR(1:9,4)= 100*values_all;
results_OC_TAYLOR(10,4)= 0;
results_OC_TAYLOR(11,4)= 0;
results_OC_TAYLOR(12,4)= PHI_G;

%% 
%E. Strict PPI
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=1
save OC_TAYLOR4 oo_ M_;
results_OC_TAYLOR(1:9,5)= 100*values_all;
results_OC_TAYLOR(10,5)= 0;
results_OC_TAYLOR(11,5)= 0;
results_OC_TAYLOR(12,5)= 0;

%%
%F. Strict CPI
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=1
save OC_TAYLOR5 oo_ M_;
results_OC_TAYLOR(1:9,6)= 100*values_all;
results_OC_TAYLOR(10,6)= 0;
results_OC_TAYLOR(11,6)= 0;
results_OC_TAYLOR(12,6)= 0;

%%
%G. Strict Growth
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=1
save OC_TAYLOR6 oo_ M_;
results_OC_TAYLOR(1:9,7)= 100*values_all;
results_OC_TAYLOR(10,7)= 0;
results_OC_TAYLOR(11,7)= 0;
results_OC_TAYLOR(12,7)= 0;

%%
%H. Strict Output gap
PHI_Y = 0.16;
PHI_PIH = 43277;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_Y;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=0 -DXgap=1
pause(1);
PHI_Y = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_Y;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=0 -DXgap=1
save OC_TAYLOR7 oo_ M_;
results_OC_TAYLOR(1:9,8)= 100*values_all;
results_OC_TAYLOR(10,8)= 0;
results_OC_TAYLOR(11,8)= 0;
results_OC_TAYLOR(12,8)= PHI_Y;
%% Tablas
options_.noprint=0;
options_.val_precis=2;
headers_string={' ';'Base-PPI';'Opt.-PPI';'Opt.-CPI';'Opt.-NGDP';'Strict PPI';'Strict CPI';'Strict NGDP';'Output gap'};
labels_string={'sigma(Y)';'sigma(C)';'sigma(L)'; 'sigma(wedge)'; 'sigma(dep)'; 'sigma(pi)'; 'sigma(pih)'; 'W unc.'; 'W cond.';'b0';'phi_pi';'phi_g'};
labels_string_tex={'\sigma(Y)';'\sigma(C)';'\sigma(L)';'\sigma(\lambda)';'\sigma(\Delta S)'; '\sigma(\pi)' ; '\sigma(\pi_H)'; 'W unc.'; 'W cond'; '\varphi_n' ;'\phi_{\pi^H}';'\phi_{\pi^H}';'\phi_G'};
dyntable(options_,'GENTAYLOR',headers_string,labels_string,results_OC_TAYLOR,size(labels_string,2)+2,5,4)
headers_string={'\\ \hline ';'\textit{Base-PPI}';'\textit{Opt.-PPI}';'\textit{Opt.-CPI}'; '\textit{Opt.-NGDP}'; '\textit{Strict PPI}';'\textit{Strict CPI}';'\textit{Strict NGDP}'; '\textit{Output gap} \hline'};
headers_string=[headers_string{:}]
second = append('\\gamma_c=', string(GGAMMAC), ', \\varepsilon_H = \\varepsilon_F = ', string(EPSH),', \\theta_H =', string(TTHETAH));
dyn_latex_table_modified(options_,char(second),'Standard Dev. and Welfare Losses under Sticky Prices (Monetary Policy Rules)','Latex/GENTAYLOR',headers_string,labels_string_tex,results_OC_TAYLOR,size(labels_string,2)+2,8,6); 


%% Graphs eps_a

TAYLOR=load('OC_TAYLOR0.mat');
PPI=load('OC_TAYLOR1.mat');
CPI=load('OC_TAYLOR2.mat');
NGDP=load('OC_TAYLOR3.mat');
var_string={'log_yh','log_c','itit','log_pitpit','atat','log_dep','ldld','log_q'};

figure
for fig_iter=1:length(var_string)
    subplot(4,2,fig_iter)
    plot(1:15,TAYLOR.oo_.irfs.([var_string{fig_iter},'_eps_a']),'b-',1:15,PPI.oo_.irfs.([var_string{fig_iter},'_eps_a']),'g--',1:15,CPI.oo_.irfs.([var_string{fig_iter},'_eps_a']),'r-x',1:15,NGDP.oo_.irfs.([var_string{fig_iter},'_eps_a']),'c-s')
    grid on
    title(TAYLOR.M_.endo_names_long(strmatch(var_string{fig_iter},TAYLOR.M_.endo_names,'exact'),:))
end
%% Graphs eps_z

TAYLOR=load('OC_TAYLOR0.mat');
PPI=load('OC_TAYLOR1.mat');
CPI=load('OC_TAYLOR2.mat');
NGDP=load('OC_TAYLOR3.mat');
var_string={'log_yh','log_c','itit','log_pitpit','x_gap','log_dep','ltlt','log_q'};

figure
for fig_iter=1:length(var_string)
    subplot(4,2,fig_iter)
    plot(1:15,TAYLOR.oo_.irfs.([var_string{fig_iter},'_eps_z']),'b-',1:15,PPI.oo_.irfs.([var_string{fig_iter},'_eps_z']),'g--',1:15,CPI.oo_.irfs.([var_string{fig_iter},'_eps_z']),'r-x',1:15,NGDP.oo_.irfs.([var_string{fig_iter},'_eps_z']),'c-s')
    grid on
    title(TAYLOR.M_.endo_names_long(strmatch(var_string{fig_iter},TAYLOR.M_.endo_names,'exact'),:))
end
