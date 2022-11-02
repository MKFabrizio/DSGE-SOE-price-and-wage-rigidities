clc; clear all;
addpath('C:\dynare\5.0\matlab')


if ~isdir('Latex/Efficient')
    mkdir('Latex/Efficient')
end

%==========================================================================
% Main Parameters
%==========================================================================
BBETA      = 0.99; % 0.995
RSTRSTSS   = 1/BBETA;
GGAMMAC    = 1; %(Above 1 - 1.4 with current paramterization explodes in NK) CO:1
GGAMMA     = 0.6;  %Between 0 and  3 homebias
CCHI       = 3; %1.5; 
GGAMMAD    = GGAMMAC;
EPS        = 6;%3.8; % ELasticity of substitution over varieties.
GGAMMAST   = GGAMMA;
MMBAR      = 1;
OOMEGA     = 500;
SSIGMA     = 1;
EPSH       = 2; % (1-3) CO:1
EPSF       = 2; % (1-3) CO:1
MMU        = 1.1;
TAUH       = (EPS*GGAMMA - EPS + 1)/(EPS*GGAMMA);
TAUH_W      = 0;
TAUH_P      = 0;
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
STD_C      = 0.0001;
STD_Z      = 0.0289;
B0=0;
B1=0;
B2=0;
B3=0;
PHI_PIH    = 2.78;
PHI_PIH_2    = 2.78;
PHI_Y = 0.5;
PHI_G   = 0;
FXIR=0;
EPSW= 6; %4.3; % ELasticity of substitution over varieties LABOR
PPSI=1/3;
TTHETAHW = 0.75;

	lambda_w = ((1-TTHETAHW)/TTHETAHW)*(1-BBETA*TTHETAHW)/(1+CCHI*EPSW);
	lambda = ((1-BBETA*TTHETAH)*(1-TTHETAH))/TTHETAH;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;

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

results_OSR_TAYLOR = zeros(7,8);
results_adhoc = zeros(7,8);
results_bank = zeros(7,8);

%% Results simulation

%A. Baseline Taylor CPI
c = 0
h = waitbar(0,'Initializing waitbar...');


%B. Optimal Taylor PPI. OSR in Taylor rule PPI
PHI_PIH    = 2.78;
%PHI_PIH_2    = 2.78;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0  -DRez=0 -DXgap=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0   -DRez=0 -DXgap=1
save OC_TAYLOR0 oo_ M_;
results_OC_TAYLOR(1:9,1)= 100*values_all;
results_OC_TAYLOR(10,1)= (values_all(7)^2 + 0.25*values_all(4)^2)*100;
results_OC_TAYLOR(11,1)= PHI_PIH;
results_OC_TAYLOR(12,1)= (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;

%CPI estricto
c=c+(100/7)
waitbar(c/100,h,sprintf('%d%% along...',c))

PHI_PIH    = 2.78;

save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TAUH TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2;
dynare nk_wr.mod -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DStrictCPI=1 -Dsubsidy=0 -DRez=0 -DXgap=1
pause(1);
save OC_TAYLOR1 oo_ M_;
results_OC_TAYLOR(1:9,2)= 100*values_all;
results_OC_TAYLOR(10,2)= (values_all(7)^2 + 0.25*values_all(4)^2)*100;
results_OC_TAYLOR(11,2)= PHI_PIH;
results_OC_TAYLOR(12,2)= (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;

c=c+(100/7)
waitbar(c/100,h,sprintf('%d%% along...',c))

%B. Output gap y Taylor
PHI_PIH    = 2.78;
%PHI_PIH_2    = 2.78;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0  -DRez=1 -DXgap=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0   -DRez=1 -DXgap=1
save OC_TAYLOR2 oo_ M_;
results_OC_TAYLOR(1:9,3)= 100*values_all;
results_OC_TAYLOR(10,3)= (values_all(7)^2 + 0.25*values_all(4)^2)*100;
results_OC_TAYLOR(11,3)= PHI_PIH;
results_OC_TAYLOR(12,3)= (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;

c=c+(100/7)
waitbar(c/100,h,sprintf('%d%% along...',c))

%C. Optimal Taylor CPI. OSR in CPI
PHI_PIH    = 1.06;

save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1  -DRez=0 -DXgap=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1   -DRez=0 -DXgap=1
save OC_TAYLOR3 oo_ M_;
results_OC_TAYLOR(1:9,4)= 100*values_all;
results_OC_TAYLOR(10,4)= (values_all(7)^2 + 0.25*values_all(4)^2)*100;
results_OC_TAYLOR(11,4)= PHI_PIH;
results_OC_TAYLOR(12,4)= (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;

c=c+(100/7)
waitbar(c/100,h,sprintf('%d%% along...',c))

%D. Optimal NGDP. OSR
PHI_PIH = 0;
PHI_G = 1.1;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1   -DRez=0 -DXgap=1
pause(1);
PHI_G = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1   -DRez=0 -DXgap=1
save OC_TAYLOR4 oo_ M_;
results_OC_TAYLOR(1:9,5)= 100*values_all;
results_OC_TAYLOR(10,5)= (values_all(7)^2 + 0.25*values_all(4)^2)*100;
results_OC_TAYLOR(11,5)= 0;
results_OC_TAYLOR(12,5)= (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;

c=c+(100/7)
waitbar(c/100,h,sprintf('%d%% along...',c))

%E. Optimal Wage inflation taylor. OSR
PHI_PIH    = 2.78;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -DtaylorW=1 -DXgap=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -DtaylorW=1  -DXgap=1
save OC_TAYLOR5 oo_ M_;
results_OC_TAYLOR(1:9,6)= 100*values_all;
results_OC_TAYLOR(10,6)= (values_all(7)^2 + 0.25*values_all(4)^2)*100;
results_OC_TAYLOR(11,6)= 0;
results_OC_TAYLOR(12,6)= (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;

c=c+(100/7)
waitbar(c/100,h,sprintf('%d%% along...',c))

%Optimal policy
PHI_Y = 0.16;

save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=0 -DXgap=1 -DOptimal=1
save OC_TAYLOR6 oo_ M_;
results_OC_TAYLOR(1:9,7)= 100*values_all;
results_OC_TAYLOR(10,7)= (values_all(7)^2 + 0.25*values_all(4)^2)*100;
results_OC_TAYLOR(11,7)= 0;
results_OC_TAYLOR(12,7)= (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;

c=c+(100/7)
waitbar(c/100,h,sprintf('%d%% along...',c))


%% 
%F. Total de trablas
Param = [0 0;0 0.75;0.75 0; 0.5 0.5; 0.25 0.75 ; 0.75 0.25 ; 0.5 0.75 ; 0.75 0.75 ]
c = 0
h = waitbar(0,'Initializing waitbar...');
for i= 1:8
   TTHETAH=Param(i,1) 
   TTHETAHW=Param(i,2)

    if (Param(i,1)==0 && Param(i,2)~=0)
        results_OSR_TAYLOR(1,i)= 0;
        c=c+(100/56)
    else
        
    %A. Optimal Taylor PPI. OSR in Taylor rule PPI
    PHI_PIH    = 2.78;
    %PHI_PIH_2    = 2.78;
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2 lambda lambda_w;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0  -DRez=0 -DXgap=1
    pause(1);
    PHI_PIH = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2 lambda lambda_w;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0   -DRez=0 -DXgap=1
    results_OSR_TAYLOR(1,i)= 100*values_all(8);
    results_adhoc(1,i) = (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;
    results_bank(1,i) =  (values_all(7)^2 + (1/16)*values_all(4)^2)*100;
    
    c=c+(100/56)
    waitbar(c/100,h,sprintf('%d%% along...',c))
    end
    
    %B. Strict CPI
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TAUH TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH PHI_PIH_2;
    dynare nk_wr.mod -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DStrictCPI=1 -Dsubsidy=0 -DRez=0 -DXgap=1
    pause(1);
    results_OSR_TAYLOR(2,i)= 100*values_all(8);
    results_adhoc(2,i) = (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;
    results_bank(2,i) =  (values_all(7)^2 + (1/16)*values_all(4)^2)*100;
    
        c=c+(100/56)
    waitbar(c/100,h,sprintf('%d%% along...',c))
    
    %  Classic Taylor
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0  -DRez=1 -DXgap=1
    pause(1);
    PHI_PIH = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0   -DRez=1 -DXgap=1
    results_OSR_TAYLOR(3,i)= 100*values_all(8);
    results_adhoc(3,i) = (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;
    results_bank(3,i) =  (values_all(7)^2 + (1/16)*values_all(4)^2)*100;
    
        c=c+(100/56)
    waitbar(c/100,h,sprintf('%d%% along...',c))
    
    %C. Optimal Taylor CPI. OSR in CPI
    PHI_PIH    = 1.06;

    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1 -DXgap=1
    pause(1);
    PHI_PIH = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1 -DXgap=1
    results_OSR_TAYLOR(4,i)= 100*values_all(8);
    results_adhoc(4,i) = (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;
    results_bank(4,i) =  (values_all(7)^2 + (1/16)*values_all(4)^2)*100;

    c=c+(100/56)
    waitbar(c/100,h,sprintf('%d%% along...',c))

    %D. Optimal NGDP. OSR
    PHI_PIH = 0;
    PHI_G = 1.1;
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1   -DRez=0  -DXgap=1
    pause(1);
    PHI_G = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1   -DRez=0 -DXgap=1
    results_OSR_TAYLOR(5,i)= 100*values_all(8);
    results_adhoc(5,i) = (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;
    results_bank(5,i) =  (values_all(7)^2 + (1/16)*values_all(4)^2)*100;

    c=c+(100/56)
    waitbar(c/100,h,sprintf('%d%% along...',c))

    %E. Optimal Wage inflation taylor. OSR
    PHI_PIH    = 2.78;
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -DtaylorW=1  -DXgap=1
    pause(1);
    PHI_PIH = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -DtaylorW=1   -DXgap=1
    results_OSR_TAYLOR(6,i)= 100*values_all(8);
    results_adhoc(6,i) = (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;
    results_bank(6,i) =  (values_all(7)^2 + (1/16)*values_all(4)^2)*100;

    c=c+(100/56)
    waitbar(c/100,h,sprintf('%d%% along...',c))

    %Optimal policy
    if (Param(i,1)==0 && Param(i,2)==0)
        results_OSR_TAYLOR(7,i)= 100*values_all(8);
        c=c+(100/56)
    elseif (Param(i,1)==0) || (Param(i,2)==0)
        results_OSR_TAYLOR(7,i)= 0;
        c=c+(100/56)
    else
        PHI_Y = 0.16;

        save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
        dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=0 -DXgap=1 -DOptimal=1
        results_OSR_TAYLOR(7,i)= 100*values_all(8);
        results_adhoc(7,i) = (GGAMMA/2)*((1+CCHI)*values_all(4)^2 + (EPS/lambda)*values_all(7)^2 + (TTHETAHW/lambda_w)*values_all(6)^2)*100;
        results_bank(7,i) =  (values_all(7)^2 + (1/16)*values_all(4)^2)*100;

        c=c+(100/56)
        waitbar(c/100,h,sprintf('%d%% along...',c))
    end
end



options_.noprint=0;
options_.val_precis=2;
headers_string={'Reglas';'\theta_p=\theta_w=0';'\theta_p=0,\theta_w=0.75';'\theta_p=0.75,\theta_w=0';'\theta_p=0.5,\theta_w=0.5';'\theta_p=0.25,\theta_w=0.75';'\theta_p=0.75,\theta_w=0.25';'\theta_p=0.5,\theta_w=0.75';'\theta_p=0.75,\theta_w=0.75'}; 
labels_string={'PPI';'Strict CPI';'Taylor';'CPI'; 'NGDP'; 'WAGE'; 'Optimal'};
labels_string_tex={'Baseline';'\pi_H';'\pi_CPI';'NGDP';'\pi_W';'Optimal'};
dyntable(options_,'OSRTAYLOR',headers_string,labels_string,results_OSR_TAYLOR,size(labels_string,4)+4,3,7)
headers_string={'\\ \hline ';'\theta_p = \theta_w = 0';'\theta_p=0, \theta_w=0.75';'\theta_p=0.75, \theta_w=0';'\theta_p=0.5, \theta_w=0.5';'\theta_p=0.25, \theta_w=0.75';'\theta_p=0.75, \theta_w=0.25';'\theta_p=0.5, \theta_w=0.75';'\theta_p=0.75, \theta_w=0.75  \hline'};
%headers_string = [headers_string{:}]
%second = append('\\gamma_c=', string(GGAMMAC), ', \\varepsilon_H = \\varepsilon_F = ', string(EPSH),', \\theta_H =', string(TTHETAH) ,', \\theta_W =', string(TTHETAHW));
dyn_latex_table_modified(options_,char(second),'Standard Dev. and Welfare Losses under Sticky Prices (Monetary Policy Rules)','Latex/OSRTAYLOR',headers_string,labels_string_tex,results_OSR_TAYLOR,size(labels_string,4)+4,3,6);
%dyn_latex_table_modified(options_,'GENTAYLOR','Latex/GENTAYLOR',headers_string,labels_string_tex,results_OC_TAYLOR,size(labels_string,2)+2,5,7); 

%%
%G. Strict CPI
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=1
save OC_TAYLOR6 oo_ M_;
results_OC_TAYLOR(1:9,7)= 100*values_all;
results_OC_TAYLOR(10,7)= 0;
results_OC_TAYLOR(11,7)= 0;
results_OC_TAYLOR(12,7)= 0;

%%
%H. Strict Growth
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=0 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=1
save OC_TAYLOR6 oo_ M_;
results_OC_TAYLOR(1:9,7)= 100*values_all;
results_OC_TAYLOR(10,7)= 0;
results_OC_TAYLOR(11,7)= 0;
results_OC_TAYLOR(12,7)= 0;

%%
%H. Output gap
PHI_Y = 2.06;
PHI_PIH = 2.7;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_Y;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=0 -DXgap=1
%pause(1);
%PHI_Y = x_opt_hat(1);
%PHI_PIH = x_opt_hat(2);
%save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA RHOA RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_Y;
%dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=0 -DXgap=1
save OC_TAYLOR5 oo_ M_;
results_OC_TAYLOR(1:9,6)= 100*values_all;
results_OC_TAYLOR(10,6)= 0;
results_OC_TAYLOR(11,6)= PHI_PIH;
results_OC_TAYLOR(12,6)= PHI_Y;


%% Tablas
%;'Strict PPI';'Strict CPI';'Strict NGDP'
% ; '\textit{Strict PPI}';'\textit{Strict CPI}';'\textit{Strict NGDP}

options_.noprint=0;
options_.val_precis=2;
headers_string={' ';'Opt-PPI';'Srict.-CPI';'Taylor';'Opt.-CPI';'Opt.-NGDP';'Wage Taylor';'Optimal Policy'}; 
labels_string={'sigma(Y)';'sigma(C)';'sigma(L)'; 'sigma(X gap)'; 'sigma(dep)'; 'sigma(piw)'; 'sigma(pih)'; 'W unc.'; 'W cond.';'b0';'phi_pi';'phi_g'};
labels_string_tex={'\sigma(Y)';'\sigma(C)';'\sigma(L)';'\sigma(X gap)';'\sigma(\Delta S)'; '\sigma(\pi_w)' ; '\sigma(\pi_H)'; 'W unc.'; 'W cond'; '\varphi_n' ;'\phi_{\pi^H}';'\phi_{\pi^H}';'\phi_G'};
dyntable(options_,'GENTAYLOR',headers_string,labels_string,results_OC_TAYLOR,size(labels_string,4)+4,3,7)
headers_string={'\\ \hline ';'\textit{Opt-PI}';'\textit{Strict.-CPI}';'\textit{Taylor}';'\textit{Opt.-CPI}'; '\textit{Opt.-NGDP}'; '\textit{Wage Taylor}' ; '\textit{Optimal Policy}  \hline'};
%headers_string = [headers_string{:}]
second = append('\\gamma_c=', string(GGAMMAC), ', \\varepsilon_H = \\varepsilon_F = ', string(EPSH),', \\theta_H =', string(TTHETAH) ,', \\theta_W =', string(TTHETAHW));
dyn_latex_table_modified(options_,char(second),'Standard Dev. and Welfare Losses under Sticky Prices (Monetary Policy Rules)','Latex/GENTAYLOR',headers_string,labels_string_tex,results_OC_TAYLOR,size(labels_string,4)+4,3,6);
%dyn_latex_table_modified(options_,'GENTAYLOR','Latex/GENTAYLOR',headers_string,labels_string_tex,results_OC_TAYLOR,size(labels_string,2)+2,5,7); 


%% Graphs eps_a

TAYLOR=load('OC_TAYLOR2.mat');
PPI=load('OC_TAYLOR3.mat');
CPI=load('OC_TAYLOR4.mat');
NGDP=load('OC_TAYLOR6.mat');
var_string={'log_yh','log_c','itit','log_pitpit','log_pihpih','log_dep','x_gap','nxnx'};

figure
for fig_iter=1:length(var_string)
    subplot(3,3,fig_iter)
    plot(1:15,TAYLOR.oo_.irfs.([var_string{fig_iter},'_eps_a']),'b-',1:15,PPI.oo_.irfs.([var_string{fig_iter},'_eps_a']),'g--',1:15,CPI.oo_.irfs.([var_string{fig_iter},'_eps_a']),'r-x',1:15,NGDP.oo_.irfs.([var_string{fig_iter},'_eps_a']),'c-s','Linewidth',1.5)
    grid on
    title(TAYLOR.M_.endo_names_long(strmatch(var_string{fig_iter},TAYLOR.M_.endo_names,'exact'),:))
end

legend('Regla de Taylor','CPI','NGDP','Óptimo','FontSize',12)

%% COPPER
TAYLOR=load('OC_TAYLOR2.mat');
PPI=load('OC_TAYLOR3.mat');
CPI=load('OC_TAYLOR4.mat');
NGDP=load('OC_TAYLOR6.mat');
var_string={'log_yh','log_c','itit','log_pitpit','log_pihpih','log_dep','log_w','cca'};

figure
for fig_iter=1:length(var_string)
    subplot(3,3,fig_iter)
    plot(1:15,TAYLOR.oo_.irfs.([var_string{fig_iter},'_eps_yc']),'b-',1:15,PPI.oo_.irfs.([var_string{fig_iter},'_eps_yc']),'g--',1:15,CPI.oo_.irfs.([var_string{fig_iter},'_eps_yc']),'r-x',1:15,NGDP.oo_.irfs.([var_string{fig_iter},'_eps_yc']),'c-s','Linewidth',1.5)
    grid on
    title(TAYLOR.M_.endo_names_long(strmatch(var_string{fig_iter},TAYLOR.M_.endo_names,'exact'),:))
end
legend('Regla de Taylor','CPI','NGDP','Óptimo','FontSize',12)

%% Graphs eps_z

TAYLOR=load('OC_TAYLOR2.mat');
PPI=load('OC_TAYLOR3.mat');
CPI=load('OC_TAYLOR4.mat');
NGDP=load('OC_TAYLOR6.mat');
var_string={'log_yh','log_c','itit','log_pitpit','log_pihpih','log_dep','log_w','log_q'};

figure
for fig_iter=1:length(var_string)
    subplot(3,3,fig_iter)
    plot(1:15,TAYLOR.oo_.irfs.([var_string{fig_iter},'_eps_z']),'b-',1:15,PPI.oo_.irfs.([var_string{fig_iter},'_eps_z']),'g--',1:15,CPI.oo_.irfs.([var_string{fig_iter},'_eps_z']),'r-x',1:15,NGDP.oo_.irfs.([var_string{fig_iter},'_eps_z']),'c-s','Linewidth',1.5)
    grid on
    title(TAYLOR.M_.endo_names_long(strmatch(var_string{fig_iter},TAYLOR.M_.endo_names,'exact'),:))
end
legend('Regla de Taylor','CPI','NGDP','Óptimo','FontSize',12)
%% Graphs eps_C

TAYLOR=load('OC_TAYLOR0.mat');
PPI=load('OC_TAYLOR3.mat');
CPI=load('OC_TAYLOR4.mat');
NGDP=load('OC_TAYLOR6.mat');
var_string={'log_yh','log_c','itit','log_pitpit','log_pihpih','log_dep','nxnx','log_q'};

figure
for fig_iter=1:length(var_string)
    subplot(3,3,fig_iter)
    plot(1:15,TAYLOR.oo_.irfs.([var_string{fig_iter},'_eps_C']),'b-',1:15,PPI.oo_.irfs.([var_string{fig_iter},'_eps_C']),'g--',1:15,CPI.oo_.irfs.([var_string{fig_iter},'_eps_C']),'r-x',1:15,NGDP.oo_.irfs.([var_string{fig_iter},'_eps_C']),'c-s','Linewidth',1.5)
    grid on
    title(TAYLOR.M_.endo_names_long(strmatch(var_string{fig_iter},TAYLOR.M_.endo_names,'exact'),:))
end
legend('Regla de Taylor','CPI','NGDP','Óptimo','FontSize',12)
%% Grafico 7 - welfare, wage rigidities & trade elasticity
GAMMA_s = 0.1:0.1:0.9
results_ = zeros(length(GAMMA_s),4);
for i = 1:1:length(GAMMA_s)

	TTHETAHW = GAMMA_s(i);
    
	PHI_G = 1.2;
    EPSF = 1;
    EPSH = 1;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1
pause(1);
PHI_G = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0  -Dngdp=1
results_(i,1)= 100*values_all(8);

    PHI_G = 1.2;
    EPSF = 1.5;
    EPSH = 1.5;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0  -Dngdp=1
pause(1);
PHI_G = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0  -Dngdp=1 
results_(i,2)= 100*values_all(8);

    PHI_G = 1.2;
    EPSF = 2;
    EPSH = 2;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1 
pause(1);
PHI_G = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1
results_(i,3)= 100*values_all(8);

    PHI_G = 1.2;
    EPSF = 2;
    EPSH = 2;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DtaylorW=1
pause(1);
PHI_G = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DtaylorW=1
results_(i,4)= 100*values_all(8);


end; 
%%
plot(GAMMA_s, results_) 
    xlabel('Rigidez de salarios');
    ylabel('Pérdida de bienestar');
    legend('\eta = 1','\eta = 2','\eta = 3')
    
%% chi
%GAMMA_s = 4:1:12;
%GAMMA_s = 0:0.1:0.9
GAMMA_s = 0.15:0.05:0.95

results2_ = zeros(length(GAMMA_s),4);

for i = 1:1:length(GAMMA_s)

	%EPSF = GAMMA_s(i);
    %EPSH = EPSF;
    
    %RHOA = GAMMA_s(i);
    %EPSW = GAMMA_s(i);
    % TTHETAHW = GAMMA_s(i);
    %CCHI =  GAMMA_s(i);
    GGAMMA =  GAMMA_s(i);
    GGAMMAST   = GGAMMA;
    %TAUH       = (EPS*GGAMMA - EPS + 1)/(EPS*GGAMMA);
    
	PHI_G = 1.2;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1  -DXgap=1
pause(1);
PHI_G = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0  -Dngdp=1 -DXgap=1
results2_(i,1)= 100*values_all(8);

    PHI_PIH = 2.78;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1  -Dngdp=0 -DXgap=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1  -Dngdp=0  -DXgap=1
results2_(i,2)= 100*values_all(8);

    PHI_PIH = 2.78;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0 -Dngdp=0 -DtaylorW=0 -DXgap=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=1 -DtaylorCPI=0 -Dngdp=0 -DtaylorW=0 -DXgap=1
results2_(i,3)= 100*values_all(8);

    PHI_PIH = 2.78;
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DtaylorW=1 -DXgap=1
pause(1);
PHI_PIH = x_opt_hat(1);
save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH PHI_G;
dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DtaylorW=1 -DXgap=1
results2_(i,4)= 100*values_all(8);


%save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
%dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=0 -DStrictPPI=0 -DStrictCPI=0 -DStrictG=0 -DXgap=1 -DOptimal=1
%results2_(i,5)= 100*values_all(8);

end; 

plot(GAMMA_s, results2_) 
    %xlabel('Elastividad de sustitución hogar-extranjero (\eta)')
    xlabel('Grado de Home Bias');
    %xlabel('Inversa de elasticidad de Frisch ($\varphi$)');
    %xlabel('Persistencia de shock')
    %xlabel('Elasticidad de sustitución de trabajo')
    ylabel('Pérdida de bienestar');
    legend('NGDP','CPI','PPI','Wage inflation')
    colormap cool;
%% Matrix NGDP

    TTHETAHWs = 0.0:0.05:0.95
    TTHETAHs  = 0.0:0.05:0.95
    c = 0
    h = waitbar(0,'Initializing waitbar...');
    %Función de pérdida
    Ls_ngdp= -Inf*ones(length(TTHETAHWs),length(TTHETAHs));
    for ii = 1:length(TTHETAHWs)
    for jj = 1:length(TTHETAHs)
        TTHETAHW = TTHETAHWs(ii);
        TTHETAH = TTHETAHs(jj);
        
     PHI_G = 1.2;
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1 -DXgap=1
    pause(1);
    PHI_G = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -Dngdp=1 -DXgap=1
    Ls_ngdp(ii,jj) = 100*values_all(8);
    
    c=c+(100/(400))
    waitbar(c/100,h,sprintf('%d%% along...',c))
    end;
    end;
% Matriz CPI

    TTHETAHWs = 0.0:0.05:0.95
    TTHETAHs  = 0.0:0.05:0.95
    c = 0
     h = waitbar(0,'Initializing waitbar...');
    Ls_cpi= -Inf*ones(length(TTHETAHWs),length(TTHETAHs));
    for ii = 1:length(TTHETAHWs)
    for jj = 1:length(TTHETAHs)
        TTHETAHW = TTHETAHWs(ii);
        TTHETAH = TTHETAHs(jj);
        
     PHI_PIH = 1.2;
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1 -Dngdp=0 -DXgap=1
    pause(1);
    PHI_PIH = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=1 -Dngdp=0 -DXgap=1
    Ls_cpi(ii,jj) = 100*values_all(8);
    
    c=c+0.25
     waitbar(c/100,h,sprintf('%d%% along...',c))
    end;
    end;


    %% Grafico 3D perdida de bienestar
  
    figure(1);
    %mesh(TTHETAHs,TTHETAHWs,Ls*10000);
    surf(TTHETAHs,TTHETAHWs,Ls_ngdp)
    xlabel('Rigidez precio');
    ylabel('Rigidez salario');
    zlabel('Pérdida de bienestar');
    FaceColor = 'flat';
    colormap cool;
    colorbar
    figure(2);
    contourf(TTHETAHs,TTHETAHWs,Ls_ngdp);
    xlabel('Rigidez precio');
    ylabel('Rigidez salario');
    colormap cool;

%%
X = zeros(size(Ls_ngdp));
Y = zeros(size(Ls_ngdp));
Z = zeros(size(Ls_ngdp));
T = zeros(size(Ls_ngdp));
for ii = 1:length(TTHETAHWs)
for jj = 1:length(TTHETAHs)
    x1 = round(Ls_ngdp(ii,jj),2);
    x2 = round(Ls_cpi(ii,jj),2);
    if x1==x2
        T(ii,jj) = 2;
    elseif Ls_ngdp(ii,jj) < Ls_cpi(ii,jj)
        T(ii,jj) = 1;
    %else
    %    T(ii,jj) = 0;
    end
    X(ii,jj) = (ii-1)/20;
    Y(ii,jj) = (jj-1)/20;
end
end
    

X=reshape(X,[],1)
Y=reshape(Y,[],1)
Z=reshape(T,[],1)

figure(1);
scatter(X,Y,60,Z,'filled');
ylabel('Rigidez precio (\theta_p)');
xlabel('Rigidez salario (\theta_w)');
figure(2); %está mal
pcolor(TTHETAHs,TTHETAHWs,T)
ylabel('Rigidez precio (\theta_p)');
xlabel('Rigidez salario (\theta_w)');

%% Matriz Wage Inflation

    TTHETAHWs = 0.0:0.05:0.95
    TTHETAHs  = 0.0:0.05:0.95
    c = 0
     h = waitbar(0,'Initializing waitbar...');
    Ls_piwpiw= -Inf*ones(length(TTHETAHWs),length(TTHETAHs));
    for ii = 16:length(TTHETAHWs)
    for jj = 1:length(TTHETAHs)
        TTHETAHW = TTHETAHWs(ii);
        TTHETAH = TTHETAHs(jj);
        
     PHI_PIH = 1.2;
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -DtaylorW=1
    pause(1);
    PHI_PIH = x_opt_hat(1);
    save PARAM_M1 BBETA GGAMMA CCHI EPS EPSH EPSF GGAMMAST MMBAR OOMEGA SSIGMA EPSW RHOA PPSI TTHETAHW RHONST STD_PSI STD_A STD_Z TAUH_P TAUH_W FXIR GGAMMAC TTHETAH PHI_PIH;
    dynare nk_wr -DCalvo=1 -Dlogutility=1 -Dosr=1 -Dramsey_policy=0 -Ddixit=0 -DtaylorPPI=0 -DtaylorCPI=0 -DtaylorW=1
    Ls_piwpiw(ii,jj) = 100*values_all(8);
    
    c=c+0.25
     waitbar(c/100,h,sprintf('%d%% along...',c))
    end;
    end;
    
    
%% Grafico 3D perdida de bienestar
  
    figure(1);
    %mesh(TTHETAHs,TTHETAHWs,Ls*10000);
    surf(TTHETAHs,TTHETAHWs,Ls_piwpiw)
    xlabel('Rigidez precio');
    ylabel('Rigidez salario');
    zlabel('Pérdida de bienestar');
    FaceColor = 'flat';
    colormap cool;
    colorbar
    figure(2);
    contourf(TTHETAHs,TTHETAHWs,Ls_piwpiw);
    xlabel('Rigidez precio');
    ylabel('Rigidez salario');
    colormap cool;
