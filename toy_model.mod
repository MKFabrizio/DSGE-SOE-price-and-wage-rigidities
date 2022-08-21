/*
 * This file implements the mod-file for the Welfare analysis in Ortiz/Montoro (2021):
 * "The Portfolio Channel of Capital Flows". 
 *
 * It computes IRFs and consumption equivalents relative to the natural allocation as described
 * in the paper. For doing so, it uses a solver to find the correct consumption equivalent lambda
 * based on the recursive lifetime welfare measures defined within the model block.
 
* THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER

* This implementation follows the code of Johannes Pfeifer. In case you spot mistakes,
* email: ma.ortizs@up.edu.de

* Key steps:
* 1: Derive the model 
* 2: Run model with bad rules, OSR and Ramsey exercises for welfare
* 3: Find consumption equivalent welfare loss with respect to natural equilibrium
*/

% Define the model used:

@#ifndef taylorPPI
    @#define taylorPPI = 0
@#endif

@#ifndef taylorCPI
    @#define taylorCPI = 0
@#endif

@#ifndef StrictPPI
    @#define StrictPPI = 0
@#endif

@#ifndef StrictCPI
    @#define StrictCPI = 0 
@#endif

@#ifndef logutility
    @#define logutility = 1
@#endif

@#ifndef ramsey_policy
    @#define ramsey_policy = 0
@#endif

@#ifndef Calvo
    @#define Calvo = 0
@#endif 

@#ifndef dixit
    @#define dixit = 0 
@#endif

@#ifndef ngdp
    @#define ngdp = 0 
@#endif

@#ifndef StrictG
    @#define StrictG = 0 
@#endif

@#ifndef Xgap
    @#define Xgap = 0 
@#endif

% @#ifndef ColeObstfeld
%    @#define ColeObstfeld = 0
% @#endif        

% @#ifndef efficient_steady_state
%    @#define efficient_steady_state = 0
% @#endif        

@#ifndef osr
    @#define osr=1
@#endif        

var qtqt        ${Q}$                       (long_name='Real Exchange Rate')
    totot       ${\tau}$                    (long_name='Terms of Trade (P^H/P^F)')
    rtrt        ${R}$                       (long_name='Real Interest Rate')
    ltlt        ${L}$                       (long_name='Labour')
    ctct        ${C}$                       (long_name='Consumption')
    wtwt        ${W}$                       (long_name='Real Wage')
    chch        ${C^{H}}$                   (long_name='Domestic Goods Consumption')
    cfcf        ${C^{F}}$                   (long_name='Domestic Foreign Consumption')
    xtxt        ${X}$                       (long_name='Exports')
    cstcst      ${C^{\ast}}$                (long_name='Foreign Consumption')
    rstrst      ${R^{\ast}}$                (long_name='Foreign Interest Rate')
    yhyh        ${Y^{H}}$                   (long_name='Output')
    ynyn        ${Y^{N}}$                   (long_name='Natural Output')
    x_gap       ${X}$                       (long_name='Output gap')
    atat        ${A}$                       (long_name='Technology Process')
    mcmc        ${MC}$                      (long_name='Marginal Cost')
    cca         ${CA}$                      (long_name='Current Account')
    dstdst      ${D^{\ast}}$                (long_name='Dealers Foreign Position')
    nxnx        ${NX}$                      (long_name='Net Exports')
    bcsbcs      ${B^{cb,\ast}}$             (long_name='Central Bank Foreign Position')
    nstnst      ${N^{\ast}}$                (long_name='Portfolio Flows')
    dep         ${\Delta S}$                (long_name='Nominal Exchange Rate Variation')
    wedge       ${\lambda}$                 (long_name='Backus-Smith Wedge')
    btbt        ${B}$                       (long_name='Household Domestic Position')
    util        ${U}$                       (long_name='Period Utility')
    welf        ${\mathbb{W}}$              (long_name='Recursive Welfare')
    welfgap     ${\mathbb{W}_{gap}}$        (long_name='Welfare gap')
    welfeq      ${\mathbb{W}^{nat}_{gap}}$  (long_name='Recursive Natural Welfare')
    log_yh      ${\log{Y}^H}$               (long_name='log output')
    log_yn      ${\log(Y)^N}$               (long_name='log natural output')
    log_w       ${\log{W/P}}$               (long_name='log real wage')
    log_l       ${\log{L}}$                 (long_name='log real labour')
    log_a       ${\log{A}}$                 (long_name='log technology level')
    log_c       ${\log{C}}$                 (long_name='log consumption level')
    log_wedge   ${\log{\lambda}}$           (long_name='log Backus-Smith wedge')
    log_q       ${\log{Q}}$                 (long_name='log real exchange rate')
    log_dep     ${\log{\Delta s}}$          (long_name='log nominal exchange rate variation')
    cnat        ${C^{CP}}$                  (long_name='Central planner consumption')
    itit       ${i}$                        (long_name='Nominal Interest rate')
    pistar     ${\pi^{\ast}}$               (long_name='Foreign Inflation rate')
    istist     ${i^{\ast}}$                 (long_name='Foreign Nominal Interest rate')
    stst        ${S}$                       (long_name='Nominal Exchange Rate')
    thth       ${T^{H}}$                    (long_name='Relative Home Price')
    tftf       ${T^{F}}$                    (long_name='Relative Foreign Price') 
    pitpit     ${\pi}$                      (long_name='CPI Inflation rate')
    ztzt       ${Z}$                        (long_name='Price Dispersion')
    demshock   ${D_S}$                      (long_name='Preference Shock')
@#if Calvo == 1
    pihpih     ${\pi^{H}}$                  (long_name='PPI Inflation rate')
    vnvn       ${V^{N}}$                    (long_name='Aux 1')
    vdvd       ${V^{D}}$                    (long_name='Aux 2')
    phtilde    ${\tilde{p}}$                (long_name='Optimal Price')
    log_pitpit ${\log{\pi}}$
    log_pihpih ${\log{\pi^{H}}}$
@#endif
@#if ngdp == 1
    gtgt     ${G}$
@#endif
@#if StrictG == 1
    gtgt     ${G}$
@#endif

;

varexo  eps_a   ${\varepsilon_a}$   (long_name='technology shock')
        eps_psi ${\varepsilon_{\psi}}$   (long_name='portfolio flows shock')
        eps_z
;

parameters      BBETA       ${\beta}$                   (long_name='discount factor')
                GGAMMA      ${\gamma}$                  (long_name='home bias')
                CCHI        ${\chi}$                    (long_name='inverse Frisch elasticity')
                EPS         ${\epsilon}$                (long_name='elasticity substitution varieties')
                GGAMMAST    ${\gamma^{\ast}}$           (long_name='foreign home bias')
                MMBAR       ${\bar{M}}$                 (long_name='dealers risk capacity')
                OOMEGA      ${\omega}$                  (long_name='CARA risk aversion')
                SSIGMA      ${\sigma}$                  (long_name='exchange rate risk')
                RHOA        ${\rho{a}}$                 (long_name='AR(1) parameter')
                RHONST      ${\rho{n^{\ast}}}$          (long_name='AR(1) parameter')
                STD_PSI     ${\sigma_{\psi}}$           (long_name='st. dev. portfolio shock')
                STD_A       ${\sigma_{a}}$              (long_name='st. dev. technology shock')
                TAUH        ${\tau^H}$                  (long_name='optimal subsidy')
                PPHI_C      ${\phi_c}$                  (long_name='Welfare consumption fraction')
                GGAMMAC     ${\gamma^c}$                (long_name='Inverse of the IES')
                TTHETAH     ${\theta_H}$                (long_name='Price stickyness')
                EPSH        ${\varepsilon_H}$           (long_name='elasticity substitution HF at home')
                EPSF        ${\varepsilon_F}$           (long_name='elasticity substitution HF at foreign')
                STD_Z 
                RHOZ
                B0
                B1
                B2
                B3
                @#if Calvo == 1
                    PHI_PIH    ${\phi_{\pi^H}}$   
                @#endif
                
                @#if ngdp == 1   
                    PHI_G    ${\phi_G}$                        
                @#endif
                @#if Xgap == 1   
                    PHI_Y    ${\phi_G}$                        
                @#endif
 ;

load PARAM_M1;
set_param_value('BBETA',BBETA)
set_param_value('GGAMMA',GGAMMA)
set_param_value('CCHI',CCHI)
set_param_value('EPS',EPS)
set_param_value('GGAMMAST',GGAMMAST)
set_param_value('MMBAR',MMBAR)
set_param_value('OOMEGA',OOMEGA)
set_param_value('SSIGMA',SSIGMA)
set_param_value('RHOA',RHOA)
set_param_value('RHOZ',RHOZ)
set_param_value('RHONST',RHONST)
set_param_value('STD_PSI',STD_PSI)
set_param_value('STD_A',STD_A)
set_param_value('STD_Z',STD_Z)
set_param_value('TAUH',TAUH)
set_param_value('GGAMMAC',GGAMMAC)
set_param_value('TTHETAH',TTHETAH)
set_param_value('EPSH',EPSH)
set_param_value('EPSF',EPSF)
set_param_value('B0',B0)
set_param_value('B1',B1)
set_param_value('B2',B2)
set_param_value('B3',B3)

@#if Calvo == 1
    set_param_value('PHI_PIH',PHI_PIH)
@#endif

@#if ngdp == 1   
    set_param_value('PHI_G',PHI_G)                       
@#endif
@#if Xgap == 1   
    set_param_value('PHI_Y',PHI_Y)                       
@#endif

PPHI_C = 0; % This is the fraction of natural consumption required to make agent indifferent.
@#if logutility == 0
    TAUH = 1/EPS;
@#endif



%----------------------------------------------------------------
% Model
%----------------------------------------------------------------

model;
    @#if Calvo == 1

        @#if dixit == 1
            [name='35. Total Inflation']
            pitpit^(1-EPSH) = GGAMMA*(pihpih^(1-EPSH)) + (1-GGAMMA) * ((dep*pistar)^(1-EPSH));
        @#else
            [name='35. Total Inflation']
            pitpit = (pihpih^(GGAMMA)) * ((dep*pistar)^(1-GGAMMA));

        @#endif

        [name='37. Home goods inflation']
        phtilde = ((1-TTHETAH*(pihpih)^(EPS-1))/(1-TTHETAH))^(1/(1-EPS));   

        [name='38. Relative optimal price phstar/phph']
        phtilde = vnvn / vdvd; 

        [name='39. Domestic goods inflation Num (K)']
        vnvn =  (EPS/(EPS-1)) * mcmc / thth +  BBETA * demshock* TTHETAH * pihpih(+1)^(EPS) * vnvn(+1);

        [name='40. Domestic goods inflation Num (F)']
        vdvd = 1 + BBETA* demshock* TTHETAH * pihpih(+1)^(EPS-1)*vdvd(+1);

        [name='43. Price dispersion']
        ztzt =  (1-TTHETAH)*  (phtilde)^(-EPS) + TTHETAH* pihpih^(EPS)*ztzt(-1);  

        [name='45. Real rate']
        1+rtrt = (1+itit)/pitpit(+1);
        
        [name='48. Definition of log inflation']
        log_pitpit = log(pitpit);
        
        [name='49. Definition of log home inflation']
        log_pihpih = log(pihpih);
        

                
             @#if taylorPPI == 1
             [name='44. Taylor Rule']    
                    itit = (1/BBETA-1) + PHI_PIH*(pihpih-1) ;
                   %// itit = (1/BBETA-1) + PHI_PIH*(pihpih-1) ;

             @#endif

             @#if taylorCPI == 1
                 [name='44. Taylor Rule']  
                    itit = (1/BBETA-1) + PHI_PIH*(pitpit-1);
                   %// itit = (1/BBETA-1) + PHI_PIH*(pitpit-1);

             @#endif
             
             
             @#if StrictPPI == 1
                 [name='44. Strict PPI']
                    pihpih = 1;
             @#endif
             
             @#if StrictCPI == 1
                 [name='44. Strict CPI']
                    pitpit = 1;
             @#endif
             
             @#if Xgap == 1
                 [name='44. Strict Output gap']
                    %x_gap = 0;
                    itit = (1/BBETA-1) + PHI_Y*x_gap +PHI_PIH*(pitpit-1);
             @#endif
             
             @#if ngdp == 1
                 [name='50. NGDP growth']
                    gtgt = (pitpit)*(yhyh/yhyh(-1)) -1 ;
                 
                 [name='51. NGDP Rule']
                    itit = (1/BBETA-1) + PHI_G*(gtgt);                    
             @#endif
             
             @#if StrictG == 1
                 [name='50. NGDP growth']
                    gtgt = (pitpit)*(yhyh/yhyh(-1)) -1 ;
                 
                 [name='51. Strict NGDPT']
                    gtgt = 0;                    
             @#endif

    @#else

        [name='23.rtrt Real rate']
        itit = rtrt;

        [name='24. Inflation']
        pitpit = 1;

        [name='14. Relative optimal price']
        thth = EPS/(EPS-1) * mcmc; 

        [name='16. Price dispersion']
        ztzt = 1; 
        
    @#endif
    
    [name='1. Real exchange rate']
    qtqt = stst * pistar;

    [name='2. Terms of trade']
    totot = qtqt/thth;

    [name='3. Intertemporal Euler']
    ctct^(-GGAMMAC) = BBETA*(demshock(+1)/demshock)*(1 + itit)* ((ctct(+1)^(-GGAMMAC))/pitpit(+1));
    // ctct^(-GGAMMAC) = BBETA*(1 + itit)* ((ctct(+1)^(-GGAMMAC))/pitpit);

    [name='4. Intratemporal Euler']
    ctct^GGAMMAC * ltlt^CCHI = wtwt;

    [name='5. Domestic goods demand']
    chch = GGAMMA * (thth)^(-EPSH) * ctct;

    [name='6. Imports demand']
    cfcf = (1-GGAMMA) * (qtqt)^(-EPSH) * ctct;

    [name='7. Exports demand']
    xtxt = (1-GGAMMAST) * (thth/qtqt)^(-EPSF) * cstcst;

    [name='8. Modified UIP']
    stst(+1)  = stst * (1 + itit)/(1+istist) + OOMEGA /(MMBAR) * (SSIGMA^2) * (dstdst);
    // log(stst(+1)) = log(stst) + log(1+itit) - log(1+istist) + OOMEGA/(MMBAR) * (SSIGMA^2) * (dstdst);
    
    [name='8B. Depreciation']
    dep = stst/stst(-1);

   [name='9. Technology']
    yhyh = atat * ltlt/ztzt;
    
    [name='52. Natural Output']
    ynyn = atat*ltlt;
    
    [name='54. Output gap']
    x_gap = yhyh - ynyn;

    [name='10. Home goods markets equilibrium']
    // ztzt* yhyh = chch + xtxt;
    yhyh = chch + xtxt;

    [name='11. Marginal cost']
    // mcmc = (EPS-1)/EPS * 1/GGAMMA * wtwt/(atat);
    mcmc = (1-TAUH) * wtwt/atat;

   [name='12. Financial Account']
    cca = stst*(dstdst - dstdst(-1) + bcsbcs - bcsbcs(-1) + nstnst - nstnst(-1));
    
    [name='13. Current Account']
    cca = nxnx + istist(-1) * stst * (dstdst(-1) + bcsbcs(-1) + nstnst(-1));
    // cca = nxnx + ((1+istist(-1))/pistar  -1) * stst * (dstdst(-1) + bcsbcs(-1)) +  ( (1+itit(-1))/pitpit -1 ) * stst(-1) * nstnst(-1)  +  (stst - stst(-1))*(dstdst(-1) + bcsbcs(-1));
    // cca = nxnx + istist(-1) * stst * (dstdst(-1) + bcsbcs(-1)) +   itit(-1) * stst(-1) * nstnst(-1);
    // cca = nxnx + istist(-1)* (dstdst(-1) + bcsbcs(-1)) +   itit * stst(-1) * nstnst(-1)  +   (stst - stst(-1))*(dstdst(-1) + bcsbcs(-1));


    [name='14. Net Exports']
    nxnx = thth*yhyh - ctct;
    
    [name='15. PTF shock']
    //atat = 1 + RHOA*(atat(-1)-1) - STD_A * eps_a;
     log(atat) = RHOA*log(atat(-1)) - STD_A * eps_a;
    // log(atat(+1)) = RHOA*log(atat) + STD_A * eps_a;

    [name='16. Central Bank FX Intervention']
    bcsbcs = B0 * (atat-1) - B1 * nstnst - B2 * (stst(+1)/stst - 1) - B3 * (qtqt - 1);
 
    [name='17. Noise traders']
    nstnst = RHONST*nstnst(-1) + STD_PSI * eps_psi;

    [name='18. Foreign demand']
        @#if logutility == 1
            cstcst = GGAMMA^(1/(1+CCHI)) ;
        @#else
            cstcst = 1 ;
        @#endif

    [name='19. Foreign interest rate']
    (1+istist) =(1+rstrst)/pistar;

    [name='20. UIP wedge']
    ctct = wedge*qtqt*cstcst;

    [name='21. Pesos bond market equilibrium']
    btbt + bcsbcs + nstnst + dstdst = 0;
    
    [name='22. Natural consumption']
    cnat = steady_state(ctct);
  
    [name='23. Recursive Welfare']
    @#if ramsey_policy==0    
        welf = util + BBETA * welf(+1);
    @#endif

    @#if logutility == 1
        [name='24. Period utility']
        util = demshock*(log(ctct) - (ltlt^(1+CCHI))/(1+CCHI));

        [name= '25. Natural Recursive Welfare']
        welfeq = demshock*(log( (1-PPHI_C) * cnat ) - ((steady_state(ltlt))^(1+CCHI))/(1+CCHI)) + BBETA * welfeq(+1);
    @#else 
        [name='24. Period utility']
        util = (ctct^(1-GGAMMAC))/(1-GGAMMAC) - (ltlt^(1+CCHI))/(1+CCHI);

        [name='25. Natural Recursive Welfare']
        welfeq =  (((1-PPHI_C) * cnat)^(1-GGAMMAC))/(1-GGAMMAC) - ((steady_state(ltlt))^(1+CCHI))/(1+CCHI) + BBETA * welfeq(+1);
    @#endif 
  
    @#if ramsey_policy==0    
    [name='26. Welfare gap']
    welfgap = welf - welfeq;
    @#endif

   [name='27. Definition log TFP']
    log_a=log(atat);

    [name='28. Definition log output']
    log_yh = log(yhyh);
    
    [name='53. Definition log natural output']
    log_yn = log(ynyn);

    [name='29. Definition log real wage']
    log_w  = log(wtwt);

    [name='30. Definition log labour']
    log_l  = log(ltlt);

    [name='31. Definition log consumption']
    log_c  = log(ctct);

    [name='32. Definition log wedge']
    log_wedge  = log(wedge);

    [name='33. Definition log real exchange rate']
    log_q  = log(qtqt);

    [name='34. Definition log nominal exchange rate variation']
    log_dep  = log(dep);

    [name='36. Foreign prices']
    pistar = 1;

     @#if dixit == 0 
        [name='41. Relative domestic price '] 
        1 = (thth^GGAMMA)*(tftf^(1-GGAMMA));
     @#else 
        [name='41. Relative domestic price '] 
        1 = GGAMMA*thth^(1-EPSH) + (1-GGAMMA)*tftf^(1-EPSH);
     @#endif
    
    [name='46. Foreign goods prices']
    tftf = qtqt;

   [name='47. Real Foreign rate']
    rstrst = 1/BBETA-1;

    [name='48. Preference shock']
    log(demshock) = RHOZ*log(demshock(-1)) + STD_Z*eps_z;
    
end;

%----------------------------------------------------------------
% Steady state values
%----------------------------------------------------------------
steady_state_model;

@#if logutility==1
    ltlt = GGAMMA^(1/(1+CCHI));
    wtwt = GGAMMA^((GGAMMAC+CCHI)/(1+CCHI));
@#else
    ltlt = 1;
    wtwt = 1;
@#endif
    bcsbcs = 0;
    mcmc = (1-TAUH)*wtwt;
    thth = EPS/(EPS-1)*mcmc;
    rtrt = 1/BBETA-1;
    chch = GGAMMA*ltlt;
    cfcf = (1-GGAMMA)*ltlt;
    xtxt = (1-GGAMMAST)*ltlt;
    ctct = ltlt;
    cstcst = ltlt;
    rstrst = 1/BBETA-1;
    yhyh = ltlt;
    ynyn = ltlt;
    x_gap = 0;
    atat = 1;
    qtqt = 1;
    totot = qtqt/thth;
    dep = 1;
    stst = 1;
    tftf = 1;
    cca = 0;
    nxnx = 0;
    nstnst = 0;
    wedge = 1;
    btbt = 0;
    dstdst = - btbt - nstnst - bcsbcs;
    pitpit = 1;
    pihpih = 1;
    pistar = 1;
    itit   = rtrt;
    istist = rstrst;
    ztzt   = 1;
    vnvn   = (EPS/(EPS-1) * mcmc / thth)/(1-BBETA*TTHETAH);
    vdvd   = 1/(1-BBETA*TTHETAH);
    phtilde = vnvn/vdvd;
    demshock = 1;
@#if ngdp==1
    %gtgt = 1;
    gtgt = 0;
@#endif   
@#if StrictG==1
    %gtgt = 1;
    gtgt = 0;
@#endif
@#if logutility == 1   
    util = demshock*(log(ctct) - (ltlt^(1+CCHI))/(1+CCHI));
    cnat = GGAMMA^(GGAMMA/(1+CCHI)) * atat^GGAMMA * cstcst^(1-GGAMMA);
    welfeq  = demshock*((log( (1-PPHI_C)*cnat ) - GGAMMA/(1+CCHI))/(1 - BBETA));
@#else
    util = (ctct^(1-GGAMMAC))/(1-GGAMMAC) - (ltlt^(1+CCHI))/(1+CCHI);
    cnat = 1;
    welfeq  = ((( (1-PPHI_C)*cnat )^(1-GGAMMAC))/(1-GGAMMAC) - (ltlt^(1+CCHI))/(1+CCHI))/(1 - BBETA);
@#endif

@#if ramsey_policy==0
    welf = util/(1-BBETA);
    welfgap = welf - welfeq;
@#endif 
    
    log_a   = log(atat);
    log_yh  = log(yhyh);
    log_w   = log(wtwt);
    log_l   = log(ltlt);
    log_c   = log(ctct);
    log_wedge  = log(wedge);
    log_q   = log(qtqt);
    log_dep = log(dep);
    log_pihpih = log(pihpih);
    log_pitpit = log(pitpit);
    log_yn = log(ynyn);
end;

write_latex_original_model;//(write_equation_tags);
write_latex_static_model;
//write_latex_steady_state_model;
// write_latex_parameter_table;

steady;
check;

%----------------------------------------------------------------
% Do All shocks and OSR
%----------------------------------------------------------------

    shocks;
    %var eps_a = 0;
    var eps_a; stderr 1;
    
    %var eps_psi; stderr 1;
    var eps_psi = 0;
    
    %var eps_z; stderr 1;
    var eps_z = 0;
    end; 

@#if Calvo == 0  
    stoch_simul(order=2, nograph, pruning) log_yh log_c log_l log_wedge log_q log_dep nxnx cca welf util bcsbcs dstdst totot thth atat stst nstnst itit;
@#else
    stoch_simul(order=2, nograph, pruning) log_yh log_c log_l log_wedge log_q log_dep nxnx cca welf util bcsbcs dstdst log_pitpit log_pihpih totot thth  atat stst nstnst itit x_gap;
@#endif

    %----------------------------------------------------------------
    % Optimal Simple Rules
    %----------------------------------------------------------------
    % Need to change limits depending of the combinations 0.86, 1.1
 
    
@#if osr == 1
    @#if taylorPPI == 1                                 
        x_start=[2.78]';
        x_opt_name={'PHI_PIH',0.86,Inf};
        options_.nomoments      = 0;
        options_.nofunctions    = 1;
        options_.nograph        = 1;
        options_.verbosity      = 0;
        options_.noprint        = 1;
        options_.TeX            = 0;
        H0      = 1e-2*eye(length(x_start));        % Initial Hessian 
        crit    = 1e-6;                             % Tolerance
        nit     = 10000;                            % Number of iterations 
 
        [fhat,x_opt_hat] = csminwel(@welfare_objective,x_start,H0,[],crit,nit,x_opt_name);
    @#endif
    
    @#if taylorCPI == 1                                 
        x_start=[2.78]';
        x_opt_name={'PHI_PIH',1.06,Inf};
        options_.nomoments      = 0;
        options_.nofunctions    = 1;
        options_.nograph        = 1;
        options_.verbosity      = 0;
        options_.noprint        = 1;
        options_.TeX            = 0;
        H0      = 1e-2*eye(length(x_start));        % Initial Hessian 
        crit    = 1e-6;                             % Tolerance
        nit     = 10000;                            % Number of iterations 
 
        [fhat,x_opt_hat] = csminwel(@welfare_objective,x_start,H0,[],crit,nit,x_opt_name);
    @#endif
    
    @#if ngdp == 1                                 
        x_start=[2.78]';
        x_opt_name={'PHI_G',1.06,Inf};
        options_.nomoments      = 0;
        options_.nofunctions    = 1;
        options_.nograph        = 1;
        options_.verbosity      = 0;
        options_.noprint        = 1;
        options_.TeX            = 0;
        H0      = 1e-2*eye(length(x_start));        % Initial Hessian 
        crit    = 1e-6;                             % Tolerance
        nit     = 10000;                            % Number of iterations 
 
        [fhat,x_opt_hat] = csminwel(@welfare_objective,x_start,H0,[],crit,nit,x_opt_name);
    @#endif
    
    @#if Xgap == 1                                 
        x_start=[2.78]';
        x_opt_name={'PHI_Y',-1.86,Inf};
        options_.nomoments      = 0;
        options_.nofunctions    = 1;
        options_.nograph        = 1;
        options_.verbosity      = 0;
        options_.noprint        = 1;
        options_.TeX            = 0;
        H0      = 1e-2*eye(length(x_start));        % Initial Hessian 
        crit    = 1e-6;                             % Tolerance
        nit     = 10000;                            % Number of iterations 
 
        [fhat,x_opt_hat] = csminwel(@welfare_objective,x_start,H0,[],crit,nit,x_opt_name);
    @#endif
 
 
@#endif
 
 
%----------------------------------------------------------------
% Compute Variances and Consumption Equivalent
%----------------------------------------------------------------
@#if Calvo == 0  
    yh_pos      =   strmatch('log_yh',var_list_ ,'exact');
    c_pos       =   strmatch('log_c',var_list_ ,'exact');
    l_pos       =   strmatch('log_l',var_list_ ,'exact');
    wedge_pos   =   strmatch('log_wedge',var_list_ ,'exact');
    q_pos       =   strmatch('log_q',var_list_ ,'exact');
    dep_pos     =   strmatch('log_dep',var_list_ ,'exact');
 
    %read out variances
    variance.yh     =   oo_.var(yh_pos,yh_pos);
    variance.c      =   oo_.var(c_pos,c_pos);
    variance.l      =   oo_.var(l_pos,l_pos);
    variance.wedge  =   oo_.var(wedge_pos,wedge_pos);
    variance.q      =   oo_.var(q_pos,q_pos);
    variance.dep    =   oo_.var(dep_pos,dep_pos);
 
 
    % Compute consumption equivalent
    % ==============================
    options_old=options_;
    options_.nocorr=1;
    options_.noprint=1;
    lambda_unconditional_all=csolve('get_consumption_equivalent_unconditional_welfare',0,[],1e-8,1000);
    lambda_conditional_all=csolve('get_consumption_equivalent_conditional_welfare',lambda_unconditional_all,[],1e-8,1000);
    options_=options_old;
   
    % display results
    labels={'sigma(yh)'; 'sigma(c)'; 'sigma(l)'; 'sigma(wedge)'; 'sigma(q)'; 'CE unc'; 'CE cond'};
    headers={'All shocks'};
    values_all = [sqrt([variance.yh; variance.c; variance.l; variance.wedge; variance.q]); 100*lambda_unconditional_all; 100*lambda_conditional_all];
    options_.noprint= 0;
    %dyntable(options_,table_title,headers,labels,100*values_all,size(labels,2)+2,4,3)
 
    stoch_simul(order=2, irf = 30, nograph) log_yh log_c log_l log_wedge log_q log_dep nxnx cca welf util bcsbcs dstdst totot thth yhyh totot ltlt ctct atat stst nstnst rtrt x_gap;
    
 
@#else
 
    yh_pos      =   strmatch('log_yh',var_list_ ,'exact');
    c_pos       =   strmatch('log_c',var_list_ ,'exact');
    l_pos       =   strmatch('log_l',var_list_ ,'exact');
    wedge_pos   =   strmatch('log_wedge',var_list_ ,'exact');
    q_pos       =   strmatch('log_q',var_list_ ,'exact');
    dep_pos     =   strmatch('log_dep',var_list_ ,'exact');
    pihpih_pos  =   strmatch('log_pihpih',var_list_ ,'exact');
    pitpit_pos  =   strmatch('log_pitpit',var_list_ ,'exact');
 
    %read out variances
    variance.yh     =   oo_.var(yh_pos,yh_pos);
    variance.c      =   oo_.var(c_pos,c_pos);
    variance.l      =   oo_.var(l_pos,l_pos);
    variance.wedge  =   oo_.var(wedge_pos,wedge_pos);
    variance.q      =   oo_.var(q_pos,q_pos);
    variance.dep    =   oo_.var(dep_pos,dep_pos);
    variance.pih    =   oo_.var(pihpih_pos,pihpih_pos);
    variance.pit    =   oo_.var(pitpit_pos,pitpit_pos);
 
    % Compute consumption equivalent
    % ==============================
    options_old=options_;
    options_.nocorr=1;
    options_.noprint=1;
    lambda_unconditional_all=csolve('get_consumption_equivalent_unconditional_welfare',0,[],1e-8,1000);
    lambda_conditional_all=csolve('get_consumption_equivalent_conditional_welfare',lambda_unconditional_all,[],1e-8,1000);
    options_=options_old;
        
    % display results
    labels={'sigma(yh)'; 'sigma(c)'; 'sigma(l)'; 'sigma(wedge)'; 'sigma(dep)'; 'sigma(pi)'; 'sigma(pih)' ; 'CE unc'; 'CE cond'};
    headers={'All shocks'};
    values_all = [sqrt([variance.yh; variance.c; variance.l; variance.wedge; variance.dep; variance.pit; variance.pih]); 100*lambda_unconditional_all; 100*lambda_conditional_all];
    options_.noprint= 0;
    %dyntable(options_,table_title,headers,labels,100*values_all,size(labels,2)+2,4,3)
 
    stoch_simul(order=2, irf = 15, nograph) log_yh log_c log_l log_wedge log_q log_dep nxnx cca welf util bcsbcs dstdst log_pitpit log_pihpih totot yhyh ltlt ctct thth atat stst nstnst x_gap itit;
 @#endif
 