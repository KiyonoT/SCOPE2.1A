function biochem_out = biochemical_MD12(leafbio,meteo,options,constants,fV,Q)

%[A,Ci,eta] = biochemical_VCM(Cs,Q,T,eb,O,p,Vcmo,m,Type,Rdparam,stress,Tyear,beta,qLs,NPQs)
% Date:     21 Sep 2012
% Update:   28 Jun 2013 Adaptation for use of Farquhar model of C3 photosynthesis (Farquhar et al 1980)
%           18 Jul 2013 Inclusion of von Caemmerer model of C4 photosynthesis (von Caemmerer 2000, 2013)
%           15 Aug 2013 Modified computation of CO2-limited electron transport in C4 species for consistency with light-limited value
%           22 Oct 2013 Included effect of qLs on Jmax and electron transport; value of kNPQs re-scaled in input as NPQs
%           08 Jan 2019 (CvdT): minor modification to adjust to SCOPE_lite
% Authors:  Federico Magnani, with contributions from Christiaan van der Tol
%
% This function calculates:
%    - CO2 concentration in intercellular spaces (umol/mol == ppmv)
%    - leaf net photosynthesis (umol/m2/s) of C3 or C4 species
%    - fluorescence yield of a leaf (fraction of reference fluorescence yield in dark-adapted and un-stressed leaf)
%
% Usage:
% function [A,Cs,eb,f,rcw] = biochemical(C,Cs,Q,T,ea,eb,O,p,Vcmo,gcparam,Type,tempcor,ra,Tparams,Rdparam,stressfactor,Tyear,beta,qLs,NPQs)
% the function was tested for Matlab 7.2.0.232 (R2006a)
%
% Input (units are important; when not otherwise specified, mol refers to mol C):
% Cs        % [umol/mol]            CO2 concentration at leaf surface
% Q         % [uE/m2/s]             photochemically active radiation absorbed by the leaf
% T         % [oC or K]             leaf temperature
% eb        % [hPa]                 vapour pressure in leaf boundary layer
% O         % [mmol/mol]            ambient O2 concentration
% p         % [Pa]                  air pressure
% Vcmo      % [umol/m2/s]           maximum carboxylation capacity
% m         % [mol/mol]             Ball-Berry coefficient 'm' for stomatal regulation
% Type      % []                    text parameter, either 'C3' for C3 or any other text for C4
% Rdparam   % [mol/mol]             respiration at reference temperature as fraction of Vcmax
% stress    % []                    optional input: stress factor to reduce Vcmax (for example soil moisture, leaf age). Default value = 1 (no stress).
% Tyear     % [oC]                  mean annual temperature
% beta      % []                    fraction of photons partitioned to PSII (0.507 for C3, 0.4 for C4; Yin et al. 2006; Yin and Struik 2012)
% qLs       % []                    fraction of functional reaction centres (Porcar-Castell 2011)
% NPQs      % [s-1]                 rate constant of sustained thermal dissipation, normalized to (kf+kD) (=kNPQs'; Porcar-Castell 2011)
%
% Note: always use the prescribed units. Temperature can be either oC or K
% Note: input can be single numbers, vectors, or n-dimensional matrices
% Note: For consistency reasons, in C4 photosynthesis electron transport rates under CO2-limited conditions are computed by inverting the equation 
%  applied for light-limited conditions(Ubierna et al 2013). A discontinuity would result when computing J from ATP requirements of Vp and Vco, as a 
%  fixed electron transport partitioning is assumed for light-limited conditions

%
% Output:
% A         % [umol/m2/s]           net assimilation rate of the leaves
% Ci        % [umol/mol]            CO2 concentration in intercellular spaces (assumed to be the same as at carboxylation sites in C3 species)
% eta       % []                    amplification factor to be applied to PSII fluorescence yield spectrum
%                                   relative to the dark-adapted, un-stressed yield calculated with either Fluspect or FluorMODleaf

%---------------------------------------------------------------------------------------------------------
%% Start

p       = meteo.p.*1e2;
BallBerrySlope       = leafbio.BallBerrySlope;
BallBerry0          = leafbio.BallBerry0;
O       = meteo.Oa;
Type    = leafbio.Type;
Tyear   = max(11, min(leafbio.Tyear, 35));
beta    = leafbio.beta;
qLs     = leafbio.qLs;
NPQs    = leafbio.kNPQs;
stress  = leafbio.stressfactor;
Cs      = meteo.Cs;
if nargin<6
    Q   = meteo.Q;
end
T       = meteo.T;
eb      = meteo.eb;
Vcmax25    = fV.*leafbio.Vcmax25;
Rdparam = leafbio.RdPerVcmax25;
Q = max(Q, 1E-9);
model = options.stomatal_model;

%% Global and site-specific constants
R             =  constants.R;                         % [J/K/mol]     universal gas constant
rhoa    = constants.rhoa;                             % [kg m-3]      specific mass of air
Mair    = constants.Mair;                             % [g mol-1]     molecular mass of dry air
if options.calc_rhoa
    rhoa = 100*(meteo.p-0.378*eb)./(R/Mair*1000)./(meteo.T+273.15);
end

%---------------------------------------------------------------------------------------------------------
%% Unit conversion and computation of environmental variables
T       = T+273.15*(T<100);                           % [K]           convert temperatures to K if not already
switch model
    case 0 % BWB
        RH = min(1, eb./satvap(T-273.15));            % []         relative humidity (decimal)
    case 1 % Leuning
        D0 = 15;                                      % [hPa] VPD sensitivity (CABLE)
        RH = 1./(1+(satvap(T-273.15)-eb)./D0);        % [] this is not relative humidity!
        % dewing and canopy-water balance are not predicted in SCOPE. here I assume stomatal pore is filled up with water (stomata is closed).
        RH(satvap(T-273.15)<eb) = 0;
    case 2 % Medlyn
        RH = 1./sqrt(0.1*satvap(T-273.15)-0.1*eb);    % [kPa^-0.5] this is not relative humidity!
        RH(satvap(T-273.15)<eb) = 0;
end
Cs      = Cs .* p .*1E-11;                            % [bar]         1E-6 to convert from ppm to fraction, 1E-5 to convert from Pa to bar
O       = O  .* p .*1E-08;                            % [bar]         1E-3 to convert from mmol/mol to fraction, 1E-5 to convert from Pa to bar
ppm2bar =  1e-6 .* (p .*1E-5);

%---------------------------------------------------------------------------------------------------------
%% Define photosynthetic parameters (at reference temperature)
SCOOP     = 2862.;                                    % [mol/mol]     Relative Rubisco specificity for CO2 vs O2 at ref temp (Cousins et al. 2010)
Rdopt     = Rdparam * Vcmax25;                           % [umol/m2/s]   dark respiration at ref temperature from correlation with Vcmax25
switch Type
    case 'C3'                                           % C3 species
        Jmo   =  Vcmax25 * 2.68;                            % [umol/m2/s]   potential e-transport at ref temp from correlation with Vcmax25 (Leuning 1997)
    otherwise                                           % C4 species
        Jmo   =  Vcmax25 * 40/6;                           % [umole-/m2/s] maximum electron transport rate (ratio as in von Caemmerer 2000)
        Vpmo  =  Vcmax25 * 2.33;                             % [umol/m2/s]   maximum PEP carboxylase activity (Yin et al. 2011)
end

%---------------------------------------------------------------------------------------------------------
%% Parameters for temperature corrections
TREF         = 25+273.15;                             % [K]            reference temperature for photosynthetic processes

HARD         = 46.39;                                 % [kJ/mol]       activation energy of Rd
HDRD         = leafbio.TDP.delHdR/1000;               % [kJ/mol]       deactivation energy of Rd (CLM4&5)
DELTASRD     = leafbio.TDP.delSR/1000;                % [kJ/mol/K]     entropy term for Rd (CLM4&5)
CRD          = 1000.*HARD/(R*TREF);                   % []             scaling factor in RD response to temperature

HAGSTAR      = 37.83;                                 % [kJ/mol]       activation energy of Gamma_star
CGSTAR       = 1000.*HAGSTAR/(R*TREF);                % []             scaling factor in GSTAR response to temperature

switch Type
    case 'C3'                                           % C3 species
        HAJ     = 49.88;                                 % [kJ/mol]       activation energy of Jm (Kattge & Knorr 2007)
        % HAJ     = leafbio.TDP.delHaJ/1000; % CLM4 & TB12
        HDJ     = 200;                                   % [kJ/mol]       deactivation energy of Jm (Kattge & Knorr 2007)
        % HDJ     = leafbio.TDP.delHdJ/1000; % CLM4 & TB12
        DELTASJ = (-0.75*Tyear+660)/1000;                % [kJ/mol/K]     entropy term for J  (Kattge and Knorr 2007)
        % DELTASJ = leafbio.TDP.delSJ/1000;  % CLM4 & TB12
        
        HAVCM   = 71.51;                                 % [kJ/mol]       activation energy of Vcm (Kattge and Knorr 2007)
        % HAVCM   = leafbio.TDP.delHaV/1000; % CLM4 & TB12
        HDVC    = 200;                                   % [kJ/mol]       deactivation energy of Vcm (Kattge & Knorr 2007)
        % HDVC    = leafbio.TDP.delHdV/1000; % CLM4 & TB12
        DELTASVC= (-1.07*Tyear+668)/1000;                % [kJ/mol/K]     entropy term for Vcmax (Kattge and Knorr 2007)
        % DELTASVC= leafbio.TDP.delSV/1000;  % CLM4 & TB12
        
        KCOP    = 404.9;                                 % [umol/mol]     Michaelis-Menten constant for CO2 at ref temp (Bernacchi et al 2001)
        HAKC    = 79.43;                                 % [kJ/mol]       activation energy of Kc (Bernacchi et al 2001)
        
        KOOP    = 278.4;                                 % [mmol/mol]     Michaelis-Menten constant for O2  at ref temp (Bernacchi et al 2001)
        HAKO    = 36.38;                                 % [kJ/mol]       activation energy of Ko (Bernacchi et al 2001)
        
    otherwise                                           % C4 species (values can be different as noted by von Caemmerer 2000)
        HAJ    = 77.9;                                   % [kJ/mol]       activation energy of Jm  (Massad et al 2007)
        HDJ     = 191.9;                                 % [kJ/mol]       deactivation energy of Jm (Massad et al 2007)
        DELTASJ = 0.627;                                 % [kJ/mol/K]     entropy term for Jm (Massad et al 2007). No data available on acclimation to temperature.
        
        HAVCM   = 67.29;                                 % [kJ/mol]       activation energy of Vcm (Massad et al 2007)
        HDVC    = 144.57;                                % [kJ/mol]       deactivation energy of Vcm (Massad et al 2007)
        DELTASVC= 0.472;                                 % [kJ/mol/K]     entropy term for Vcm (Massad et al 2007). No data available on acclimation to temperature.
        
        HAVPM   = 70.37;                                 % [kJ/mol]       activation energy of Vpm  (Massad et al 2007)
        HDVP    = 117.93;                                % [kJ/mol]       deactivation energy of Vpm (Massad et al 2007)
        DELTASVP= 0.376;                                 % [kJ/mol/K]     entropy term for Vpm (Massad et al 2007). No data available on acclimation to temperature.
        
        KCOP    = 944.;                                  % [umol/mol]     Michaelis-Menten constant for CO2 at ref temp (Chen et al 1994; Massad et al 2007)
        Q10KC   = 2.1;                                   % []             Q10 for temperature response of Kc (Chen et al 1994; Massad et al 2007)
        
        KOOP    = 633.;                                  % [mmol/mol]     Michaelis-Menten constant for O2 at ref temp (Chen et al 1994; Massad et al 2007)
        Q10KO   = 1.2;                                   % []             Q10 for temperature response of Ko (Chen et al 1994; Massad et al 2007)
        
        KPOP    = 82.;                                   % [umol/mol]     Michaelis-Menten constant of PEP carboxylase at ref temp (Chen et al 1994; Massad et al 2007)
        Q10KP   = 2.1;                                   % []             Q10 for temperature response of Kp (Chen et al 1994; Massad et al 2007)
        
end


%---------------------------------------------------------------------------------------------------------
%% Corrections for effects of temperature and non-stomatal limitations
dum1   = R./1000.*T;                                  % [kJ/mol]
dum2   = R./1000.*TREF;                               % [kJ/mol]

Rd     = Rdopt.*exp(CRD-HARD./dum1);                  % [umol/m2/s]    mitochondrial respiration rates adjusted for temperature (Bernacchi et al. 2001)
if strcmp('C3',Type) % CLM5
    Rd     = Rd.*(1.+exp((TREF*DELTASRD-HDRD)./dum2));
    Rd     = Rd./(1.+exp((T.*DELTASRD-HDRD)./dum1));
end
SCO    = SCOOP./exp(CGSTAR-HAGSTAR./dum1);            % []             Rubisco specificity for CO2 adjusted for temperature (Bernacchi et al. 2001)
Gamma  = 0.5*O./SCO;                                  % [bar]          CO2 compensation point in the absence of mitochondrial respiration

Jmax   = Jmo .* exp(HAJ.*(T-TREF)./(TREF*dum1));
Jmax   = Jmax.*(1.+exp((TREF*DELTASJ-HDJ)./dum2));
Jmax   = Jmax./(1.+exp((T.*DELTASJ-HDJ)./dum1));     % [umol e-/m2/s] max electron transport rate at leaf temperature (Kattge and Knorr 2007; Massad et al. 2007)

Vcmax  = Vcmax25 .* exp(HAVCM.*(T-TREF)./(TREF*dum1)) .* stress; % stress factor expresses mesophyll resistance here
Vcmax  = Vcmax.*(1+exp((TREF*DELTASVC-HDVC)/dum2));
Vcmax  = Vcmax./(1+exp((T.*DELTASVC-HDVC)./dum1));    % [umol/m2/s]    max carboxylation rate at leaf temperature (Kattge and Knorr 2007; Massad et al. 2007)

switch Type
    case 'C3'                                           % C3 species
        CKC    = 1000.*HAKC/(R*TREF);                     % []             scaling factor in KC response to temperature
        Kc     = KCOP.*exp(CKC-HAKC./dum1).*1e-11.*p;     % [bar]          Michaelis constant of carboxylation adjusted for temperature (Bernacchi et al. 2001)
        
        CKO    = 1000.*HAKO/(R*TREF);                     % []             scaling factor in KO response to temperature
        Ko     = KOOP.*exp(CKO-HAKO./dum1).*1e-8.*p;      % [bar]          Michaelis constant of oxygenation adjusted for temperature (Bernacchi et al. 2001)
        Vpmax  = 0; % dummy variable, to be consistent with C4
    otherwise                                           % C4 species
        Vpmax  = Vpmo .* exp(HAVPM.*(T-TREF)./(TREF*dum1));
        Vpmax  = Vpmax.*(1+exp((TREF*DELTASVP-HDVP)/dum2));
        Vpmax  = Vpmax./(1+exp((T.*DELTASVP-HDVP)./dum1));% [umol/m2/s]    max carboxylation rate at leaf temperature (Massad et al. 2007)
        
        Kc     = KCOP.*Q10KC .^ ((T-TREF)/10.)*1e-11*p;    % [bar]          Michaelis constant of carboxylation temperature corrected (Chen et al 1994; Massad et al 2007)
        
        Ko     = KOOP.*Q10KO .^ ((T-TREF)/10.)*1e-8*p;     % [bar]          Michaelis constant of oxygenation  temperature corrected (Chen et al 1994; Massad et al 2007)
        
        Kp     = KPOP.*Q10KP .^ ((T-TREF)/10.)*1e-11*p;    % [bar]          Michaelis constant of PEP carboxyl temperature corrected (Chen et al 1994; Massad et al 2007)
        MM_consts.Kp = Kp;
end
MM_consts.Kc = Kc; % Michaelis-Menten constants
MM_consts.Ko = Ko;

%---------------------------------------------------------------------------------------------------------
%% Define electron transport and fluorescence parameters
kf        = 6.7E7;                                    % [s-1]         rate constant for fluorescence (Rabinowich & Govindjee 1969, p136)
kD        = 6.03E8;                                   % [s-1]         rate constant for thermal deactivation at Fm (Porcar-Castell et al. 2006)
kd        = 1.005E9;                                  % [s-1]         rate constant of energy dissipation in closed RCs (for theta=0.7 under un-stressed conditions)  
po0max    = 0.88;                                     % [mol e-/E]    maximum PSII quantum yield, dark-acclimated in the absence of stress (Pfundel 1998)
kPSII     = (kD+kf) * po0max/(1.-po0max);             % [s-1]         rate constant for photochemisty (Genty et al. 1989)
fo0       = kf./(kf+kPSII+kD);                        % [E/E]         reference dark-adapted PSII fluorescence yield under un-stressed conditions

kps       = kPSII * qLs;                              % [s-1]         rate constant for photochemisty under stressed conditions (Porcar-Castell 2011)
kNPQs     = NPQs * (kf+kD);                           % [s-1]         rate constant of sustained thermal dissipation (Porcar-Castell 2011)
kds       = kd * qLs;
kDs       = kD + kNPQs;
Jms       = Jmax * qLs;                               % [umol e-/m2/s] potential e-transport rate reduced for PSII photodamage
po0       = kps ./(kps+kf+kDs);                       % [mol e-/E]    maximum PSII quantum yield, dark-acclimated in the presence of stress
THETA     = (kps-kds)./(kps+kf+kDs);                  % []            convexity factor in J response to PAR

%---------------------------------------------------------------------------------------------------------
%% Calculation of electron transport rate
Q2     = beta * Q * po0;
J      = (Q2+Jms-sqrt((Q2+Jms).^2-4*THETA.*Q2.*Jms))./(2*THETA); % [umol e-/m2/s]    electron transport rate under light-limiting conditions

%---------------------------------------------------------------------------------------------------------
%% Calculation of Ci & net photosynthesis
switch Type
    case 'C3'
        minCi = 0.3;
    otherwise
        minCi = 0.1;
end
computeA()  % clears persistent fcount
computeA_fun = @(x) computeA(x, Type, MM_consts, Rd, Vcmax, Gamma, O, J, Vpmax);

Ci = BallBerry(model, Cs, RH, [], BallBerrySlope, BallBerry0, Gamma, minCi); % initial value

if all(BallBerry0 == 0)
    % b = 0: no need to iterate:
    %     A =  computeA_fun(Ci);   
else
    % compute Ci using iteration (JAK)
    % it would be nice to use a built-in root-seeking function but fzero requires scalar inputs and outputs,
    % Here I use a fully vectorized method based on Brent's method (like fzero) with some optimizations.
    tol = 1e-7;  % 0.1 ppm more-or-less
    % Setting the "corner" argument to Gamma may be useful for low Ci cases, but not very useful for atmospheric CO2, so it's ignored.
    %                     (fn,                           x0, corner, tolerance)
    [Ci] = fixedp_brent_ari(@(x) Ci_next(x, Cs, RH, minCi, model, BallBerrySlope, BallBerry0, Gamma, computeA_fun, ppm2bar), Ci, [], tol); % [] in place of Gamma: it didn't make much difference
    %NOTE: A is computed in Ci_next on the final returned Ci. fixedp_brent_ari() guarantees that it was done on the returned values.
    %     A =  computeA_fun(Ci);
end

[A, biochem_out]    = computeA_fun(Ci);
Ja = biochem_out.Ja;

%---------------------------------------------------------------------------------------------------------
%% Calculation of PSII quantum yield and fluorescence
ps     = Ja ./(beta.*Q);                            % [mol e-/E]    PSII photochemical quantum yield
[fs,phase]   = MD12(ps,Ja,Jms,kps,kf,kds,kDs);      % [E/E]         PSII fluorescence yield
eta    = fs./fo0;                                   % []            scaled PSII fluorescence yield
qLr = ps./fs * kf/kPSII / qLs; % fraction of open reaction centres
kp  = kps * qLr;
Phi_N  = 1 - ps - fs - fs*kD/kf;

%rcw         = 0.625*(Cs-Ci)./A *rhoa/Mair*1E3    * 1e6 ./ p .* 1E5;
gs = 1.6 * A* ppm2bar./ (Cs-Ci);
rcw      =  (rhoa./(Mair*1E-3))./gs;
rcw(A<=0 & rcw~=0)   = 0.625*1E6;

%% convert back to ppm
biochem_out.CiCa = Ci./Cs;
Ci          = Ci*1e6 ./ p .* 1E5;

%%
biochem_out.Ci = Ci;
biochem_out.ps = ps;
biochem_out.eta = eta;
biochem_out.fs  = fs;
biochem_out.rcw = rcw;
biochem_out.qE  = rcw*NaN; % dummy output, to be consistent with SCOPE
biochem_out.Kn  = Phi_N.*(kp+kf+kD)./(1-Phi_N); %
biochem_out.Phi_N = Phi_N;
biochem_out.phase = phase;

return;

end
%%% end of function biochemical


%---------------------------------------------------------------------------------------------------------
%% MD12 algorithm for the computation of fluorescence yield

function [fs, phase] = MD12(ps,Ja,Jms,kps,kf,kds,kDs)

fs1    = ps .* (kf./kps) ./ (1. - Ja./Jms);         % [E/E]   PSII fluorescence yield under CO2-limited conditions

par1   = kps./(kps-kds);                            % [E/E]   empirical parameter in the relationship under light-limited conditions
par2   = par1.* (kf+kDs+kds)./kf;                  % [E/E]   empirical parameter in the relationship under light-limited conditions
fs2    = (par1-ps)./par2;                           % [E/E]   PSII fluorescence yield under light-limited conditions

fs     = min(fs1,fs2);                              % [E/E]   PSII fluorescence yield
phase  = fs1 < fs2;
end



%% Ball Berry Model
function [Ci, gs] = BallBerry(model, Cs, RH, A, BallBerrySlope, BallBerry0, Gamma, minCi, Ci_input)
%  Cs  : CO2 at leaf surface
%  RH  : relative humidity
%  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
% BallBerrySlope, BallBerry0, 
% minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
% Ci_input : will calculate gs if A is specified.
if nargin > 8 && ~isempty(Ci_input)
    % Ci is given: try and compute gs
    Ci = Ci_input;
    gs = [];
    if ~isempty(A) && nargout > 1
        gs = gsFun(model, Cs, RH, A, BallBerrySlope, BallBerry0, Gamma);
    end
elseif all(BallBerry0 == 0) || isempty(A)
    % EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
    %  so Ci = Cs - A/gs (at equilibrium)                                 [2]
    %  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
    %  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
    %  Substituting into [2]
    %  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
    Ci      = max(minCi .* Cs,  Cs.*(1-1.6./(BallBerrySlope .* RH)));
    gs      = [];
else
    %  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
    % if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
    % note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
    %   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
    %  don't let gs go below its minimum value (i.e. when A goes negative)
    gs = gsFun(model, Cs, RH, A, BallBerrySlope, BallBerry0, Gamma);
    Ci = max(minCi .* Cs,  Cs - 1.6 * A./gs) ;
end

end % function

function gs = gsFun(model, Cs, RH, A, BallBerrySlope, BallBerry0, Gamma)
% add in a bit just to avoid div zero. 1 ppm = 1e-6 (note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
switch model
    case 0 % BWB
        gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs+1e-9)  + BallBerry0);
    case 1 % Leuning
        gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs-Gamma) + BallBerry0);
    case 2 % Medlyn
        gs = max(BallBerry0, (1+BallBerrySlope.*RH).* A ./ (Cs+1e-9) + BallBerry0);
end
% clean it up:
%gs( Cs == 0 ) = would need to be max gs here;  % eliminate infinities
gs( isnan(Cs) ) = NaN;  % max(NaN, X) = X  (MATLAB 2013b) so fix it here
end

%% Test-function for iteration
%   (note that it assigns A in the function's context.)
%   As with the next section, this code can be read as if the function body executed at this point.
%    (if iteration was used). In other words, A is assigned at this point in the file (when iterating).
function [err, Ci_out] = Ci_next(Ci_in, Cs, RH, minCi, model, BallBerrySlope, BallBerry0, Gamma, A_fun, ppm2bar)
% compute the difference between "guessed" Ci (Ci_in) and Ci computed using BB after computing A
A = A_fun(Ci_in);
A_bar = A .* ppm2bar;
Ci_out = BallBerry(model, Cs, RH, A_bar, BallBerrySlope, BallBerry0, Gamma, minCi); %[Ci_out, gs]

err = Ci_out - Ci_in; % f(x) - x
end

%% Compute Assimilation.
%  Note: even though computeA() is written as a separate function,
%    the code is, in fact, executed exactly this point in the file (i.e. between the previous if clause and the next section
function [A, biochem_out] = computeA(Ci, Type, MM_consts, Rd, Vcmax, Gamma, O, J, Vpmax)
persistent fcount
if nargin == 0
    fcount = 0;
    return
end
Kc = MM_consts.Kc;
Ko = MM_consts.Ko;
switch Type
    case 'C3'                                           % C3 species, based on Farquhar model (Farquhar et al. 1980)
        Cc  = Ci;                                        % [bar]             CO2 concentration at carboxylation sites (neglecting mesophyll resistance)        
        Wc  = Vcmax .* Cc ./ (Cc + Kc .* (1+O./Ko));     % [umol/m2/s]       RuBP-limited carboxylation
        Wj  = J.*Cc ./ (4.5*Cc + 10.5*Gamma);            % [umol/m2/s]       electr transp-limited carboxyl
        
        W   = min(Wc,Wj);                                % [umol/m2/s]       carboxylation rate
        Ag  = (1 - Gamma./Cc) .*W;                       % [umol/m2/s]       gross photosynthesis rate
        A   = Ag - Rd;                                   % [umol/m2/s]       net photosynthesis rate
        Ja  = J.*W ./Wj;                                 % [umole-/m2/s]     actual linear electron transport rate

    otherwise                                           % C4 species, based on von Caemmerer model (von Caemmerer 2000)
        Kp = MM_consts.Kp;
        Vpr   =  80;                                     % [umol/m2/s]   PEP regeneration rate, constant (von Caemmerer 2000)
        gbs   =  (0.0207*Vcmax25+0.4806)*1000.;           % [umol/m2/s]   bundle sheath conductance to CO2 (Yin et al. 2011)
        x     =  0.4;                                     % []            partitioning of electron transport to mesophyll (von Caemmerer 2013)
        alpha =  0;                                      % []            bundle sheath PSII activity (=0 in maize/sorghum; >=0.5 in other cases; von Caemmerer 2000)

        %Ci    =  max(9.9e-6*(p*1e-5),Cs.*(1-1.6./(m.*RH*stress)));
        % [bar]             intercellular CO2 concentration from Ball-Berry model (Ball et al. 1987)
        Cm    =  Ci;                                     % [bar]             mesophyll CO2 concentration (neglecting mesophyll resistance)
        Rs    =  0.5 .* Rd;                               % [umol/m2/s]       bundle sheath mitochondrial respiration (von Caemmerer 2000)
        Rm    =  Rs;                                     % [umol/m2/s]       mesophyll mitochondrial respiration
        gam   =  Gamma./O;                               % []                half the reciprocal of Rubisco specificity for CO2
        
        Vpc   = Vpmax .* Cm./(Cm+Kp);                     % [umol/m2/s]       PEP carboxylation rate under limiting CO2 (saturating PEP)
        Vp    = min(Vpc,Vpr);                            % [umol/m2/s]       PEP carboxylation rate
        
        % Complete model proposed by von Caemmerer (2000)
        dum1  =  alpha/0.047;                           % dummy variables, to reduce computation time
        dum2  =  Kc./Ko;
        dum3  =  Vp-Rm+gbs.*Cm;
        dum4  =  Vcmax-Rd;
        dum5  =  gbs.*Kc.*(1+O./Ko);
        dum6  =  gam.*Vcmax;
        dum7  =  x*J./2. - Rm + gbs.*Cm;
        dum8  =  (1.-x).*J./3.;
        dum9  =  dum8 - Rd;
        dum10 =  dum8 + Rd * 7/3;
        
        a     =  1. - dum1 .* dum2;
        b     =  -(dum3+dum4+dum5+dum1.*(dum6+Rd.*dum2));
        c     =  dum4.*dum3-dum6.*gbs*O+Rd.*dum5;
        Ac    =  (-b - sqrt(b.^2-4.*a.*c))./(2.*a);           % [umol/m2/s]       CO2-limited net photosynthesis
        
        a     =  1.- 7./3.*gam.*dum1;
        b     =  -(dum7+dum9 + gbs.*gam.*O.*7./3. + dum1.*gam.*dum10);
        c     =  dum7.*dum9 - gbs.*gam.*O.*dum10;
        Aj    =  (-b - sqrt(b.^2-4.*a.*c))./(2.*a);           % [umol/m2/s]       light-limited net photosynthesis (assuming that an obligatory Q cycle operates)
        
        A     =  min(Ac,Aj);                             % [umol/m2/s]       net photosynthesis
        
        Ja    =  J;                                      % [umole-/m2/s]     actual electron transport rate, CO2-limited
               
        if any(A==Ac) %IPL 03/09/2013
            
            ind=A==Ac;
            a(ind)   =  x.*(1-x)./6./A(ind);
            b(ind)   =  (1-x)/3.*(gbs(ind)./A(ind).*(Cm(ind)-Rm(ind)./gbs(ind)-gam(ind).*O)-1-alpha.*gam(ind)./0.047)-x./2.*(1.+Rd(ind)./A(ind));
            c(ind)   =  (1+Rd(ind)./A(ind)).*(Rm(ind)-gbs(ind).*Cm(ind)-7.*gbs(ind).*gam(ind).*O./3)+(Rd(ind)+A(ind)).*(1-7.*alpha.*gam(ind)./3./0.047);
            Ja(ind)  =  (-b(ind) + sqrt(b(ind).^2-4.*a(ind).*c(ind)))./(2.*a(ind));            % [umole-/m2/s]     actual electron transport rate, CO2-limited
        end
end

fcount = fcount + 1; % # of times we called computeA

if nargout > 1
    biochem_out.A = A;
    biochem_out.Ag = Ag;
    biochem_out.Ja = Ja;
    biochem_out.fcount = fcount;
end

end


% Sources:
%  Ball J. T., I. E. Woodrow and J. A. Berry. (1987) A model predicting stomatal conductance and its contribution to the control of photosynthesis
%    under different environmental conditions. In: Progress in Photosynthesis Research (Ed. J. Biggens), p. 221-224, The Netherlands:Martinus Nijhoff.
%  Bernacchi C.J., E.L. Singsaas, C. Pimentel, A.R. Portis and S.P. Long (2001) Improved temperature response functions for models of Rubisco-limited
%    photosynthesis. Plant Cell Envir 24:253-259.
%  Bernacchi C.J., C. Pimentel and S.P. Long (2003) In vivo temperature response functions of parameters required to model RuBP-limited photosynthesis.
%    Plant Cell Envir 26 (9):1419-1430.
%  Chen D.X., M.B. Coughenour, A.K. Knapp, and C.E. Owensby (1994) Mathematical simulation of C4 grass photosynthesis in ambient and elevated CO2.
%    Ecol.Model. 73:63-80, 1994.
%  Cousins A.B., O. Ghannoum, S. von Caemmerer, and M.R. Badger (2010) Simultaneous determination of Rubisco carboxylase and oxygenase kinetic parameters
%    in Triticum aestivum and Zea mays using membrane inlet mass spectrometry. Plant Cell Envir 33:444-452.
%  Farquhar G.D., S. von Caemmerer and J.A. Berry (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta 149:78-90.
%  Genty B., J.-M. Briantais and N. R. Baker (1989) The relationship between quantum yield of photosynthetic electron transport and quenching of
%    chlorophyll fluorescence. Biochimica et Biophysica Acta 990:87-92.
%  Kattge  J. and W. Knorr (2007) Temperature acclimation in a biochemical model of photosynthesis: a reanalysis of data from 36 species.
%    Plant Cell Envir 30:1176-1190.
%  Leuning R. (1997) Scaling to a common temperature improves the correlation between the photosynthesis parameters Jmax and Vcmax.
%    J.Exp.Bot. 48 (307):345-347.
%  Massad R.S., A. Tuzet and O. Bethenod (2007) The effect of temperature on C4-type leaf photosynthesis parameters. Plant Cell Envir 30:1191-1204.
%  Pfundel E. (1998) Estimating the contribution of Photosystem I to total leaf chlorophyll fluorescence. Photosynthesis Research 56:185-195.
%  Porcar-Castell A. (2011) A high-resolution portrait of the annual dynamics of photochemical and non-photochemical quenching in needles of  Pinus sylvestris.
%    Physiol.Plant. 143:139-153.
%  von Caemmerer S. (2000) Biochemical Models of Leaf Photosynthesis, Canberra:CSIRO Publishing.
%  von Caemmerer S. (2013) Steady-state models of photosynthesis. Plant Cell Envir, in press.
%  Yin X., Z. Sun, P.C. Struik, P.E.L. van der Putten, W. van Ieperen and J. Harbinson (2011) Using a biochemical C4 photosynthesis model and combined
%    gas exchange and chlorophyll fluorescence measurements to estimate bundle-sheath conductance of maize leaves differing in age and nitrogen content.
%    Plant Cell Envir 34:2183-2199.
%
