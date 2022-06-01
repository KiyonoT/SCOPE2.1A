function [resist_out] = resistances(constants,soil,canopy,meteo,options)
%
%   function resistances calculates aerodynamic and boundary resistances
%   for soil and vegetation
%
%   Date:       01 Feb 2008
%   Authors:    Anne Verhoef            (a.verhoef@reading.ac.uk)
%               Christiaan van der Tol  (tol@itc.nl)
%               Joris Timmermans        (j_timmermans@itc.nl)
%   Source:     Wallace and Verhoef (2000) 'Modelling interactions in
%               mixed-plant communities: light, water and carbon dioxide', in: Bruce
%               Marshall, Jeremy A. Roberts (ed), 'Leaf Development and Canopy Growth',
%               Sheffield Academic Press, UK. ISBN 0849397693
%               
%               ustar:  Tennekes, H. (1973) 'The logaritmic wind profile', J.
%               Atmospheric Science, 30, 234-238
%               Psih:   Paulson, C.A. (1970), The mathematical
%               representation of wind speed and temperature in the
%               unstable atmospheric surface layer. J. Applied Meteorol. 9,
%               857-861
%       
% Note: Equation numbers refer to equation numbers in Wallace and Verhoef (2000)
% 
% Usage:
%   [resist_out] = resistances(resist_in)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input: 
%   resist_in   aerodynamic resistance parameters and wind speed
%
% The strucutre resist_in contains the following elements:
% u         =   windspeed
% L         =   stability
% LAI       =   Leaf Area Index

% rbs       =   Boundary Resistance of soil                         [s m-1]
% rss       =   Surface resistance of soil for vapour transport     [s m-1]
% rwc       =   Within canopy Aerodynamic Resistance canopy         [s m-1]

% z0m       =   Roughness lenght for momentum for the vegetation    [m]
% d         =   Displacement height (Zero plane)                    [m]
% z         =   Measurement height                                  [m]
% h         =   Vegetation height                                   [m]

%
% Output:
%   resist_out  aeorodynamic resistances
%
% The strucutre resist_out contains the following elements:
% ustar     =   Friction velocity                                   [m s-1]
% raa       =   Aerodynamic resistance above the canopy             [s m-1]                     
% rawc      =   Total resistance within the canopy (canopy)         [s m-1]
% raws      =   Total resistance within the canopy (soil)           [s m-1]

% rai       =   Aerodynamic resistance in inertial sublayer         [s m-1]
% rar       =   Aerodynamic resistance in roughness sublayer        [s m-1]
% rac       =   Aerodynamic resistance in canopy layer (above z0+d) [s m-1]

% rbc       =   Boundary layer resistance (canopy)                  [s m-1]
% rwc       =   Aerodynamic Resistance within canopy(canopy)(Update)[s m-1]

% rbs       =   Boundary layer resistance (soil) (Update)           [s m-1]
% rws       =   Aerodynamic resistance within canopy(soil)          [s m-1] 

% rss       =   Surface resistance vapour transport(soil)(Update)   [s m-1]

% uz0       =   windspeed at z0                                     [m s-1]
% Kh        =   Diffusivity for heat                                [m2s-1]

%% parameters
%global constants
kappa   = constants.kappa;
Cd      =  canopy.Cd;
LAI     =  canopy.LAI;
rwc     =  canopy.rwc;
z0m     =  canopy.zo;
d       =  canopy.d;
h       =  canopy.hc;
w       =  canopy.leafwidth;
z       =  meteo.z;
u       =  max(0.3,meteo.u);
L       =  meteo.L;
rbs     =  soil.rbs;
%rss       =  resist_in.rss;

% Canopy aerodynamic resistance for heat
if options.MoninObukhov < 2
    %% SCOPE's built-in function (W&V)
    if options.MoninObukhov==0, L = -1E6; end
    % derived parameters
    %zr: top of roughness sublayer, bottom of intertial sublayer
    zr			= 2.5*h;                   %                            [m]			
    %n: dimensionless wind extinction coefficient                       W&V Eq 33
    n			= Cd*LAI/(2*kappa^2);      %                            [] 

    %% stability correction for non-neutral conditions
    unst        = (L <  0 & L>-500);
    st          = (L >  0 & L<500);  
    x       	= (1-16*z./L).^(1/4); % only used for unstable

    % stability correction functions, friction velocity and Kh=Km=Kv
    pm_z    	= psim(z -d,L,unst,st,x);
    ph_z    	= psih(z -d,L,unst,st,x);
    pm_h        = psim(h -d,L,unst,st,x);
    %ph_h       = psih(h -d,L,unst,st);
    ph_zr       = psih(zr-d,L,unst,st,x).*(z>=zr) + ph_z.*(z<zr);
    phs_zr      = phstar(zr,zr,d,L,st,unst,x);
    phs_h		= phstar(h ,zr,d,L,st,unst,x);

    ustar   	= max(.001,kappa*u./(log((z-d)/z0m) - pm_z));%          W&V Eq 30
    Kh          = kappa*ustar*(zr-d);                  %                W&V Eq 35

    if unst
        resist_out.Kh	= Kh*(1-16*(h-d)./L).^.5;% W&V Eq 35
    elseif st
        resist_out.Kh   = Kh*(1+ 5*(h-d)./L  ).^-1;% W&V Eq 35
    else
        resist_out.Kh = Kh;
    end

    %% wind speed at height h and z0m
    uh			= max(ustar/kappa .* (log((h-d)/z0m) - pm_h     ),.01);
    uz0 		= uh*exp(n*((z0m+d)/h-1));                      %       W&V Eq 32

    %% resistances

    resist_out.uz0 = uz0; 
    rai = (z>zr).*(1./(kappa*ustar).*(log((z-d) /(zr-d))  - ph_z   + ph_zr));% W&V Eq 41 
    rar = 1./(kappa*ustar).*((zr-h)/(zr-d)) 	 - phs_zr + phs_h;% W&V Eq 39
    rac = h*sinh(n)./(n*Kh)*(log((exp(n)-1)/(exp(n)+1)) - log((exp(n*(z0m+ d )/h)-1)/(exp(n*(z0m +d )/h)+1))); % W&V Eq 42
    rws = h*sinh(n)./(n*Kh)*(log((exp(n*(z0m+d)/h)-1)/(exp(n*(z0m+d)/h)+1)) - log((exp(n*(.01    )/h)-1)/(exp(n*(.01    )/h)+1))); % W&V Eq 43
    %rbc = 70/LAI * sqrt(w./uz0);						%		W&V Eq 31, but slightly different

    resist_out.rai = rai;
    resist_out.rar = rar;
    resist_out.rac = rac;
    resist_out.rws = rws;
    %resist_out.rbc = rbc;
    resist_out.RiB = (z-d)/L * (log((z-d)/z0m)-ph_z) * (log((z-d)/z0m)-pm_z)^-2;

    raa  = rai + rar + rac;
    rawc = rwc;% + rbc;
    raws = rws + rbs;
else
    %% CLM5 
    % Equation numbers refer to those in CLM5 tech note
    % https://www.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
    % moisture-induced air-density change is omitted
    g   = constants.g;
    Va  = meteo.Va;
    rv  = 0.622*meteo.ea/(meteo.p-meteo.ea); % mixing ratio
    Ta  = (meteo.Ta + 273.15*(meteo.Ta<100))*(1+0.61*rv); % virtual temperature
    Tc  = canopy.Tave + 273.15*(canopy.Tave<100);
    Ts  = soil.Tave + 273.15*(soil.Tave<100);
    rhoa = 100*meteo.p/(constants.R/constants.Mair*1000)/(meteo.Ta+273.15);
    nu  = 1.821E-5*((293.15+117)/(Ta+117))*(Ta/293.15)^1.5 / rhoa; % 1.5E-5 [m2 s-1] @ 20 degC and 1000 hPa
    z0h = z0m; % Canopy roughness length for heat; assumed to be equal to that for momentum (Eq 5.125)
    if L==0, L = Inf; end % to avoid zero division (neutral condition)
    zeta = (z-d)/L;
    %% Stability correction (Eqs 5.32--5.40)
    % Tplus is identical to 'theta* / (theta_atm - theta_s)' in CLM5 tech note
    if zeta < -1.574
        ustar = kappa*Va / (log(-1.574*L/z0m) - psi_m(-1.574) ...
            + 1.14*((-zeta)^(1/3) - 1.574^(1/3)) + psi_m(z0m/L));
        Tplus = kappa / (log(-0.465*L/z0h) - psi_h(-0.465) ...
            + 0.8*(0.465^(-1/3) - (-zeta)^(-1/3)) + psi_h(z0h/L));
    elseif zeta <= -0.465
        ustar = kappa*Va / (log((z-d)/z0m) - psi_m(zeta) + psi_m(z0m/L));
        Tplus = kappa / (log(-0.465*L/z0h) - psi_h(-0.465) ...
            + 0.8*(0.465^(-1/3) - (-zeta)^(-1/3)) + psi_h(z0h/L));
    elseif zeta <= 0
        ustar = kappa*Va / (log((z-d)/z0m) - psi_m(zeta) + psi_m(z0m/L));
        Tplus = kappa / (log((z-d)/z0h) - psi_h(zeta) + psi_h(z0h/L));
    elseif zeta <= 1
        ustar = kappa*Va / (log((z-d)/z0m) + 5*zeta - 5*z0m/L);
        Tplus = kappa / (log((z-d)/z0h) + 5*zeta - 5*z0h/L);
    else
        ustar = kappa*Va / (log(L/z0m) + 5 + 5*log(zeta) + zeta - 1 - 5*z0m/L);
        Tplus = kappa / (log(L/z0h) + 5 + 5*log(zeta) + zeta - 1 - 5*z0h/L);
    end
    
    %% Above-canopy aerodynamic resistance
    raa = 1/Tplus/ustar; % Eq 5.56
    
    %% Bulk leaf boundary layer resistance
    rawc = 100/LAI * sqrt(w./ustar); % Eqs 5.117 & 5.122
    
    %% Soil resistance
    z0mg = 0.01; % Bare-soil roughness length for momentum
    Cs_bare = kappa/0.13*(z0mg*ustar/nu)^(-0.45); % Eq 5.121
    Cs = Cs_bare*exp(-LAI) + 0.004*(1-exp(-LAI)); % Eq 5.118--5.120
    raws = 1 / (Cs*ustar); % Eq 5.116
    %% Calculate within-canopy air temperature and update the mean wind speed (Va)
    Ta_surf = (Ta/raa + Ts/raws + Tc/rawc) / (1/raa + 1/raws + 1/rawc); % Eqs 5.93--5.96 (see Fig 5.1b)
    Tstar = Tplus * (Ta - Ta_surf);
    if Tstar > 0 || L == Inf
        Va = u;
    else
        wstar = (-g*ustar*Tstar*1000/Ta)^(1/3); % convective velocity scale (Eq 5.29)
        beta = 1; % Eq 5.28
        Va = sqrt(u^2 + (beta*wstar)^2); % Eq 5.24
    end
    resist_out.Va = Va;
    resist_out.RiB = zeta * (log((z-d)/z0h)-psi_h(zeta)) * (log((z-d)/z0m)-psi_m(zeta))^-2;
end

resist_out.ustar = ustar;
resist_out.raa  = raa;          % aerodynamic resistance above the canopy
resist_out.rawc	= rawc;			% aerodynamic resistance within the canopy (canopy)
resist_out.raws	= raws;			% aerodynamic resistance within the canopy (soil)

return


%% subfunction pm for stability correction (eg. Paulson, 1970)
function pm = psim(z,L,unst,st,x)
pm      	= 0;
if unst
    pm          = 2*log((1+x)/2)+log((1+x.^2)/2) -2*atan(x)+pi/2;   %   unstable
elseif st
    pm          = -5*z./L;                                      %   stable
end
return

%% subfunction ph for stability correction (eg. Paulson, 1970)
function ph = psih(z,L,unst,st,x)
ph        = 0;
if unst
    ph      = 2*log((1+x.^2)/2);                                %   unstable
elseif st
    ph      = -5*z./L;                                      %   stable
end
return

%% subfunction ph for stability correction (eg. Paulson, 1970)
function phs = phstar(z,zR,d,L,st,unst,x)
phs         = 0;
if unst
    phs     = (z-d)/(zR-d)*(x.^2-1)./(x.^2+1);
elseif st
    phs     = -5*z./L;
end
return

%% Alternative subfunction for pm (Zeng 1998)
function pm = psi_m(zeta)
x = (1-16*zeta).^(1/4);
pm = 2*log((1+x)/2) + log((1+x.^2)/2) -2*atan(x) + pi/2;
return

%% Alternative subfunction for ph (Zeng 1998)
function ph = psi_h(zeta)
x = (1-16*zeta).^(1/4);
ph = 2*log((1+x.^2)/2);
return