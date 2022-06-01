function [lE, H, ec, Cc, lambda, s]  = heatfluxes2(ra,rs,Tc,ea,Ta,e_to_q,Ca,Ci,rhoa,cp, es_fun, s_fun, alpha)    
   
    lambda      = (2.501-0.002361*Tc)*1E6;  %      [J kg-1]  Evapor. heat (J kg-1)
    ei = es_fun(Tc);
    s = s_fun(ei, Tc);

    if nargin < 13 % leaf
        alpha = 1;
        rav = ra*0.92; % laminar boundary layer (heat -> H2O)
    else % soil
        rav = ra;      % turbulent boundary layer (heat -> H2O)
    end
    qi = ei .* e_to_q .* alpha; % water-potential correction
    qa = ea .* e_to_q;
    
    rac = ra*1.23;  % laminar boundary layer (heat -> CO2)
    rsc = rs*1.6;  % stomata - i.e. still air (water -> CO2)
    
    lE = rhoa./(rav+rs).*lambda.*(qi-qa); % [W m-2]   Latent heat flux
    H  = (rhoa*cp)./ra.*(Tc-Ta);          % [W m-2]   Sensible heat flux
    ec = ea + (ei-ea).*rav./(rav+rs);     % [hPa] vapour pressure at the leaf surface
    Cc = Ca - (Ca-Ci).*rac./(rac+rsc);    % [ppm] CO2 concentration at the leaf surface
end