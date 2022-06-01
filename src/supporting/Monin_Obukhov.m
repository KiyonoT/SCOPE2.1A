function L = Monin_Obukhov(constants,meteo,H,options)

if options.calc_rhoa
    rhoa = 100*(meteo.p-0.378*meteo.ea)/(constants.R/constants.Mair*1000)/(meteo.Ta+273.15);
else
    rhoa = constants.rhoa;
end
L           = -rhoa*constants.cp*meteo.ustar.^3.* ...
            (meteo.Ta+273.15)./(constants.kappa*constants.g*H);           % [1]
     
L(isnan(L)) = -1E6;
