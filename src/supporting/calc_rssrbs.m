function [rss,rbs] = calc_rssrbs(soil,LAI,Ts_old1,options)

if Ts_old1 == -999
    Ts = 20 + 273.15;
else
    Ts = Ts_old1 + 273.15*(Ts_old1<100);
end
SMC = soil.SMC;

if options.calc_rss_rbs == 1
    % SCOPE's built-in function
    rss        = 11.2*exp(42*(0.22-SMC));
    rbs        = soil.rbs*LAI/3.3;
elseif options.calc_rss_rbs == 2
    % Dry Surface Layer (DSL) model used in CLM5 (Swenson & Lawrence, JGR 2014)
    theta_s = soil.porosity;                          % [m3 m-3]
    psi_s = soil.psi_sat;                             % saturated (air entry) hydraulic head [J kg-1] or [mm]
    b = soil.b;                                       % Campbell's "b"
    theta_i = 0.8*theta_s;                            % SMC value at which the DSL initiates [m3 m-3]
    theta_a = theta_s * (psi_s*10^-7)^(1/b);          % air dry SMC [m3 m-3]
    phi_a = theta_s - theta_a;                        % air filled pore space [m3 m-3]
    Dv = 2.12 *10^-5 * (Ts/273.15)^1.75;              % molecular diffusivity of water vapour [m2 s-2]
    tau = phi_a^2 * (phi_a/theta_s)^(3/soil.b);       % gas diffusivity factor in soil
    if SMC < theta_i
        Dmax = 0.015;                                 % maximum depth of DSL [m]
        DSL = Dmax * (theta_i-SMC)/(theta_i-theta_a); % [m]
    else
        DSL = 0;
    end
    rss = DSL/Dv/tau;
    rbs = soil.rbs;
else
    rss = NaN;
    rbs = NaN;
end