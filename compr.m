function [teta_c, tPr, minhubr, mintipr, M, maxTt] =...
    compr(mdotc, Nc, Tt2, Pt2, geom, Mdat)

% This function takes following input:
%               1. mdotc: Corrected mass flow rate (kg/s)
%               2. Nc: Corrected rotational speed (RPM)
%               3. Tt2: Total temperature at the inlet of compressor (K)
%               4. Pt2: Total Pressure at the inlet of compressor (Pa)
%               5. geom: Structure datatype containing geometric parameters
%               of the compressor
%               6. Mdat: Table of corrected massflux as a function of Mach
%               no.
%
%
% This function returns following outputs
%               1. teta_c: Total compressor efficiency
%               2. tPr: Total Pressure ratio across the comrpessor
%               3. minhubr: Exit hub radius (m)
%               4. mintipr: Exit tip radius (m)
%               7. M: Exit Mach no.
%               6. maxTt: Exit total temperature (K)
%
%   
% This code returns NaN for impossible operating conditions like compressor
% inlet choke, windmilling, or compressor stalling.
%
%
% This code is developed by Pamraat Parmar and Elizabeth Fisher for
% educational purpose to aid in instruction of MAE 4510/5510 of Cornell
% University. This code open to use for academic purposes.

%% Thermodynamic Properties
Pref = 101325;
Tref = 288.15;
Ru = 8.314;
M_air = 28.96;
R = Ru*1000/M_air;
gamma = 1.4;
cp = gamma*R/(gamma - 1);

%% Compressor geometric parameter extraction
ns = geom.comp.ns;                                                          % Number of compressor stages
rt = geom.comp.rt;                                                          % Stage 1 tip radius (m)
rh = geom.comp.rh;                                                          % Stage 1 hub radius (m)
Gam = 47;                                                                   % Blade stacking angle Gamma (degrees) according to Sforza Fig. 7.41
alpha2 = geom.comp.alpha2;                                                  % Angle alpha_2 (degrees) according to Sforza Fig. 7.41
c = 20e-3;                                                                  % Chord length of blade airfoil (m)
sigma = 0.5;                                                                % Blade solidity

%% Extraction of relevant parameters
theta = Tt2/Tref;                                                           % Corrected temperature at compressor inlet
delta = Pt2/Pref;                                                           % Corrected pressure at compressor inlet
Gam = Gam/180*pi;
alpha2 = alpha2/180*pi;
s = c/sigma;                                                                % Space between two blades at mean radius (m)
A2(1) = pi*(rt^2 - rh^2);                                                   % Inlet area of compressor (m^2)
rm = (rt + rh)/2;                                                           % Mean radius
theta_s = s/rm;                                                             % Angle subtended by 2 blades (rad)
nb = round(2*pi/theta_s, 0);                                                % Number of blades
Sb(1) = c*(rt - rh);                                                        % Blade surface area (m^2)

%% Deriving actual quantities from corrected quantities
N = Nc*theta^0.5;
mdot = mdotc*delta/theta^0.5;
U = N*2*pi*rm/60;

%% Definitions of Initial Conditions
Pt(1) = Pt2;
Tt(1) = Tt2;
mdotcorr(1) = mdotc;
Ncorr(1) = Nc;
hubr(1) = rh;
tipr(1) = rt;
M2(1) = interp1(Mdat.Mdat(:,2), Mdat.Mdat(:,1), mdotc/A2(1));               % Determining mach no. at 2
exitflag = 1;
if(1 - M2(1) < 0.001)
    exitflag = -1;                                                          % Exit with NaN if inlet of compressor is choked
end
T2(1) = Tt(1)/(1 + (gamma - 1)/2*M2(1)^2);
P2(1) = Pt(1)/(1 + (gamma - 1)/2*M2(1)^2)^(gamma/(gamma - 1));
V2 = M2(1)*sqrt(gamma*R*T2(1));                                             % Determine V2 from M2
beta2 = atan(sin(alpha2)/(U/V2 - cos(alpha2)));                             % Determine beta2 from knowledge of N and v2.
phi = V2*sin(alpha2)/U;                                                     % Flow rate coefficient

[cl, eps] = lift(beta2, Gam);    
eta_c = 2*phi*(1 - 2*phi*eps)/(eps + 2*phi);                                % Equation 8.75 of Sforza 2nd ed.

tPr(1) = 1;

% This solver is based on section 7.8 of Sforza, 2nd ed. The work equations
% 7.84 and 7.85 on pg. 364 are the equations of interest. Turbine/compressor
% work is substituted with expression in terms of pressure ratio. The lift
% model on pg. 369, equation 7.97 is followed here. The geometric
% parameters are defined based on fig 7.41 on pg. 363 for compressor blade.
% The blade solidity is at mean radius, used to estimate the number of
% blades in the compressor.


if(exitflag == 1)
    for i=1:ns
        Kc(i) = pi*tipr(i)*nb*Sb(i)*R*Tref/(60*A2(i)^2*Pref)*(T2(i)/Tt(i))*(Pt(i)/P2(i));
    
        Pr(i) = (1 + eta_c/(cp*Tref)*Kc(i)*mdotcorr(i)*Ncorr(i)*cl*(1 + eps*cot(beta2))/sin(beta2))^(gamma/(gamma - 1));
        if(Pr(i) < 1)
            tPr = 1;
            break;
        end
        
        Pt(i + 1) = Pr(i)*Pt(i);
        Tt(i + 1) = Tt(i)*(1 + 1/eta_c*(Pr(i)^(1 - 1/gamma) - 1));
        mdotcorr(i + 1) = mdot*sqrt(Tt(i + 1)/Tref)/(Pt(i + 1)/Pref);
        Ncorr(i + 1) = N/sqrt(Tt(i + 1)/Tref);
        M2(i + 1) = 1/sqrt(gamma*R*Tt(i + 1)/V2^2 - (gamma - 1)/2);
        T2(i + 1) = Tt(i + 1)/(1 + (gamma - 1)/2*M2(i + 1)^2);
        P2(i + 1) = Pt(i + 1)/(1 + (gamma - 1)/2*M2(i + 1)^2)^(gamma/(gamma - 1));
        A2(i + 1) = mdot/V2*sqrt(R*T2(i + 1)/P2(i + 1));
        hubr(i + 1) = rm - A2(i + 1)/(4*pi*rm);
        tipr(i + 1) = rm + A2(i + 1)/(4*pi*rm);
        Sb(i + 1) = c*(tipr(i + 1) - hubr(i + 1));
        tPr = Pr(i)*tPr;
    end
else
    tPr = 1;
end

% Return NaN if compressor pressure ratio is less than 1, else return
% determined values.
if(tPr == 1)
    tPr = NaN;
    teta_c = NaN;
    mintipr = NaN;
    minhubr = NaN;
    maxTt = NaN;
    M = NaN;
else
    teta_c = ((tPr)^((gamma - 1)/gamma) - 1)/(Tt(end)/Tt(1) - 1);
    mintipr = tipr(end - 1);
    minhubr = hubr(end - 1);
    maxTt = Tt(end);
    M = M2(end);
end

