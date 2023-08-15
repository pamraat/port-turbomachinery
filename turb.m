function [teta_t, tPr, maxhubr, maxtipr, M, minTt] = ...
    turb(mdotc, Nc, Tt4, Pt4, geom, Mdat)

% This function takes following input:
%               1. mdotc: Corrected mass flow rate (kg/s)
%               2. Nc: Corrected rotational speed (RPM)
%               3. Tt4: Total temperature at the inlet of Turbine (K)
%               4. Pt4: Total Pressure at the inlet of Turbine (Pa)
%               5. geom: Structure datatype containing geometric parameters
%               6. Mdat: Table of corrected massflux as a function of Mach
%               no.
%
%
% This function returns following outputs
%               1. teta_t: Total Turbine efficiency
%               2. tPr: Total Pressure ratio across the Turbine
%               3. maxhubr: Exit hub radius (m)
%               4. maxtipr: Exit tip radius (m)
%               5. M: Exit Mach no.
%               6. minTt: Exit total temperature (K)
%
%   
% This code returns NaN for impossible operating conditions like turbine
% inlet choke.
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
gammat = 1.4;
cp = gammat*R/(gammat - 1);

%% Compressor geometric parameter extraction
ns = geom.turb.ns;                                                          % Number of turbine stages
rt = geom.turb.rt;                                                          % Stage 1 tip radius (m)
rh = geom.turb.rh;                                                          % Stage 1 hub radius (m)
alpha4 = geom.turb.alpha4;                                                  % Angle alpha_4 (degrees) according to Sforza Fig. 7.41
c = 30e-3;                                                                  % Chord length of blade airfoil (m)
sigma = 1;                                                                  % Blade solidity

%% Extraction of relevant parameters
theta = Tt4/Tref;                                                           % Corrected temperature at turbine inlet
delta = Pt4/Pref;                                                           % Corrected pressure at compressor inlet
alpha4 = alpha4/180*pi;
s = c/sigma;                                                                % Space between two blades at mean radius (m)
A4(1) = pi*(rt^2 - rh^2);                                                   % Inlet area of compressor (m^2)
rm = (rt + rh)/2;                                                           % Mean radius
theta_s = s/rm;                                                             % Angle subtended by 2 blades (rad)
nb = round(2*pi/theta_s, 0);                                                % Number of blades
Sb(1) = c*(rt - rh);                                                        % Blade surface area (m^2)

%% Deriving actual quantities from corrected quantities
N = Nc*theta^0.5;
mdot = mdotc*delta/theta^0.5;
U = N*2*pi*rm/60;

%% Definitions of Initial Conditions
Pt(1) = Pt4;
Tt(1) = Tt4;
mdotcorr(1) = mdotc;
Ncorr(1) = Nc;
hubr(1) = rh;
tipr(1) = rt;
M4(1) = interp1(Mdat.Mdat(:,2), Mdat.Mdat(:,1), mdotc/A4(1));
exitflag = 1;
if(1 - M4(1) < 0.001 || isnan(M4(1)))
    exitflag = -1;
end
T4(1) = Tt(1)/(1 + (gammat - 1)/2*M4(1)^2);
P4(1) = Pt(1)/(1 + (gammat - 1)/2*M4(1)^2)^(gammat/(gammat - 1));
V4 = M4(1)*sqrt(gammat*R*T4(1));
beta4 = atan(sin(alpha4)/(U/V4 - cos(alpha4)));
phi = V4*sin(alpha4)/U;

cl = 1.2;
eps = 0.04;
eta_t = 2*phi*(1 - 2*phi*eps)/(eps + 2*phi);

tPr(1) = 1;

% This solver is based on section 7.8 of Sforza, 2nd ed. The work equations
% 7.84 and 7.85 on pg. 364 are the equations of interest. Turbine/compressor
% work is substituted with expression in terms of pressure ratio. The lift
% model on pg. 369, equation 7.97 is followed here. The geometric
% parameters are defined based on fig 7.41 on pg. 363 for turbine blade.
% The blade solidity is at mean radius, used to estimate the number of
% blades in the turbine.

if(exitflag == 1)
    for i=1:ns
        Kt(i) = pi*tipr(i)*nb*Sb(i)*R*Tref/(60*A4(i)^2*Pref)*(T4(i)/Tt(i))*(Pt(i)/P4(i));
    
        Pr(i) = (1 - Kt(i)/(eta_t*cp*Tref)*mdotcorr(i)*Ncorr(i)*cl*(1 - eps*cot(beta4))/sin(beta4))^(gammat/(gammat - 1));
        if(Pr(i) > 1)
            tPr = 1;
            break;
        end
        
        Pt(i + 1) = Pr(i)*Pt(i);
        Tt(i + 1) = Tt(i)*(1 - eta_t*(1 - Pr(i)^(1 - 1/gammat)));
        mdotcorr(i + 1) = mdot*sqrt(Tt(i + 1)/Tref)/(Pt(i + 1)/Pref);
        Ncorr(i + 1) = N/sqrt(Tt(i + 1)/Tref);
        M4(i + 1) = 1/sqrt(gammat*R*Tt(i + 1)/V4^2 - (gammat - 1)/2);
        T4(i + 1) = Tt(i + 1)/(1 + (gammat - 1)/2*M4(i + 1)^2);
        P4(i + 1) = Pt(i + 1)/(1 + (gammat - 1)/2*M4(i + 1)^2)^(gammat/(gammat - 1));
        A4(i + 1) = mdot/V4*sqrt(R*T4(i + 1)/P4(i + 1));
        hubr(i + 1) = rm - A4(i + 1)/(4*pi*rm);
        tipr(i + 1) = rm + A4(i + 1)/(4*pi*rm);
        Sb(i + 1) = c*(tipr(i + 1) - hubr(i + 1));
        tPr = Pr(i)*tPr;
    end
else
    tPr = 1;
end

% Return NaN if Turbine Pressure ratio is greater than 1, else return
% determined values.
if(tPr == 1 || tPr < 0.01)
    tPr = NaN;
    teta_t = NaN;
    maxtipr = NaN;
    maxhubr = NaN;
    minTt = NaN;
    M = NaN;
else
    teta_t = (1 - Tt(end)/Tt(1))/(1 - (tPr)^((gammat - 1)/gammat));
    maxtipr = tipr(end);
    maxhubr = hubr(end);
    minTt = Tt(end);
    M = M4(end);
end
