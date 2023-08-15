function Wdotdiff = wfunc(mdot, N, Tt2, Pt2, Tt4, geom, Mdat)

% Read the description of wfunc in problem statement.
% This function takes following input:
%               1. mdot: Mass flow rate (kg/s)
%               2. N: Rotational speed (RPM)
%               3. Tt2: Total temperature at the inlet of compressor (K)
%               4. Pt2: Total Pressure at the inlet of compressor (Pa)
%               5. geom: Variable of structure datatype with geometric
%               information about compressor and turbine
%               6. Mdat: Table of corrected massflux as a function of
%               Mach no. for gamma = 1.4
%
%
% This function returns following outputs
%               1. Wdiff: Absolute value of difference between turbine and compressor work


Pref = 101325;
Tref = 288.15;
Ru = 8.314;
M_air = 28.96;
R = Ru*1000/M_air;
gamma = 1.4;
cp = gamma*R/(gamma - 1);

%% Compressor
theta2 = Tt2/Tref;
delta2 = Pt2/Pref;
mdotc2 = mdot/delta2*(theta2)^0.5;
Nc2 = N/(theta2)^0.5;

[eta_c, tPr_c, ~, ~, ~, ~] = compr(mdotc2, Nc2, Tt2, Pt2, geom, Mdat);
Pt3 = tPr_c*Pt2;

%% Combustor
Pt4 = Pt3;

%% Turbine
theta4 = Tt4/Tref;
delta4 = Pt4/Pref;
mdotc4 = mdot/delta4*(theta4)^0.5;
Nc4 = N/(theta4)^0.5;

[eta_t, tPr_t, ~, ~, ~, ~] = turb(mdotc4, Nc4, Tt4, Pt4, geom, Mdat);

%% Work Difference
Wdotdiff = abs(mdot*cp*Tt2/eta_c*(tPr_c^(1 - 1/gamma) - 1) - ...
    mdot*eta_t*cp*Tt4*(1 - (tPr_t)^(1 - 1/gamma)));
end