function [r_planet,v_planet,ephem] = meeus2( JD )
% This function will take an input Julian Date and output the ephemeris
% JD is assumed UTI
% calculate T_TDB, equ. 3-42 pg 184
T_TDB = (JD-2451545)/(36525);
ephem ={'Planet';'Semimajor-axis(a)';'Eccentricity(e)';'Inclination(i)';'Lon. of Ascending Node(Omega)';'Longitude of Perihelion(w_true)';'Mean Longitude(lambda)';'Mean Anomaly(M)';'Argument of Perihelion(w)'}; 
% Mercury- Appendix D pg 1046
AU = 149597870.691; %conversion AU to km
a_mer = 0.387098310*AU; %(km)
e_mer = 0.20563175+0.000020406*T_TDB-0.0000000284*T_TDB^2-0.00000000017*T_TDB^3; % no unit
i_mer = 7.004986-0.0059516*T_TDB+0.00000081*T_TDB^2+0.000000041*T_TDB^3; % (deg)
Omega_mer = 48.330893-0.125422*T_TDB-0.00008833*T_TDB^2-0.000000196*T_TDB^3; %(deg)
w_true_mer = 77.456119+0.1588643*T_TDB-0.00001343*T_TDB^2+ 0.000000039*T_TDB^3; %(deg)
lam_M_mer = 252.250906+149472.6746358*T_TDB-0.00000535*T_TDB^2+0.000000002*T_TDB^3; %(deg)
% bound lam_M [0-360]
aa = floor(lam_M_mer/360);
bb = aa*360;
lam_M_mer = lam_M_mer-bb;

M_mer = lam_M_mer-w_true_mer; %(deg)
%bound M [0-360]
aa = floor(M_mer/360);
bb = aa*360;
M_mer = M_mer-bb;

w_mer = w_true_mer-Omega_mer; %(deg)
%bound w [0-360]
aa = floor(w_mer/360);
bb = aa*360;
w_mer = w_mer-bb;

[nu_mer,~,~,~] = Orbit_Anomaly(a_mer,e_mer,[],[],M_mer,[]);
ephem(:,2) = {'Mercury';a_mer;e_mer;i_mer;Omega_mer;w_true_mer;lam_M_mer;M_mer;w_mer};
[r_mer,v_mer] = helioCOE2RV(a_mer,i_mer,Omega_mer,w_mer,nu_mer,e_mer);

% Venus - Appendis D pg 1046
a_ven = 0.723329820*AU; %(km)
e_ven = 0.00677188-0.000047766*T_TDB+0.0000000975*T_TDB^2-0.00000000044*T_TDB^3; % no unit
i_ven = 3.394662-0.0008568*T_TDB-0.00003244*T_TDB^2+0.000000010*T_TDB^3; % (deg)
Omega_ven = 76.679920-0.2780080*T_TDB-0.00014256*T_TDB^2-0.000000198*T_TDB^3; %(deg)
w_true_ven = 131.563707+0.0048646*T_TDB-0.00138232*T_TDB^2-0.000005332*T_TDB^3; %(deg)
lam_M_ven = 181.979801+58517.8156760*T_TDB+0.00000165*T_TDB^2-0.000000002*T_TDB^3; %(deg)

% bound lam_M [0-360]
aa = floor(lam_M_ven/360);
bb = aa*360;
lam_M_ven = lam_M_ven-bb;

M_ven = lam_M_ven-w_true_ven; %(deg)
%bound M [0-360]
aa = floor(M_ven/360);
bb = aa*360;
M_ven = M_ven-bb;

w_ven = w_true_ven-Omega_ven; %(deg)
%bound w [0-360]
aa = floor(w_ven/360);
bb = aa*360;
w_ven = w_ven-bb;

[nu_ven,~,~,~] = Orbit_Anomaly(a_ven,e_ven,[],[],M_ven,[]);
ephem(:,3) = {'Venus';a_ven;e_ven;i_ven;Omega_ven;w_true_ven;lam_M_ven;M_ven;w_ven};
[r_ven,v_ven] = helioCOE2RV(a_ven,i_ven,Omega_ven,w_ven,nu_ven,e_ven);

% Earth - Appen. D pg 1046
a_ear = 1.000001018*AU; %(km)
e_ear = 0.01670862-0.000042037*T_TDB-0.0000001236*T_TDB^2+0.00000000004*T_TDB^3; % no unit
i_ear = 0.0000000+0.0130546*T_TDB-0.00000931*T_TDB^2-0.000000034*T_TDB^3; % (deg)
Omega_ear = 174.873174-0.2410908*T_TDB+0.00004067*T_TDB^2-0.000001327*T_TDB^3; %(deg)
w_true_ear = 102.937348+0.3225557*T_TDB+0.00015026*T_TDB^2+ 0.000000478*T_TDB^3; %(deg)
lam_M_ear = 100.466449+35999.3728519*T_TDB-0.00000568*T_TDB^2+0.000000000*T_TDB^3; %(deg)

% bound lam_M [0-360]
aa = floor(lam_M_ear/360);
bb = aa*360;
lam_M_ear = lam_M_ear-bb;

M_ear = lam_M_ear-w_true_ear; %(deg)
%bound M [0-360]
aa = floor(M_ear/360);
bb = aa*360;
M_ear = M_ear-bb;

w_ear = w_true_ear-Omega_ear; %(deg)
%bound w [0-360]
aa = floor(w_ear/360);
bb = aa*360;
w_ear = w_ear-bb;

[nu_ear,~,~,~] = Orbit_Anomaly(a_ear,e_ear,[],[],M_ear,[]);
ephem(:,4) = {'Earth';a_ear;e_ear;i_ear;Omega_ear;w_true_ear;lam_M_ear;M_ear;w_ear};
[r_ear,v_ear] = helioCOE2RV(a_ear,i_ear,Omega_ear,w_ear,nu_ear,e_ear);

% Mars - Appen. D pg 1047
a_mar = 1.523679342*AU; %(km)
e_mar = 0.09340062+0.000090483*T_TDB-0.0000000806*T_TDB^2-0.00000000035*T_TDB^3; % no unit
i_mar = 1.849726-0.0081479*T_TDB-0.00002255*T_TDB^2-0.000000027*T_TDB^3; % (deg)
Omega_mar = 49.558093-0.2949846*T_TDB-0.00063993*T_TDB^2-0.000002143*T_TDB^3; %(deg)
w_true_mar = 336.060234+0.4438898*T_TDB-0.00017321*T_TDB^2+0.000000300*T_TDB^3; %(deg)
lam_M_mar = 355.433275+19140.2993313*T_TDB+0.00000261*T_TDB^2-0.000000003*T_TDB^3; %(deg)

% bound lam_M [0-360]
aa = floor(lam_M_mar/360);
bb = aa*360;
lam_M_mar = lam_M_mar-bb;

M_mar = lam_M_mar-w_true_mar; %(deg)
%bound M [0-360]
aa = floor(M_mar/360);
bb = aa*360;
M_mar = M_mar-bb;

w_mar = w_true_mar-Omega_mar; %(deg)
%bound w [0-360]
aa = floor(w_mar/360);
bb = aa*360;
w_mar = w_mar-bb;

[nu_mar,~,~,~] = Orbit_Anomaly(a_mar,e_mar,[],[],M_mar,[]);
ephem(:,5) = {'Mars';a_mar;e_mar;i_mar;Omega_mar;w_true_mar;lam_M_mar;M_mar;w_mar};
[r_mar,v_mar] = helioCOE2RV(a_mar,i_mar,Omega_mar,w_mar,nu_mar,e_mar);

% Jupiter - Appen D pg 1047
a_jup = (5.202603191+0.0000001913*T_TDB)*AU; %(km)
e_jup = 0.04849485+0.000163244*T_TDB-0.0000004719*T_TDB^2-0.00000000197*T_TDB^3; % no unit
i_jup = 1.303270-0.0019872*T_TDB+0.00003318*T_TDB^2+0.000000092*T_TDB^3; % (deg)
Omega_jup = 100.464441+0.1766828*T_TDB+0.00090387*T_TDB^2-0.000007032*T_TDB^3; %(deg)
w_true_jup = 14.331309+0.2155525*T_TDB+0.00072252*T_TDB^2-0.000004590*T_TDB^3; %(deg)
lam_M_jup = 34.351484+3034.9056746*T_TDB-0.00008501*T_TDB^2+0.000000004*T_TDB^3; %(deg)

% bound lam_M [0-360]
aa = floor(lam_M_jup/360);
bb = aa*360;
lam_M_jup = lam_M_jup-bb;

M_jup = lam_M_jup-w_true_jup; %(deg)
%bound M [0-360]
aa = floor(M_jup/360);
bb = aa*360;
M_jup = M_jup-bb;

w_jup = w_true_jup-Omega_jup; %(deg)
%bound w [0-360]
aa = floor(w_jup/360);
bb = aa*360;
w_jup = w_jup-bb;

[nu_jup,~,~,~] = Orbit_Anomaly(a_jup,e_jup,[],[],M_jup,[]);
ephem(:,6) = {'Jupiter';a_jup;e_jup;i_jup;Omega_jup;w_true_jup;lam_M_jup;M_jup;w_jup};
[r_jup,v_jup] = helioCOE2RV(a_jup,i_jup,Omega_jup,w_jup,nu_jup,e_jup);

% Saturn - Appen D. pg 1047
a_sat = (9.554909596-.0000021389*T_TDB)*AU; %(km)
e_sat = 0.05550862-0.000346818*T_TDB-0.0000006456*T_TDB^2+0.00000000338*T_TDB^3; % no unit
i_sat = 2.488878+0.0025515*T_TDB-0.00004903*T_TDB^2+0.000000018*T_TDB^3; % (deg)
Omega_sat = 113.665524-0.2566649*T_TDB-0.00018345*T_TDB^2+0.000000357*T_TDB^3; %(deg)
w_true_sat = 93.056787+0.5665496*T_TDB+0.00052809*T_TDB^2+0.000004882*T_TDB^3; %(deg)
lam_M_sat = 50.077471+1222.1137943*T_TDB+0.00021004*T_TDB^2-0.000000019*T_TDB^3; %(deg)

% bound lam_M [0-360]
aa = floor(lam_M_sat/360);
bb = aa*360;
lam_M_sat = lam_M_sat-bb;

M_sat = lam_M_sat-w_true_sat; %(deg)
%bound M [0-360]
aa = floor(M_sat/360);
bb = aa*360;
M_sat = M_sat-bb;

w_sat = w_true_sat-Omega_sat; %(deg)
%bound w [0-360]
aa = floor(w_sat/360);
bb = aa*360;
w_sat = w_sat-bb;

[nu_sat,~,~,~] = Orbit_Anomaly(a_sat,e_sat,[],[],M_sat,[]);
ephem(:,7) = {'Saturn';a_sat;e_sat;i_sat;Omega_sat;w_true_sat;lam_M_sat;M_sat;w_sat};
[r_sat,v_sat] = helioCOE2RV(a_sat,i_sat,Omega_sat,w_sat,nu_sat,e_sat);

% Uranus - Appen D. pg 1047-1048
a_ura = (19.218446062-0.0000000372*T_TDB+0.000000000098*T_TDB^2)*AU; %(km)
e_ura = 0.04629590-0.000027337*T_TDB+0.0000000790*T_TDB^2+0.00000000025*T_TDB^3; % no unit
i_ura = 0.773196-0.0016869*T_TDB+0.00000349*T_TDB^2+0.000000016*T_TDB^3; % (deg)
Omega_ura = 74.005947+0.0741461*T_TDB+0.00040540*T_TDB^2+0.000000104*T_TDB^3; %(deg)
w_true_ura = 173.005159+0.0893206*T_TDB-0.00009470*T_TDB^2+0.000000413*T_TDB^3; %(deg)
lam_M_ura = 314.055005+428.4669983*T_TDB-0.00000486*T_TDB^2+0.000000006*T_TDB^3; %(deg)

% bound lam_M [0-360]
aa = floor(lam_M_ura/360);
bb = aa*360;
lam_M_ura = lam_M_ura-bb;

M_ura = lam_M_ura-w_true_ura; %(deg)
%bound M [0-360]
aa = floor(M_ura/360);
bb = aa*360;
M_ura = M_ura-bb;

w_ura = w_true_ura-Omega_ura; %(deg)
%bound w [0-360]
aa = floor(w_ura/360);
bb = aa*360;
w_ura = w_ura-bb;

[nu_ura,~,~,~] = Orbit_Anomaly(a_ura,e_ura,[],[],M_ura,[]);
ephem(:,8) = {'Uranus';a_ura;e_ura;i_ura;Omega_ura;w_true_ura;lam_M_ura;M_ura;w_ura};
[r_ura,v_ura] = helioCOE2RV(a_ura,i_ura,Omega_ura,w_ura,nu_ura,e_ura);

% Neptune - Appen D pg 1048
a_nep = (30.110386869-0.0000001663*T_TDB+0.00000000069*T_TDB^2)*AU; %(km)
e_nep = 0.00898809+0.000006408*T_TDB-0.0000000008*T_TDB^2; % no unit
i_nep = 1.769952+0.0002257*T_TDB+0.00000023*T_TDB^2-0.000000000*T_TDB^3; % (deg)
Omega_nep = 131.784057-0.0061651*T_TDB-0.00000219*T_TDB^2-0.000000078*T_TDB^3; %(deg)
w_true_nep = 48.123691+0.0291587*T_TDB+0.00007051*T_TDB^2-0.000000000*T_TDB^3; %(deg)
lam_M_nep = 304.348665+218.4862002*T_TDB+0.00000059*T_TDB^2-0.000000002*T_TDB^3; %(deg)

% bound lam_M [0-360]
aa = floor(lam_M_nep/360);
bb = aa*360;
lam_M_nep = lam_M_nep-bb;

M_nep = lam_M_nep-w_true_nep; %(deg)
%bound M [0-360]
aa = floor(M_nep/360);
bb = aa*360;
M_nep = M_nep-bb;

w_nep = w_true_nep-Omega_nep; %(deg)
%bound w [0-360]
aa = floor(w_nep/360);
bb = aa*360;
w_nep = w_nep-bb;

[nu_nep,~,~,~] = Orbit_Anomaly(a_nep,e_nep,[],[],M_nep,[]);
ephem(:,9) = {'Neptune';a_nep;e_nep;i_nep;Omega_nep;w_true_nep;lam_M_nep;M_nep;w_nep};
[r_nep,v_nep] = helioCOE2RV(a_nep,i_nep,Omega_nep,w_nep,nu_nep,e_nep);

% Final planet parameters matrix
r_planet = [r_mer,r_ven,r_ear,r_mar,r_jup,r_sat,r_ura,r_nep];
v_planet = [v_mer,v_ven,v_ear,v_mar,v_jup,v_sat,v_ura,v_nep];
end