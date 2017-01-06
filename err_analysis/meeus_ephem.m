
function [ephem] = meeus_ephem(JDE)
%Name: meeus_ephem.m
%Author: Michael Von Hendy, edited by Martin Heaney
%Date: 1/28/15
%Purpose: Compute planetary ephemerides from meeus approximation.
%
%Inputs: JDE = julian date epoch
%
%Outputs: ephem = planetary ephemeris data object
%
%-----Output Format------
%   ephem.planet_L = mean longitude of planet [deg]
%   ephem.planet_a = semimajor axis of planet [km], (yes, output is km)
%   ephem.planet_e = eccentircity
%   ephem.planet_i = inclination [rad]
%   ephem.planet_Omega = Longitude of ascending node [rad]
%   ephem.planet_LOP = longitude of perihelion [deg]
%   ephem.planet_w = argument of perihelion [rad]
%   ephem.planet_M = mean anomaly [deg]
%   ephem.planet_nu = true anomaly [rad]
%------------------------


%Conversion:
AU2km = 149597870.691; %km

%% Compute T:

T = (JDE - 2451545.0)/36525; %Julian centuries

%% Mercury
L.a0 = 252.250906;
L.a1 = 149472.6746358;
L.a2 = -0.00000535;
L.a3 = 0.000000002;

a.a0 = 0.387098310;

e.a0 = 0.20563175;
e.a1 = 0.000020406;
e.a2 = -0.0000000284;
e.a3 = -0.00000000017;

i.a0 = 7.004986;
i.a1 = -0.0059516;
i.a2 = 0.00000081;
i.a3 = 0.000000041;

Omega.a0 = 48.330893;
Omega.a1 = -0.125422;
Omega.a2 = -0.00008833;
Omega.a3 = -0.000000196;

LOP.a0 = 77.456119;
LOP.a1 = 0.1588643;
LOP.a2 = -0.00001343;
LOP.a3 = 0.000000039;

ephem.mercury_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.mercury_a = a.a0 * AU2km; %km
ephem.mercury_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.mercury_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.mercury_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.mercury_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.mercury_w = deg2rad(ephem.mercury_LOP) - ephem.mercury_Omega; %rad
ephem.mercury_M = mod(rad2deg(deg2rad(ephem.mercury_L)-deg2rad(ephem.mercury_LOP)),360); %deg
Ccen = (2*ephem.mercury_e-ephem.mercury_e^3/4+5/96*ephem.mercury_e^5)*sind(ephem.mercury_M)+...
        (5/4*ephem.mercury_e^2-11/24*ephem.mercury_e^4)*sind(2*ephem.mercury_M) +...
         (13/12*ephem.mercury_e^3 - 43/64*ephem.mercury_e^5)*sind(3*ephem.mercury_M) + ...
          103/96*ephem.mercury_e^4*sind(4*ephem.mercury_M)+1097/960*ephem.mercury_e^5*sind(5*ephem.mercury_M);
ephem.mercury_nu = deg2rad(mod(rad2deg(deg2rad(ephem.mercury_M)+Ccen),360)); %rad

%% Venus
L.a0 = 181.979801;
L.a1 = 58517.8156760;
L.a2 = 0.00000165;
L.a3 = -0.000000002;

a.a0 = 0.72332982;

e.a0 = 0.00677188;
e.a1 = -0.000047766;
e.a2 = 0.0000000975;
e.a3 = 0.00000000044;

i.a0 = 3.394662;
i.a1 = -0.0008568;
i.a2 = -0.00003244;
i.a3 = 0.000000010;

Omega.a0 = 76.679920;
Omega.a1 = -0.2780080;
Omega.a2 = -0.00014256;
Omega.a3 = -0.000000198;

LOP.a0 = 131.563707;
LOP.a1 = 0.0048646;
LOP.a2 = -0.00138232;
LOP.a3 = -0.000005332;

ephem.venus_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.venus_a = a.a0 * AU2km; %km
ephem.venus_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.venus_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.venus_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.venus_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.venus_w = deg2rad(ephem.venus_LOP) - ephem.venus_Omega; %rad
ephem.venus_M = mod(rad2deg(deg2rad(ephem.venus_L)-deg2rad(ephem.venus_LOP)),360); %deg
Ccen = (2*ephem.venus_e-ephem.venus_e^3/4+5/96*ephem.venus_e^5)*sind(ephem.venus_M)+...
        (5/4*ephem.venus_e^2-11/24*ephem.venus_e^4)*sind(2*ephem.venus_M) +...
         (13/12*ephem.venus_e^3 - 43/64*ephem.venus_e^5)*sind(3*ephem.venus_M) + ...
          103/96*ephem.venus_e^4*sind(4*ephem.venus_M)+1097/960*ephem.venus_e^5*sind(5*ephem.venus_M);
ephem.venus_nu = deg2rad(mod(rad2deg(deg2rad(ephem.venus_M)+Ccen),360)); %rad

%% Earth
L.a0 = 100.466449;
L.a1 = 35999.3728519;
L.a2 = -0.00000568;
L.a3 = 0;

a.a0 = 1.000001018;

e.a0 = 0.01670862;
e.a1 = -0.000042037;
e.a2 = -0.0000001236;
e.a3 = 0.00000000004;

i.a0 = 0;
i.a1 = 0.0130546;
i.a2 = -0.00000931;
i.a3 = -0.000000034;

Omega.a0 = 174.873174;
Omega.a1 = -0.2410908;
Omega.a2 = 0.00004067;
Omega.a3 = -0.000001327;

LOP.a0 = 102.937348; 
LOP.a1 = 0.3225557;
LOP.a2 = 0.00015026;
LOP.a3 = 0.000000478;

ephem.earth_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.earth_a = a.a0 * AU2km; %km
ephem.earth_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.earth_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.earth_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.earth_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.earth_w = deg2rad(ephem.earth_LOP) - ephem.earth_Omega; %deg
ephem.earth_M = mod(rad2deg(deg2rad(ephem.earth_L)-deg2rad(ephem.earth_LOP)),360); %deg
Ccen = (2*ephem.earth_e-ephem.earth_e^3/4+5/96*ephem.earth_e^5)*sind(ephem.earth_M)+...
        (5/4*ephem.earth_e^2-11/24*ephem.earth_e^4)*sind(2*ephem.earth_M) +...
         (13/12*ephem.earth_e^3 - 43/64*ephem.earth_e^5)*sind(3*ephem.earth_M) + ...
          103/96*ephem.earth_e^4*sind(4*ephem.earth_M)+1097/960*ephem.earth_e^5*sind(5*ephem.earth_M);
ephem.earth_nu = deg2rad(mod(rad2deg(deg2rad(ephem.earth_M)+Ccen),360)); %rad

%% Mars
L.a0 = 355.433275;
L.a1 = 19140.2993313;
L.a2 = 0.00000261;
L.a3 = -0.000000003;

a.a0 = 1.523679342;

e.a0 = 0.09340062;
e.a1 = 0.000090483;
e.a2 = -0.0000000806;
e.a3 = -0.00000000035;

i.a0 = 1.849726;
i.a1 = -0.0081479;
i.a2 = -0.00002255;
i.a3 = -0.000000027;

Omega.a0 = 49.558093;
Omega.a1 = -0.2949846;
Omega.a2 = -0.00063993;
Omega.a3 = -0.000002143;

LOP.a0 = 336.060234;
LOP.a1 = 0.4438898;
LOP.a2 = -0.00017321;
LOP.a3 = 0.000000300;

ephem.mars_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.mars_a = a.a0 * AU2km; %km
ephem.mars_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.mars_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.mars_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.mars_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.mars_w = deg2rad(ephem.mars_LOP) - ephem.mars_Omega; %deg
ephem.mars_M = mod(rad2deg(deg2rad(ephem.mars_L)-deg2rad(ephem.mars_LOP)),360); %deg
Ccen = (2*ephem.mars_e-ephem.mars_e^3/4+5/96*ephem.mars_e^5)*sind(ephem.mars_M)+...
        (5/4*ephem.mars_e^2-11/24*ephem.mars_e^4)*sind(2*ephem.mars_M) +...
         (13/12*ephem.mars_e^3 - 43/64*ephem.mars_e^5)*sind(3*ephem.mars_M) + ...
          103/96*ephem.mars_e^4*sind(4*ephem.mars_M)+1097/960*ephem.mars_e^5*sind(5*ephem.mars_M);
ephem.mars_nu = deg2rad(mod(rad2deg(deg2rad(ephem.mars_M)+Ccen),360)); %rad

%% Jupiter
L.a0 = 34.351484;
L.a1 = 3034.9056746;
L.a2 = -0.00008501;
L.a3 = 0.000000004;

a.a0 = 5.202603191;
a.a1 = 0.0000001913;

e.a0 = 0.04849485;
e.a1 = 0.000163244;
e.a2 = -0.0000004719;
e.a3 = -0.00000000197;

i.a0 = 1.303270;
i.a1 = -0.0019872;
i.a2 = 0.00003318;
i.a3 = 0.000000092;

Omega.a0 = 100.464441;
Omega.a1 = 0.1766828;
Omega.a2 = 0.00090387;
Omega.a3 = -0.000007032;

LOP.a0 = 14.331309;
LOP.a1 = 0.2155525;
LOP.a2 = 0.00072252;
LOP.a3 = -0.000004590;

ephem.jupiter_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.jupiter_a = (a.a0 +a.a1*T)* AU2km; %km
ephem.jupiter_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.jupiter_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.jupiter_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.jupiter_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.jupiter_w = deg2rad(ephem.jupiter_LOP) - ephem.jupiter_Omega; %deg
ephem.jupiter_M = mod(rad2deg(deg2rad(ephem.jupiter_L)-deg2rad(ephem.jupiter_LOP)),360); %deg
Ccen = (2*ephem.jupiter_e-ephem.jupiter_e^3/4+5/96*ephem.jupiter_e^5)*sind(ephem.jupiter_M)+...
        (5/4*ephem.jupiter_e^2-11/24*ephem.jupiter_e^4)*sind(2*ephem.jupiter_M) +...
         (13/12*ephem.jupiter_e^3 - 43/64*ephem.jupiter_e^5)*sind(3*ephem.jupiter_M) + ...
          103/96*ephem.jupiter_e^4*sind(4*ephem.jupiter_M)+1097/960*ephem.jupiter_e^5*sind(5*ephem.jupiter_M);
ephem.jupiter_nu = deg2rad(mod(rad2deg(deg2rad(ephem.jupiter_M)+Ccen),360)); %rad


%% Saturn
L.a0 = 50.077471;
L.a1 = 1222.1137943;
L.a2 = 0.00021004;
L.a3 = -0.000000019;

a.a0 = 9.554909596;
a.a1 = -0.0000021389;

e.a0 = 0.05550862;
e.a1 = -0.000346818;
e.a2 = -0.0000006456;
e.a3 = 0.00000000338;

i.a0 = 2.488878;
i.a1 = 0.0025515;
i.a2 = -0.00004903;
i.a3 = 0.000000018;

Omega.a0 = 113.665524;
Omega.a1 = -0.2566649;
Omega.a2 = -0.00018345;
Omega.a3 = 0.000000357;

LOP.a0 = 93.056787;
LOP.a1 = 0.5665496;
LOP.a2 = 0.00052809;
LOP.a3 = 0.000004882;

ephem.saturn_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.saturn_a = (a.a0 +a.a1*T)* AU2km; %km
ephem.saturn_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.saturn_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.saturn_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.saturn_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.saturn_w = deg2rad(ephem.saturn_LOP) - ephem.saturn_Omega; %deg
ephem.saturn_M = mod(rad2deg(deg2rad(ephem.saturn_L)-deg2rad(ephem.saturn_LOP)),360); %deg
Ccen = (2*ephem.saturn_e-ephem.saturn_e^3/4+5/96*ephem.saturn_e^5)*sind(ephem.saturn_M)+...
        (5/4*ephem.saturn_e^2-11/24*ephem.saturn_e^4)*sind(2*ephem.saturn_M) +...
         (13/12*ephem.saturn_e^3 - 43/64*ephem.saturn_e^5)*sind(3*ephem.saturn_M) + ...
          103/96*ephem.saturn_e^4*sind(4*ephem.saturn_M)+1097/960*ephem.saturn_e^5*sind(5*ephem.saturn_M);
ephem.saturn_nu = deg2rad(mod(rad2deg(deg2rad(ephem.saturn_M)+Ccen),360)); %rad

%% Uranus

L.a0 = 314.055005;
L.a1 = 428.4669983;
L.a2 = -0.00000486;
L.a3 = 0.000000006;

a.a0 = 19.218446062;
a.a1 = -0.0000000372;
a.a2 = 0.000000000098;

e.a0 = 0.04629590;
e.a1 = -0.000027337;
e.a2 = 0.0000000790;
e.a3 = 0.00000000025;

i.a0 = 0.773196;
i.a1 = -0.0016869;
i.a2 = 0.00000349;
i.a3 = 0.000000016;

Omega.a0 = 74.005947;
Omega.a1 = 0.0741461;
Omega.a2 = 0.00040540;
Omega.a3 = 0.000000104;

LOP.a0 = 173.005159;
LOP.a1 = 0.0893206;
LOP.a2 = -0.00009470;
LOP.a3 = 0.000000413;

ephem.uranus_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.uranus_a = (a.a0 +a.a1*T+a.a2*T^2)* AU2km; %km
ephem.uranus_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.uranus_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.uranus_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.uranus_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.uranus_w = deg2rad(ephem.uranus_LOP) - ephem.uranus_Omega; %deg
ephem.uranus_M = mod(rad2deg(deg2rad(ephem.uranus_L)-deg2rad(ephem.uranus_LOP)),360); %deg
Ccen = (2*ephem.uranus_e-ephem.uranus_e^3/4+5/96*ephem.uranus_e^5)*sind(ephem.uranus_M)+...
        (5/4*ephem.uranus_e^2-11/24*ephem.uranus_e^4)*sind(2*ephem.uranus_M) +...
         (13/12*ephem.uranus_e^3 - 43/64*ephem.uranus_e^5)*sind(3*ephem.uranus_M) + ...
          103/96*ephem.uranus_e^4*sind(4*ephem.uranus_M)+1097/960*ephem.uranus_e^5*sind(5*ephem.uranus_M);
ephem.uranus_nu = deg2rad(mod(rad2deg(deg2rad(ephem.uranus_M)+Ccen),360)); %rad



%% Neptune 

L.a0 = 304.348665;
L.a1 = 218.4862002;
L.a2 = 0.00000059;
L.a3 = -0.000000002;

a.a0 = 30.110386869;
a.a1 = -0.0000001663;
a.a2 = 0.00000000069;

e.a0 = 0.00898809;
e.a1 = 0.000006408;
e.a2 = -0.0000000008;
e.a3 = 0;

i.a0 = 1.769952;
i.a1 = 0.0002257;
i.a2 = 0;
i.a3 = 0;

Omega.a0 = 131.784057;
Omega.a1 = -0.0061651;
Omega.a2 = -0.00000219;
Omega.a3 = -0.000000078;

LOP.a0 = 48.123691;
LOP.a1 = 0.0291587;
LOP.a2 = 0.00007051;
LOP.a3 = 0;

ephem.neptune_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.neptune_a = (a.a0 +a.a1*T+a.a2*T^2)* AU2km; %km
ephem.neptune_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.neptune_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.neptune_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.neptune_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.neptune_w = deg2rad(ephem.neptune_LOP) - ephem.neptune_Omega; %deg
ephem.neptune_M = mod(rad2deg(deg2rad(ephem.neptune_L)-deg2rad(ephem.neptune_LOP)),360); %deg
Ccen = (2*ephem.neptune_e-ephem.neptune_e^3/4+5/96*ephem.neptune_e^5)*sind(ephem.neptune_M)+...
        (5/4*ephem.neptune_e^2-11/24*ephem.neptune_e^4)*sind(2*ephem.neptune_M) +...
         (13/12*ephem.neptune_e^3 - 43/64*ephem.neptune_e^5)*sind(3*ephem.neptune_M) + ...
          103/96*ephem.neptune_e^4*sind(4*ephem.neptune_M)+1097/960*ephem.neptune_e^5*sind(5*ephem.neptune_M);
ephem.neptune_nu = deg2rad(mod(rad2deg(deg2rad(ephem.neptune_M)+Ccen),360)); %rad


%% Pluto
L.a0 = 238.92903833;
L.a1 = 145.20780515;
L.a2 = 0.0;
L.a3 = 0.0;

a.a0 = 39.48211675;
a.a1 = -0.00031596;

e.a0 = 0.24882730;
e.a1 = 0.00005170;
e.a2 = 0;
e.a3 = 0;

i.a0 = 17.14001206;
i.a1 = 0.00004818;
i.a2 = 0;
i.a3 = 0;

Omega.a0 = 110.30393684;
Omega.a1 = -0.01183482;
Omega.a2 = 0;
Omega.a3 = 0;

LOP.a0 = 224.06891629;
LOP.a1 = -0.04062942;
LOP.a2 = 0;
LOP.a3 = 0;

ephem.pluto_L = mod(L.a0 + L.a1*T + L.a2*T^2 + L.a3*T^3,360); %deg
ephem.pluto_a = (a.a0 +a.a1*T)* AU2km; %km
ephem.pluto_e = e.a0 + e.a1*T + e.a2*T^2 + e.a3*T^3;
ephem.pluto_i = deg2rad(i.a0 + i.a1*T + i.a2*T^2 + i.a3*T^3); %rad
ephem.pluto_Omega = deg2rad(Omega.a0 + Omega.a1*T + Omega.a2*T^2 + Omega.a3*T^3); %rad
ephem.pluto_LOP = LOP.a0 + LOP.a1*T + LOP.a2*T^2 + LOP.a3*T^3; %deg
ephem.pluto_w = deg2rad(ephem.pluto_LOP) - ephem.pluto_Omega; %deg
ephem.pluto_M = mod(rad2deg(deg2rad(ephem.pluto_L)-deg2rad(ephem.pluto_LOP)),360); %deg
Ccen = (2*ephem.pluto_e-ephem.pluto_e^3/4+5/96*ephem.pluto_e^5)*sind(ephem.pluto_M)+...
        (5/4*ephem.pluto_e^2-11/24*ephem.pluto_e^4)*sind(2*ephem.pluto_M) +...
         (13/12*ephem.pluto_e^3 - 43/64*ephem.pluto_e^5)*sind(3*ephem.pluto_M) + ...
          103/96*ephem.pluto_e^4*sind(4*ephem.pluto_M)+1097/960*ephem.pluto_e^5*sind(5*ephem.pluto_M);
ephem.pluto_nu = deg2rad(mod(rad2deg(deg2rad(ephem.pluto_M)+Ccen),360)); %rad

end