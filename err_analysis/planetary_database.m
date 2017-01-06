function  planetary_database
%planetary_database.m uploads planetary data for the solar system into
%structures for use in MATLAB

global sun moon mercury venus earth mars jupiter saturn uranus neptune pluto;

%Definitions:
    %a = semi-major axis of orbit relative to Heliocenter [km]
    %e = eccentricity of orbit
    %i = inclination of orbit [deg]
    %RAAN = Right Ascension of Ascending Node of Orbit [deg]
    %LOP = Longitude of Perihelion [deg]
    %T = orbital period [yr]
    %v = orbital velocity [km/s]
    %mu = body graviational parameter [km^3/s^2]
    %L = mean longitude at J2000 epoch [deg]

%Lunar:
%Note lunar body orbits earth, not sun
moon.a = 384400; %km
moon.e = 0.05490;
moon.i = 5.145396; %deg
moon.RAAN = NaN;
moon.LOP = NaN;
moon.T = 0.0748; %yr
moon.v = 1.0232; %km/s
moon.r = 1738.0; %km
moon.mu = 4902.799; %km^3/s^2

%Mercury:
mercury.a = 57909083; %km
mercury.e = 0.205631752; 
mercury.i = 7.00498625; %deg
mercury.RAAN = 48.33089304; %deg
mercury.LOP = 77.45611904; %deg
mercury.T = 0.24084445; %yr
mercury.v = 47.8725; %km/s
mercury.r = 2439.0; %km
mercury.mu = 2.2032e4; %km^3/s^2
mercury.L = 252.25084; %deg

%Venus:
venus.a = 108208601; %km
venus.e = 0.006771882;
venus.i = 3.39446619; %deg
venus.RAAN = 76.67992019; %deg
venus.LOP = 131.56370724; %deg
venus.T = 0.61518257; %yr
venus.v = 35.0214; %km/s
venus.r = 6052.0; %km
venus.mu = 3.257e5; %km^3/s^2
venus.L = 181.97973; %deg

%Earth:
earth.a = 149598023; %km
earth.e = 0.016708617; 
earth.i = 0; %deg
earth.RAAN = 0; %deg
earth.LOP = 102.93734808; %deg
earth.T = 0.99997862; %yr
earth.v = 29.7859; %km/s
earth.r = 6378.1363; %km
earth.mu = 3.986004415e5; %km^3/s^2
earth.L = 100.46435; %deg

%Mars:
mars.a = 227939186*8/9;
mars.e = 0.093400620;
mars.i = 1.84972648;
mars.RAAN = 49.55809321;
mars.LOP = 336.06023398;
mars.T = 1.8871105; 
mars.v = 24.1309;
mars.r = 3397.2;
mars.mu = 4.305e4;
mars.L = 355.45332;

%Jupiter:
jupiter.a = 778298361;
jupiter.e = 0.048494851;
jupiter.i = 1.30326966;
jupiter.RAAN = 100.46444064;
jupiter.LOP = 14.33130924;
jupiter.T = 11.856525;
jupiter.v = 13.0697;
jupiter.r = 71492.0;
jupiter.mu = 1.268e8;
jupiter.L = 34.40438;

%Saturn:
saturn.a = 1429394133; %km
saturn.e = 0.055508622;
saturn.i = 2.48887810; %deg
saturn.RAAN = 113.66552370; %deg
saturn.LOP = 93.05678728; %deg
saturn.T = 29.423519; %yr
saturn.r = 60268; %km
saturn.mu = 3.794e7; %km^3/s^2
saturn.L = 49.94432;

%Uranus:
uranus.a = 2875038615; %km
uranus.e = 0.046295898; 
uranus.i = 0.77319617; %deg
uranus.RAAN = 74.00594723; %deg
uranus.LOP = 173.00515922; %deg
uranus.T = 83.747406; %yr
uranus.r = 25559; %km
uranus.mu = 5.794e7; %km^3/s^2
uranus.L = 313.23218; %deg

%Neptune:
neptune.a = 4504449769;
neptune.e = 0.008988095; 
neptune.i = 1.76995221; 
neptune.RAAN = 131.78405702;
neptune.LOP = 48.12369050;
neptune.T = 163.7232045;
neptune.r = 24764;
neptune.mu = 6.809e6;
neptune.L = 304.88003; %deg

%Pluto:
pluto.a = 5915799000;
pluto.e = 0.249050;
pluto.i = 17.14216667;
pluto.RAAN = 110.29713889;
pluto.LOP = 224.13486111;
pluto.T = 248.0208;
pluto.r = 1151;
pluto.mu = 900;
pluto.L = 238.92881;

end

