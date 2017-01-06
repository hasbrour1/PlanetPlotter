function planet_ellipse( planet,normalizer,col,ls )
%palnet_ellipse.m plots the 2-D orbit of a planet in 3-D space 
%
%Inputs: 
%   planet = planetary data structure for current body
%   normalizer = normalization factor. Use Radius of Earth orbit for AU
%   col = color specificaiton eg 'r' for red
%   ls = line specification eg '--' for dashed


%Calculate Ellipse Values:
RAAN = planet.RAAN;
i = deg2rad(planet.i);
LOP = deg2rad(planet.LOP);
a = planet.a;
e = planet.e;
b = a*sqrt(1-e^2);
c = a*e;

%Parametrization:
t0 = 0;
t = linspace(t0,2*pi+t0,1000);
X = a*cos(t);%+c;
Y = b*sin(t);
Z = zeros(1,length(X)); 
Rpqw2ijk = [cos(RAAN)*cos(LOP) - sin(RAAN)*sin(LOP)*cos(i),...
    -cos(RAAN)*sin(LOP)-sin(RAAN)*cos(LOP)*cos(i),...
    sin(RAAN)*sin(i);...
    sin(RAAN)*cos(LOP)+cos(RAAN)*sin(LOP)*cos(i),...
    -sin(RAAN)*sin(LOP)+cos(RAAN)*cos(LOP)*cos(i),...
    -cos(RAAN)*sin(i);...
    sin(LOP)*sin(i), cos(LOP)*sin(i),cos(i)];
for i = 1:length(t)
    Pos(i,:) = Rpqw2ijk*[X(i);Y(i);Z(i)];
end


%Ellipse centered at sun:
x = Pos(:,1)/normalizer;
y = Pos(:,2)/normalizer;
z = Pos(:,3)/normalizer;

%Plot body:
plot3(x,y,z,'Color',col,'LineStyle',ls)
end

