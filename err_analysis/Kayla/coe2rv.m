function [r_ijk,v_ijk] = coe2rv(a,e,i,Omega,w,nu)
% This function takes  the inputs a, e, i, Omega, w, and nu and
% converts it to position and velocity vectors.

% Variables:
% a - semi-major axis (km)
% e - eccentricity (unitless)
% r - position vector (km)
% v - velocity vector (km/s)
% h_bar - angular momentum vector (km^2/s)
% n_bar - vector pointing to the node
% k_bar - unit vector [0,0,1]
% mu - Gravitational parameter of Earth! (km^3/s^2)
% xi - energy (km^2/s^2)

% Rotation matricies

% ROT1(x) = [1   0    0  ]
%           [0 cosX  sinX]  
%           [0 -sinX cosX] 

% ROT2(x) = [cosX 0  -sinX],
%           [0    1     0 ]  
%           [sinX 0   cosX]

% ROT3(x) = [ cosX sinX  0 ]
%           [-sinX cosX  0 ]  
%           [ 0     0    1 ]


mu = 3.986004415e5;
p = a*(1-e^2);

% convert to radians 
nu = deg2rad(nu);
i = deg2rad(i);
Omega = deg2rad(Omega);
w = deg2rad(w);

% compute rpqw and vpqw
rpqw = [p*cos(nu)/(1+e*cos(nu));p*sin(nu)/(1+e*cos(nu));0];
vpqw = [-sqrt(mu/p)*sin(nu);sqrt(mu/p)*(e+cos(nu));0];

% Convert to rijk and vijk
ijk_pqw = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i),...
           -cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i),...
            sin(Omega)*sin(i);...
            sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i),...
            -sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i),...
            -cos(Omega)*sin(i);...
            sin(w)*sin(i),cos(w)*sin(i),cos(i)];
        
r_ijk = ijk_pqw*rpqw;
v_ijk = ijk_pqw*vpqw;



