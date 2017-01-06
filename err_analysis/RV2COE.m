%Creaotr:Michael Von Hendy
%Date: 2/3/2013
%Name: get_COEs.m
%Purpose: ASEN 3200 homework #3.
%Modified: 9/15/2014 for ASEN 5050 homework #3.



function [a,e,i,w,Omega,nu] = RV2COE(r,v,mu)
%
%Input: r = position row vecotr in geocentric
%       v = velocity row vector in geocentric
%       mu = body gravitational parameter
%Output: Orbital parameters in [km, or degrees]

%Parameters:
khat = [0,0,1];
ihat = [1,0,0];
jhat = [0,1,0];
rhat = r/norm(r);

%Angular Momentum:
hvec = cross(r,v);
h = norm(hvec);
hhat = hvec/h;

%Energy:
epp = .5 * dot(v,v) - mu/norm(r);

%Geometry:
a = -mu/(2*epp);
e = sqrt(1+2*epp*h^2/mu^2);

%Orbit Orientation:
i = rad2deg(acos(dot(hhat,khat))); 
Omega = atan2(dot(hhat,ihat),-dot(hhat,jhat));
evec = 1/mu*cross(v,hvec) - rhat;
ehat = evec/e;
nhat = [cos(Omega),sin(Omega),0];
nhatperp = cross(hhat,nhat);
w = atan2(dot(nhatperp,ehat),dot(nhat,ehat));

%Position in Orbit:
ehatperp = cross(hhat,ehat);
nu = atan2(dot(r,ehatperp),dot(r,ehat));

%Convert to Degrees:
w = rad2deg(w)+(w < 0)*(360);
Omega = rad2deg(Omega)+(Omega < 0)*(360);
nu = rad2deg(nu)+(nu < 0)*(360);
if (nu < 0)
    nu = nu+360;
end

end