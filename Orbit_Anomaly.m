function [nu,E,M,t] = Orbit_Anomaly(a,e,nu,E,M,t)
% This function solves for a desired anomaly given one known input of 
% time, Ecentric anomaly, Mean anomaly, or True anomaly. Assume the 
% orbit is elliptical. 
%
% Variables:
% a - semi-major axis (km)
% e - eccentricity (unitless)
% nu - True anomaly (deg)
% E - Eccentric Anomaly (deg)
% M - Mean Anomaly (deg)
% t - Time since pariapse(t-tp)(s)
% n - mean angular rate of change (rad/s)
%
% mu - Gravitational parameter of Earth! (km^3/s^2)

mu = 1.32712428e11;
n = sqrt(mu/a^3);

% Given Nu
if ~isempty(nu) 
    nu = deg2rad(nu);
    % calculate E
    E = 2*atan2(tan(nu/2),sqrt((1+e)/(1-e)));
    
    if E<0
        E = E+2*pi;
    end
    %calculate M
    M = E-e*sin(E);
    if M<0
        M = M+2*pi;
    end
    %calculate t
    t = M/n;
    
    nu = rad2deg(nu);
    E = rad2deg(E);
    M= rad2deg(M);
end 

%Given E
if ~isempty(E)
    E = deg2rad(E);
    if E<0
        E = E+2*pi;
    end
   % calculate nu
    nu = 2*atan2(sqrt((1+e)/(1-e))*tan(E/2),1); 
    if nu<0
        nu = nu+2*pi;
    end
    %calculate M
    M = E-e*sin(E);
    if M<0
        M = M+2*pi;
    end
    %calculate t
    t = M/n;
    
    nu = rad2deg(nu);
    E = rad2deg(E);
    M= rad2deg(M);
end 

%Given M
if ~isempty(M)
    M= deg2rad(M);
    if M<0
        M = M+2*pi;
    end
    %calculate t
    t = M/n;
    
    %calculate E
    E_equ = @(EE) EE-e*sin(EE)-M;
    E = fzero(E_equ,M);
    if E<0
        E = E+2*pi;
    end
    %calculate nu
    nu = 2*atan2(sqrt((1+e)/(1-e))*tan(E/2),1);
    if nu<0
        nu = nu+2*pi;
    end
    nu = rad2deg(nu);
    E = rad2deg(E);
    M= rad2deg(M);
end 

%Given t
if ~isempty(t)
    %calcualte M
    M = n*t;
    if M<0
        M = M+2*pi;
    end
    %calculate E
     E_equ = @(EE) EE-e*sin(EE)-M;
     E = fzero(E_equ,M);
    if E<0
        E = E+2*pi;
    end
   % x = true;

   %if -pi < M < 0 or M > pi
    % E = M-e;
% else
%     E = M + e;
% end
 
 %count = 0;
 %while x == true
 %    En = M;
 %    En1 = En + (M-En+e*sin(En))/(1-e*cos(En));
  %   count = count +1;
  %   if (En1-En) < 0.0001
  %       x = false;
  %   end
  %   if count > 10000
  %       break;
  %end 
     
 %end
 %    E = En1;
     if E<0
        E = E+2*pi;
     end
    
    %calculate Nu
    nu = 2*atan2(sqrt((1+e)/(1-e))*tan(E/2),1);
    if nu<0
        nu = nu+2*pi;
    end
    nu = rad2deg(nu);
    E = rad2deg(E);
    M= rad2deg(M);
end

end

