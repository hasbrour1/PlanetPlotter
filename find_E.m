
e = .321;
M= deg2rad(50);

x = true;

if -pi < M < 0 || M > pi
    En = M-e;
else
    En = M + e;
end

 count = 0;
 while x == true
    En1 = En + (M-En+e*sin(En))/(1-e*cos(En));
    count = count +1;
    if abs(En1-En) < 0.0001
        x = false;
    end
    if count > 10000
        break;
    end 
    En = En1
 end
 E1 = mod(rad2deg(En1),360)
