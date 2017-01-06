function plot_handle = plot_sphere( r,position,normalizer,planet_scale,col )
%plot_sphere.m plots a sphere.
%
%Inputs: 
%   r = sphere radius [km]
%   position = Helicentric Inertial position vector [x,y,z] [km]
%   normalizer = position normalization factor. Use Radius of Earth orbit for AU
%   planet_scale = scale factor for sphere size
%   col = color specification eg 'y' for solid yellow sphere

%Calculate spherical surface:
[x,y,z] = sphere(81);
x = x*r/planet_scale;
y = y*r/planet_scale;
z = z*r/planet_scale;

X=position(1);
Y=position(2);
Z = position(3);
%Plot body:
plot_handle = surf(x+X/normalizer,y+Y/normalizer,z+Z/normalizer,...
    'XDataSource','X',...
    'YDataSource','Y',...
    'ZDataSource','Z')
set(plot_handle,'FaceColor',col)
set(plot_handle,'linestyle','none')

end

