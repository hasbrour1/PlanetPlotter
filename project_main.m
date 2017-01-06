function [r,v,ephem] = project_main(JD,doy,date,date_time)
% This script takes a user input date and outputs the planets position and
% velocity. Pluto is excluded from the planets. Year bounded between
% 1900-2100.

% Input Variables
% JD - Julian Date (UTI) An integer input assumes zeros after decimal
% doy - day of year (1-365),year
% date - MM/DD/YYYY (assume time is noon)
% date_time - MM/DD/YYYY/hh/mm/ss

% Output Variables
% r - position vector of all planets
% v - velocity vector of all planets
% ephem - ephemeris, planetary position data 

% For an input of JD (xxxxxxx.xx)
if ~isempty(JD)
    [r_planet,v_planet,ephem] = meeus2(JD)
end 

% For an input of days (xxx.xxxx,yyyy)
if ~isempty(doy) 
    fid = fopen('data.txt')
    %data = csvread('data.csv');
    data = fread(fid,100);
    yr = data(2);
    num_days = data(1);
    % Check if Leap Year
    if rem(yr,4) == 0
    	tmp = [31 60 91 121 152 182 213 244 274 305 335 366.25];
          if yr == 2000
                tmp = [31 59 90 120 151 181 212 243 273 304 334 365];
          end
    else
        tmp = [31 59 90 120 151 181 212 243 273 304 334 365.25];
    end
    %Find month and days:
    tmp1 = (floor(num_days)./tmp);
    month = find((tmp1== max(tmp1(tmp1<=1)))==1);
    if month == 1
        day = floor(num_days);
    else
        day = floor(num_days-tmp(month-1));
    end
    hours = (num_days-floor(num_days))*24; % GMT
    minutes = (hours-floor(hours))*60;
    seconds = (minutes-floor(minutes))*60;
    
    % convert to JD
    JD = 367*yr-floor((7*(yr+floor((month+9)/12)))/(4))+...
    floor(275*month/9)+day+1721013.5+((seconds/60+minutes)/60+hours)/24;
    
    % compute terms
    [r_planet,v_planet,ephem] = meeus2(JD)
end 

% For an input of date (MM/DD/YYYY)
if ~isempty(date)
    data = ('data.csv');
    yr = str2num(data(7:10));
    month = str2num(data(1:2));
    day = str2num(data(4:5));
    hr = 12;
    min = 0;
    sec = 0;
    
    % convert to JD
    JD = 367*yr-floor((7*(yr+floor((month+9)/12)))/(4))+...
    floor(275*month/9)+day+1721013.5+((sec/60+min)/60+hr)/24;

    % compute terms
    [r_planet,v_planet,ephem] = meeus2(JD)
end 

% For an input of date_time (MM/DD/YYYY/hh/mm/ss)
if ~isempty(date_time)
    yr = ;
    month = ;
    day = ;
    hr = ;
    min = ;
    sec = ;
    
    % convert to JD
    JD = 367*yr-floor((7*(yr+floor((month+9)/12)))/(4))+...
    floor(275*month/9)+day+1721013.5+((sec/60+min)/60+hr)/24;
end 
end