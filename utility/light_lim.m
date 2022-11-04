function [pc,I_lim,dl] = light_lim(a_Chl, psat, Isurf, att0, lat, yd, zw)
%LIGHT_LIM calculates the depth- and time averaged light limitation,
% considering the law of Lambert & Beer (for vertical average) while  
% also assuming a triangular diurnal light cycle. 
% It requires the 'expint' (exponential-integral) function for the 
% vertical integral, as well as daylength (dl) for the approximate 
% temporal average. 
% INPUT:
%  psat      = light saturated photosynthetic rate (1 / day)
%  a_Chl     = Chl-dependent initial slope of photon capture rate 
%              (unit = inverse irradiance / day)
%  psat      = maximum (light saturated) growth rate 
%  Isurf     = irradiance (light) at surface
%  att0      = light attenuation within individual layers (1 / m)
%  lat       = latitude
%  yd        = day of year (year day)
%  zw        = widths of the vertical boxes (m) 
%
% OUTPUT:
% pc         = photosynthetic (carbon production) rate (1 / day) 
% I_lim      = light limitation (depth- and time averaged light limitation) 
% dl         = daylength (fraction of 24 hours)

% light profile (at upper bounds of depth layers)
nz      = length(att0);
att     = att0 .* zw';
att     = [0 cumsum(att(1:nz))];
e(1,:)  = exp(-att(1:nz));
e(2,:)  = exp(-att(2:nz+1));
I       = Isurf .* e;  
% I       = max(Isurf,1e3) .* e; % low threshold is approx. 1e-2 muEin/s/m^2  

% determine daylength
dl      = daylen(lat*pi/180, yd);

% factor for exponent
fa      = 2 .* a_Chl ./ psat;

% Convert from daily mean (I) to mean daytime irradiance (mI)
if dl>0
    mI = I / dl;
    dl_ = dl;
else
    mI = I;
    dl_ = 1; % this is for the back conversion from daytime mean to daily mean (see line 55)
end

% mI

% NOTE, Ei( x) = -expint(-x), for real x > 0 
% -->   Ei(-x) = -expint( x)  (?)  
I_lim1  = 1 - ( expint(fa.* mI(2,:)) - expint(fa.* mI(1,:)) ) ./ (att0 .* zw' );
I_lim2  = ( (1-exp(-fa.*mI(2,:)))./mI(2,:) - (1-exp(-fa.*mI(1,:)))./mI(1,:) ) ./ ( fa .* att0 .* zw' );

% Since the I_lim1 and I_lim2 calculations represent daytime averages (means),  
% we have to convert them back to daily means by multiplying with the daylength 
% I_lim   = dl * (I_lim1 - I_lim2);
I_lim   = max(2e-16, dl_ * (I_lim1 - I_lim2));   % 2e-16
pc      = psat .* I_lim;

    % nested function
    function dl = daylen (lat, doy)
    %DAYLEN: daylength in decimal days after Brock (1981)
    % DL = DAYLEN(LAT, DOY)
    %   LAT     latitude in rad
    %   DOY     day of year
    % This function assumes that a year has exactly 365 days
      rad = pi/180;
      fd = rad*360/365;
      tilt = rad*23.45;
      d1 = tilt.*sin(fd.*(284 + doy));            % declination
      cw1 = -tan(lat).*tan(d1);                   % cos(sunset hour-angle)
      dl = acos(sign(cw1).*min(abs(cw1), 1))./pi;
    end
end
