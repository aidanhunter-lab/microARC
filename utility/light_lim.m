function [pc,I_lim,dl] = light_lim(a_Chl, psat, Isurf, att0, lat, yd, zw, varargin)
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
%  Setting optional argument 'dist'='flat' calculates only the
%  depth-averaged light limitation using daily mean irradiance.
%
% OUTPUT:
% pc         = photosynthetic (carbon production) rate (1 / day) 
% I_lim      = light limitation (depth- and time averaged light limitation) 
% dl         = daylength (fraction of 24 hours)

extractVarargin(varargin)
if ~exist('dist', 'var')
    dist = 'triangular'; % triangular light distribution by default
end

% light profile (at upper bounds of depth layers)
nz      = length(att0);
att     = att0 .* zw';
att     = [0 cumsum(att(1:nz))];
e(1,:)  = exp(-att(1:nz));
e(2,:)  = exp(-att(2:nz+1));
I       = Isurf .* e;  
% I       = max(Isurf,1e3) .* e; % low threshold is approx. 1e-2 muEin/s/m^2  

% factor for exponent
fa      = a_Chl ./ psat;

% determine daylength
dl      = daylen(lat*pi/180, yd);
dl_ = dl;
switch dist
    case 'flat'
        dl_ = 1;
    otherwise
        fa = 2 .* fa;
end

% Convert from daily mean (I) to mean daytime irradiance (mI) (unless dist=flat)
if dl>0
    mI = I / dl_;
else
    mI = I;
    dl_ = 1; % this is for the back conversion from daytime mean to daily mean (see line 73)
end

% Depth integral
I_lim1  = 1 - (expint(fa.* mI(2,:)) - expint(fa.* mI(1,:))) ./ (att0 .* zw');
% NOTE, Ei( x) = -expint(-x), for real x > 0 
% -->   Ei(-x) = -expint( x)  (?)

switch dist, case 'flat'
    I_lim = min(max(I_lim1, 0), 1);
    otherwise
        % Time integral
        I_lim2  = ((1-exp(-fa.*mI(2,:)))./mI(2,:) - (1-exp(-fa.*mI(1,:)))./mI(1,:)) ./ (fa .* att0 .* zw');
        % Since the I_lim1 and I_lim2 calculations represent daytime averages (means),
        % we have to convert them back to daily means by multiplying with the daylength
        I_lim   = min(max(dl_ * (I_lim1 - I_lim2), 0), 1);
end

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
