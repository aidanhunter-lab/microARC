function tc = yearday(t)
[yr,mo,dy,hr,mn,sc] = datevec(t);
ly = mod(yr,4) == 0; % leapyears
tc = datenum(ones(size(mo)),mo,dy,hr,mn,sc)-366;
tc_ly = datenum(zeros(size(mo)),mo,dy,hr,mn,sc);
tc(ly) = tc_ly(ly);
