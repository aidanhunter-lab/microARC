function Forc = latestSampleTime(Forc, Data)

% Find the latest sampling time associated to each forcing data trajectory.
% This allows us to reduce integration times during parameter optimisation
% by eliminating unneccessary integration steps beyond the observation dates.

EventTraj = Data.scalar.EventTraj; % indexes which sampling events are associated to each trajectory
Forc.sampleTime = false(size(Forc.t));

for i = 1:Forc.nTraj
     j = find(EventTraj(:,i));
     sampleYearday = unique(Data.scalar.Yearday(ismember(Data.scalar.Event, j)));
     Forc.sampleTime(ismember(Forc.Yearday(:,i), sampleYearday), i) = true;
end

Forc.integrateFullTrajectory = true; % may change to false for efficiency when optimising parameters
