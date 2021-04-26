function Forc = latestSampleTime(Forc, Data)

% Find the latest sampling time associated to each forcing data trajectory.
% This allows us to reduce integration times during parameter optimisation
% by eliminating unneccessary integration steps beyond the observation dates.

EventTraj = Data.scalar.EventTraj; % indexes which sampling events are associated to each trajectory
Forc.sampleTime = false(size(Forc.t));

for i = 1:Forc.nTraj
     j = find(EventTraj(:,i));
%      sampleYearday = unique(Data.scalar.Yearday(ismember(Data.scalar.Event, j)));
     sampleYearday = max(Data.scalar.Yearday(ismember(Data.scalar.Event, j)));
     Forc.sampleTime(ismember(Forc.Yearday(:,i), sampleYearday), i) = true;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The size data may be grouped into trajectories originating from Arctic or
% Atlantic waters, so adjust the final sample times to accomodate the
% latest times for each water mass. This is only a slight adjustment.
if isfield(Forc, 'waterMass')
    waterMasses = unique(Forc.waterMass);
    for i = 1:length(waterMasses)
        waterMass = waterMasses{i};
        ind = strcmp(Forc.waterMass, waterMass);
        sampleTimes = Forc.sampleTime(:,ind);
        latest = find(any(sampleTimes, 2), 1, 'last');
        sampleTimes = false(size(sampleTimes));
        sampleTimes(latest,:) = true;
        Forc.sampleTime(:,ind) = sampleTimes;
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Forc.integrateFullTrajectory = true; % may change to false for efficiency when optimising parameters
