function forc = filterTrajectoryJumps(forc, maxVel)
    % =====================================================================
    % function forc = filterTrajectoryJumps(forc, maxJump)
    % =====================================================================
    % This function filters trajectories that have significant jumps
    % between subsequent locations, which can occur in 'particulator'
    %
    % INPUT : forc   ... forcing structure as produced by 'particulator' or
    %                    by sandbox models
    %         maxVel ... maximum velocity ((in m/s) between two subsequent
    %                    locations of each trajectory; default: 1.5 m/s
    %
    % OUTPUT: forc   ... updated forcing structure; only including
    %                    trajectories without jumps
    % =====================================================================
    % Fabian Grosse, U Strathclyde, 6 February 2020
    % =====================================================================
    
    if nargin<nargin(@filterTrajectoryJumps), maxVel = 1.5; end
    if nargin==0, error('Need to provide a forcing structure.'); end
    
    % calculate maximum distance based on time step and maximum speed
    dt = diff(forc.t(1:2)) * 86400; % in seconds
    maxDist = maxVel*dt;            % distance in metres
    
    % find trajectories to be kept
    x = squeeze(forc.x);
    y = squeeze(forc.y);
    nx = size(x,2);
    iKeep = false(1,nx);
    for ix = 1:nx
        d = distance_on_earth(deg2rad([x(1:end-1,ix), y(1:end-1,ix)]), ...
                              deg2rad([x(2:end  ,ix), y(2:end  ,ix)]));
        iKeep(ix) = max(d)<=maxDist;
    end
    
    % filter trajectories
    flds = fieldnames(forc);
    for i = 1:length(flds)
        if  strcmpi(flds{i}, 'profiles')
            % treat data in 'profiles' sub-structure
            pFlds = fieldnames(forc.(flds{i}));
            for j = 1:length(pFlds)
                xDim = size(forc.(flds{i}).(pFlds{j}))==nx;
                if ~any(xDim), continue, end
                if  ismatrix(forc.(flds{i}).(pFlds{j}))
                    forc.(flds{i}).(pFlds{j}) = forc.(flds{i}).(pFlds{j})(:,iKeep);
                else
                    forc.(flds{i}).(pFlds{j}) = forc.(flds{i}).(pFlds{j})(:,iKeep,:);
                end
            end
        else
            % treat data in main structure
            xDim = size(forc.(flds{i}))==nx;
            if ~any(xDim), continue, end
            if  ismatrix(forc.(flds{i}))
                forc.(flds{i}) = forc.(flds{i})(:,iKeep);
            else
                if  find(xDim)==2
                    forc.(flds{i}) = forc.(flds{i})(:,iKeep,:);
                else
                    forc.(flds{i}) = forc.(flds{i})(:,:,iKeep);
                end
            end
        end
    end
    
end