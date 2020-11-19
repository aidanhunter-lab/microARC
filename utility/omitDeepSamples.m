function dat = omitDeepSamples(Data, FixedParams)

% Discard in-situ samples collected from below the maximum modelled depth

dat = Data;
scalarData = Data.scalar;

remove = scalarData.Depth > max(abs(FixedParams.z));
if any(remove)
    fields = fieldnames(scalarData);
    for i = 1:length(fields)
        if size(scalarData.(fields{i}), 1) > 1
            scalarData.(fields{i})(remove) = [];
        end
    end
    scalarData.nSamples = length(scalarData.Year);
    scalarData.nEvents = length(unique(scalarData.Event));
end

dat.scalar = scalarData;
