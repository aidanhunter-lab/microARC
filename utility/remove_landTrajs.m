function struc = remove_landTrajs(struc, iMask)
% remove trajectories that go on land at any time
if nargin<nargin(@remove_landTrajs), error('Insufficient input arguments.'); end
fNames = fieldnames(struc);
% get dimensions
for i = 1:length(fNames)
    if  ndims(struc.(fNames{i}))==3
        nx = size(struc.(fNames{i}),3);
        break
    end
end
% filter data
for i = 1:length(fNames)
    if  size(struc.(fNames{i}),1)==nx
        struc.(fNames{i}) = struc.(fNames{i})(iMask,:,:);
    elseif size(struc.(fNames{i}),2)==nx
        struc.(fNames{i}) = struc.(fNames{i})(:,iMask);
    elseif size(struc.(fNames{i}),3)==nx
        struc.(fNames{i}) = struc.(fNames{i})(:,:,iMask);
    end
end
