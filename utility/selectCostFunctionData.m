function Data = selectCostFunctionData(Data, obsScalar, obsSize)

Data.scalar.obsInCostFunction = obsScalar;
Data.scalar.inCostFunction = ismember(Data.scalar.Variable, obsScalar);

% Omit any events where data used in the cost function was not collected
fields = fieldnames(Data.scalar);
uev = unique(Data.scalar.Event);
for i = 1:length(uev)
    ind = Data.scalar.Event == uev(i);
    x = unique(Data.scalar.Variable(ind));
    y = ismember(x, Data.scalar.obsInCostFunction);
    if ~any(y)
        for j = 1:length(fields)
            if size(Data.scalar.(fields{j}), 1) > 1
                Data.scalar.(fields{j}) = Data.scalar.(fields{j})(~ind);
            end
        end
    end
end
if length(unique(Data.scalar.Event)) < Data.scalar.nEvents
    Data.scalar.nSamples = size(Data.scalar.Value,1);    
    Data.scalar.nEvents = length(unique(Data.scalar.Event));    
    oldEventLabs = unique(Data.scalar.Event);
    newEventLabs = (1:length(oldEventLabs))';    
    x = table(oldEventLabs, newEventLabs);
    for i = 1:height(x)
        Data.scalar.Event(Data.scalar.Event == x.oldEventLabs(i)) = x.newEventLabs(i);
    end
end


Data.size.obsInCostFunction = obsSize;
Data.size.dataBinned.inCostFunction = ismember(Data.size.dataBinned.Variable, obsSize);

Data.sizeFull.obsInCostFunction = obsSize;
Data.sizeFull.dataBinned.inCostFunction = ismember(Data.sizeFull.dataBinned.Variable, obsSize);

