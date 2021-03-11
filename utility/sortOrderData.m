function out = sortOrderData(dat)

obs = dat.obsInCostFunction;

% loop through data types
for i = 1:length(obs)
    ind = strcmp(dat.Variable, obs{i});
    y = dat.scaled_Value(ind);
    [~, yo] = sort(y);
    dat.(['sortOrder_' obs{i}]) = yo;
end

out = dat;
