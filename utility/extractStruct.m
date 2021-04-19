function [] = extractStruct(z, varargin)
% Assign each variable stored in struct z to the worksapce
workspace = 'caller';
if ~isempty(varargin)
    i = strcmp(varargin, 'workspace');
    if any(i), workspace = varargin{find(i) + 1}; end
    if ~ismember(workspace, {'base','caller'})
        warning('The "workspace" argument in extractStruct.m should be either "base" or the default "caller". It has been reset to the default value.')
        workspace = 'caller';
    end
end
fields = fieldnames(z);
for j = 1:length(fields)
    assignin(workspace, fields{j}, z.(fields{j}))
end