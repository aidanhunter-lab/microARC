function var = gapFill_forc(var, mask)
% fill gaps in forcing data that remained near seafloor after interpolation
% (use last value from layers above)
if nargin<nargin(@gapFill_forc), error('Insufficient input arguments.'); end
iNaN = isnan(var(:)) & mask(:);
if ~any(iNaN), return, end % no data gaps
while any(iNaN)
    v_ext = cat(2,var(:,1,:),var);
    iLast = ~isnan(var) & mask & isnan(v_ext(:,1:end-1,:));
    iFill = circshift(iLast,-1,2);
    var(iFill) = var(iLast);
    iNaN = isnan(var(:)) & mask(:);
end
