function varOut = interp_forc(varIn, xIn, xOut, tIn, tOut, zIn, zOut)
% interpolate forcing data to model dimensions
if  ismatrix(varIn)
    % 2D data: horizontal space and time
    if  nargin<nargin(@interp_forc)-3
        error('Insufficient input arguments for 2D data. Need exactly %i.', nargin(@interp_forc)-2);
    end
    varOut = reshape(interpn(tIn, xIn, varIn, tOut, xOut), [length(tOut), 1, length(xOut)]);
else
    if  nargin<nargin(@interp_forc)
        error('Insufficient input arguments for 3D data. Need exactly %i.', nargin(@interp_forc));
    end
    % 3D data: horizontal/vertical space and time
    varOut = interpn(tIn, zIn, xIn, varIn, tOut, zOut, xOut);
end
