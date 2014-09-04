function varargout = xformparam(varargin)
%[paramVector, nFixedP, nCols] = xformparam(paramMatrix, fixedParam);
%   Transforms a matrix and vector into a single vector and returns the number
%   of elements in the original vector, nFixedP, and the number of columns in
%   the original matrix, nCols.
%
%[paramMatrix, fixedParam] = xformparam(paramVector, nFixedP, nCols);
%   Transforms a single vector into a matrix, paramMatrix, with nCols columns, 
%   and a vector, fixedParam with nFixedP elements.

if nargin == 2
paramMatrix = varargin{1};
fixedParam  = varargin{2};

[~, nCols]     = size(paramMatrix);
paramMatrix    = paramMatrix';
nFixedP        = length(fixedParam);
paramVector    = [fixedParam(:); paramMatrix(:)];

varargout = {paramVector, nFixedP, nCols};

elseif nargin == 3
paramVector = varargin{1};
nFixedP     = varargin{2};
nCols       = varargin{3};

fixedParam  = paramVector(1:nFixedP)';
paramVector = paramVector(nFixedP+1:end);
paramMatrix = vec2mat(paramVector, nCols);

varargout   = {paramMatrix, fixedParam};

else
    error('Wrong number inputs')
end %if


end %xformparam