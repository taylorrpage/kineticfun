function validInputs = checkRows(varargin)
%Checks if the number of rows in each input is equal.
%
%Returns boolean.

validInputs = true;
currentInput = varargin{1};
s1 = size(currentInput);

for i = 2:nargin
    s2 = size(varargin{i});
    validInputs = s1(1) == s2(1) & validInputs;
end %for