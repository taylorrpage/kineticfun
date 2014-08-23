function newArray = populateconstants(array)
%Replaces NaN values in array with the first value in the corresponding column.
%
%Example:
%arr = 
%   3    4    10
%   2  NaN     8
%  10    7   NaN
% NaN  NaN     5
%
%ans = populateconstants(arr)
%ans = 
%   3    4    10
%   2    4     8
%  10    7    10
%   3    4     5

for i = 1:numel(array)
    if isnan(array(i))
        [~, n] = ind2sub(size(array), i);
        array(i) = array(1, n);
    end %if
end %for

newArray = array;