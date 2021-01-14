%
% REMOVE_TRANSITIONS This function removes data from a data sequence that
% is part of transitional data (i.e. during state changes or changes in class)
%
% function [data_out, class_out] =
%   remove_transitions(data_in, class_in, before, after)
%
% Author Adrian Chan
%
% This function removes transitional data. For example, if we have a series
% of feature vectors with class labels:
%
% [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3]'
%                      ^                 ^
% The ^ indicates where we have a transition. Data associated with this
% transition, along with the feature vectors before and after the
% transition will be removed (the number before and after are user
% specified)
%
% Inputs
%    data_in: rows of data (e.g. feature vectors)
%    class_in: class labels (column)
%    before: number of rows of data to remove before transitions
%            (default 8)
%    after: number of rows of data to remove after transitions
%            (default same as before)
%
% Outputs
%    data_out: rows of data with transitional data removed
%    class_out: class labels with transitional data removed
%
% Modifications
% 06/12/07 AC First created.
function [data_out, class_out] = remove_transitions(data_in, class_in, before, after)

if nargin < 4
    if nargin < 3
        before = 8;
    end
    after = before;
end

transition_index = diff([class_in(1);class_in]);
transition_index = find(transition_index ~= 0);
transition_index = [transition_index; length(class_in)];

N = length(transition_index);
Ndata = size(data_in,1);

data_out = [];
class_out = [];
for i = 1:(N-1)
    data_out = [data_out; data_in(min(transition_index(i)+after,Ndata):max(transition_index(i+1)-before,1),:)];
    class_out = [class_out; class_in(min(transition_index(i)+after,Ndata):max(transition_index(i+1)-before,1),:)];
end
