function output = qStats(y,x,opts)
% output = qStats(y,x,opts)
%
% This is a wrapper function for grpstats(), which I *think* is in the
% statistics toolbox.
% 
% y = DV, x = IVs (multiple IVs in separate columns), opts = 1: just mean,
% opts = 2: mean and conditions, opts = 3: mean, conds, and sem, etc -- see
% the code below for the full set of options
% 
% WJH willjharri@gmail.com

if nargin < 3
    opts = 3;
end

if opts == 1
    output = grpstats(y, x, {'mean'});
elseif opts == 2
    [m,g] = grpstats(y, x, {'mean','gname'});
    output = [str2double(g) m];
elseif opts == 3
    [m,g,s] = grpstats(y, x, {'mean','gname','sem'});
    output = [str2double(g) m s];
elseif opts == 4
    [m,g,s,sd] = grpstats(y, x, {'mean','gname','sem','std'});
    output = [str2double(g) m s sd];
elseif opts == 5
    [m,g,s,n] = grpstats(y, x, {'mean','gname','sem','numel'});
    output = [str2double(g) m s n];
elseif opts == 6
    [m,g,s,n,sd] = grpstats(y, x, {'mean','gname','sem','numel','std'});
    output = [str2double(g) m s n sd];
end