%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Before you use call this function, 
%       I recommand you to set your figure's xlim and ylim (i.e., axis range)
%       and size of the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p] = plotDash(h, x, y, spec_, num_, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   h: Handle of the figure
%       ex) h = figure(1)
%   x: data x
%   y: data y
%   spec_: spec of line - (relative length of)(solid, blank, ...)
%          Default option: -1 --> dashed line
%       ex) spec_ = [1, 0.5] || spec_ =  [10, 5] || ...
%       ex) spec_ = [1, 0.1] || ...
%       ex) spec_ = [1, 0.1, 0.1, 0.1] || ...
%   num_: How many times that pattern will be repeated for xlim or ylim
%   vargargin: "Pair(s)" of options for plot
%       ex) varargin{1} = 'color', varargin{2} = 'g'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   User input - Points per pattern
    n_node = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Size of x
    nx = length(x);
    ny = length(y);
    if (nx ~= ny)
        p = plot(0, 0, 'o');
        delete(p);
        return;
    else
        n = nx;
    end
    
    %   Size of plot
    tmp = get(h, 'InnerPosition');
    W = tmp(3);
    D = tmp(4);
    SF = W / 560;
    R  = D / W * 28.55 / 27.1;
    if (sum(size(h.CurrentAxes)) == 0)
        p = plot(0, 0);
        delete(p);
        DX = max(x) - min(x);
        DY = max(y) - min(y);
    else
        tmp = h.CurrentAxes.XLim;
        DX = tmp(2) - tmp(1);
        tmp = h.CurrentAxes.YLim;
        DY = tmp(2) - tmp(1);
        clear tmp
    end
    
    %   Normalization
    r = zeros(n, 2);
    for i = 1:n
        r(i, 1) = x(i) / DX * SF;
        r(i, 2) = y(i) / DY * SF * R;
    end
       
    %   Curve length
    s  = zeros(n    , 1);
    ds = zeros(n - 1, 1);
    for i = 1:(n - 1)
        ds(i)    = norm(r(i + 1, :) - r(i, :));
        s(i + 1) = s(i) + ds(i);
    end
    
    %   Linearly splined data - User input: 50
    ds_min = min(ds);
    ds_min = min(ds_min, 1/(n_node*num_));
    s_s = 0:ds_min:s(end);
    x_s = interp1(s, x, s_s);
    y_s = interp1(s, y, s_s);
    
    %   Default line spec = dashed
    if (spec_(1) < 0)
        spec_ = [10, 6];
        num_ = 26;
    end
    
    %   spec_: absolute --> relative
    spec_ = spec_ / sum(spec_) / num_;
    m = length(spec_);
    l_spec = sum(spec_);
    
    %   Set T/F flag for index
    k = length(s_s);
    flag_idx = true(k, 1);
    for i = 1:k
        curr = rem(s_s(i), l_spec);
        for j = 2:2:m
            flag = (curr > sum(spec_(1:(j-1))));
            flag = (curr < sum(spec_(1:(j-0)))) && flag;
            
            if (flag)
                flag_idx(i) = false;
            end
        end
        flag_idx(i);
    end
    
    %   Set NaN to x_s and y_s
    x_s(~flag_idx) = NaN;
    y_s(~flag_idx) = NaN;
        
    %   Plot graph
    p = plot(h.CurrentAxes, x_s, y_s);
    
    %   Apply nargin
    for i = 1:2:(nargin - 5)
        set(p, varargin{i}, varargin{i + 1});
    end
    
end

