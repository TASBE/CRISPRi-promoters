%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Before you use call this function, 
%       DisplayName should be given for p_
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p] = plotLegendBoundary(h, x1_, x2_, y1_, y2_, varargin)
%   h: Handle of figure
%   x_1: Relative x percent location of rectangle
%   x_2: Relative x percent location of rectangle
%   x_1: Relative y percent location of rectangle
%   x_2: Relative y percent location of rectangle
%   varargin: "Pair(s)" of options for text object

    %   Initialization
    hold(h.CurrentAxes, 'on');

    %   Determine positions and length
    tmp = h.CurrentAxes.XLim;
    x0 = tmp(1);
    DX = tmp(2) - tmp(1);
    tmp = h.CurrentAxes.YLim;
    y0 = tmp(1);
    DY = tmp(2) - tmp(1);
    x1 = x0 + DX*x1_;
    x2 = x0 + DX*x2_;
    y1 = y0 + DY*y1_;
    y2 = y0 + DY*y2_;
    clear tmp
    
    %   Draw line with spec_, num_
%     p = rectangle('Position', [x1, y1, x2 - x1, y2 - y1]);
    p = rectangle('Position', [min(x1, x2), min(y1, y2), max(x1, x2) - min(x1, x2), max(y1, y2) - min(y1, y2)]);
    
    %   Apply nargin
    for i = 1:2:(nargin - 5)
        set(p, varargin{i}, varargin{i + 1});
    end
    hold(h.CurrentAxes, 'off');

end

