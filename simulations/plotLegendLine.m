%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Before you use call this function, 
%       DisplayName should be given for p_
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, p] = plotLegendLine(h, x_, y_, p_, spec_, num_, varargin)
%   h: Handle of figure
%   x_: Relative x location of line legend
%   y_: Relative y location of line legend
%   p_: Handle of line object
%   spec_: spec_ of plot_dash
%   num_: num_ of plot_dash
%   varargin: "Pair(s)" of options for text object
%   To use original line style, spec_ = [1, 0], num_ = 1
%   User option: Length of space for line
%     len = 0.1;          %   10 %
    len = 0.2;          %   20 %

    %   Initialization
    hold(h.CurrentAxes, 'on');
    p = plot(0, 0, 'o');
    t = text(0, 0, 'a');
    delete(p);
    delete(t);
    if (p_.DisplayName == "")
        return;
    end

    %   Determine positions and length
    tmp = get(h, 'InnerPosition');
    W = tmp(3) / 560;
    tmp = h.CurrentAxes.XLim;
    x0 = tmp(1);
    DX = tmp(2) - tmp(1);
    tmp = h.CurrentAxes.YLim;
    y0 = tmp(1);
    DY = tmp(2) - tmp(1);
    len = len * DX / W;
    x = x0 + DX*x_;
    y = y0 + DY*y_;
    clear tmp
    
    
    %   Draw line with spec_, num_
    p = plotDash(h, x + [0, len], y*[1, 1], spec_, num_); 
    %   Copy content - color, LineWidth
    p.Color = p_.Color;
    p.LineWidth = p_.LineWidth;
    p.LineStyle = p_.LineStyle;
    
    
    %   Legend text
    t = text(x + 1.1 * len, y, p_.DisplayName);
    %   Apply nargin
    for i = 1:2:(nargin - 6)
        set(t, varargin{i}, varargin{i + 1});
    end
    hold(h.CurrentAxes, 'off');

end

