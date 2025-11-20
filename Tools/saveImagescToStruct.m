function s = saveImagescToStruct(h)
s = struct();

% Extract image data and properties
s.CData = h.CData;               % actual image data matrix
s.XData = h.XData;               % x-axis limits
s.YData = h.YData;               % y-axis limits
s.CDataMapping = h.CDataMapping; % 'scaled' or 'direct'

% Extract color limits and colormap
ax = ancestor(h, 'axes');
s.XLim = ax.XLim;
s.YLim = ax.YLim;
s.CLim = ax.CLim;
s.Colormap = colormap(ax);

% Optional: store labels and title
s.Title = get(get(ax, 'Title'), 'String');
s.XLabel = get(get(ax, 'XLabel'), 'String');
s.YLabel = get(get(ax, 'YLabel'), 'String');
end