function outlineMask = drawPolygonOutline(maskSize, polyCoords, lineWidth)

outlineMask = zeros(maskSize);
x = polyCoords(1:2:end);
y = polyCoords(2:2:end);
x(end+1) = x(1);
y(end+1) = y(1);

for k = 1:numel(x)-1
    x1 = x(k);
    y1 = y(k);
    x2 = x(k+1);
    y2 = y(k+1);

    t = linspace(0, 1, max(abs([x2 - x1, y2 - y1])) + 1);
    xs = round(x1 + t * (x2 - x1));
    ys = round(y1 + t * (y2 - y1));

    for n = 1:numel(xs)
        for dx = -floor(lineWidth/2):floor(lineWidth/2)
            for dy = -floor(lineWidth/2):floor(lineWidth/2)
                xx = xs(n) + dx;
                yy = ys(n) + dy;
                if xx >= 1 && xx <= maskSize(2) && yy >= 1 && yy <= maskSize(1)
                    outlineMask(yy, xx) = 1;
                end
            end
        end
    end
end
