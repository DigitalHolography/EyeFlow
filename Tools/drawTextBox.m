function out = drawTextBox(img, position, txt, fontSize, boxColor, boxOpacity, textColor)

if ndims(img)==2
    img = repmat(img,1,1,3);
end

boxColor = normalizeColor(boxColor);
textColor = normalizeColor(textColor);

hFig = figure('Visible','off');
imshow(img, 'Border','tight')
axis off
hold on

t = text(position(1), position(2), txt, ...
    'FontSize', fontSize, ...
    'Color', textColor, ...
    'Units','data', ...
    'VerticalAlignment','top');

drawnow

ext = t.Extent;
x = ext(1);
y = ext(2);
w = ext(3);
hBox = ext(4);

patch([x x+w x+w x], [y y y+hBox y+hBox], boxColor, ...
    'EdgeColor','none', ...
    'FaceAlpha', boxOpacity)

t = text(position(1), position(2), txt, ...
    'FontSize', fontSize, ...
    'Color', textColor, ...
    'Units','data', ...
    'VerticalAlignment','top');

drawnow

frame = getframe(gca);
out = frame2im(frame);
close(hFig)

end


function c = normalizeColor(c)
if ischar(c)
    rgb = getColorByName(c);
    c = rgb / 255;
elseif isnumeric(c)
    if max(c) > 1
        c = c / 255;
    end
end
end


function rgb = getColorByName(name)
switch lower(name)
    case 'white',     rgb = [255 255 255];
    case 'black',     rgb = [0 0 0];
    case 'red',       rgb = [255 0 0];
    case 'green',     rgb = [0 255 0];
    case 'blue',      rgb = [0 0 255];
    case 'yellow',    rgb = [255 255 0];
    case 'cyan',      rgb = [0 255 255];
    case 'magenta',   rgb = [255 0 255];
    otherwise
        error('Unknown color name: %s', name)
end
end
