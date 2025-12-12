function out = drawTextBox(img, position, txt, fontSize, boxColor, boxOpacity, textColor)

if ndims(img)==2
    img = repmat(img,1,1,3);
end

boxColor = normalizeColor(boxColor);
textColor = normalizeColor(textColor);

hFig = figure('Visible','off');
% Ensure the figure size matches the image to avoid resizing artifacts during getframe
hFig.Position(3:4) = [size(img,2), size(img,1)]; 

imshow(img, 'Border','tight')
axis off
hold on

% 1. Draw text invisibly first to calculate extent
t = text(position(1), position(2), txt, ...
    'FontSize', fontSize, ...
    'Color', textColor, ...
    'Units','data', ...
    'VerticalAlignment','top', ...
    'Visible', 'on'); % Must be visible for Extent to be accurate in some versions

drawnow

% 2. Get Extent
ext = t.Extent;
x = ext(1);
w = ext(3);
hBox = ext(4);

% In imshow (YDir='reverse'), ext(2) is the visual bottom (max Y value).
% We need the visual top (min Y value) for the patch start.
y = ext(2) - hBox; 

% Optional: Add a small margin so the text isn't flush against the box edge
margin = 2; 
x = x - margin;
y = y - margin;
w = w + 2*margin;
hBox = hBox + 2*margin;

% 3. Draw the background box
patch([x x+w x+w x], [y y y+hBox y+hBox], boxColor, ...
    'EdgeColor','none', ...
    'FaceAlpha', boxOpacity)

% 4. Redraw text on top of the box
% (We can actually just reuse 't' by bringing it to front, strictly speaking, 
% but deleting and recreating ensures layer order)
delete(t);
text(position(1), position(2), txt, ...
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
