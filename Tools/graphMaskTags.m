function himage = graphMaskTags(figId, Image, mask, etiquettes_locs, etiquettes_values, x_center, y_center, cmap, opt)
% Plots on an existing fig the image combination of a raw image and a mask displayed in red
arguments
    figId
    Image
    mask
    etiquettes_locs
    etiquettes_values
    x_center
    y_center
    cmap
    opt.Fontsize double = 14
    opt.Title = []
    opt.Visible = false
    opt.circles = [];

end

ratio_etiquette = 1.2;


% image_RGB = repmat(Image - Image .* mask, 1, 1, 3) + reshape(NameValueArgs.Color, 1, 1, 3) .* mask .* Image; % adding the Red value to the mask pixels
% Get the image dimensions

[imgHeight, imgWidth, ~] = size(Image);
t_matrix = zeros(1, size(etiquettes_locs, 1));

if ~isempty(etiquettes_locs)

    for etIdx = 1:size(etiquettes_locs, 1)

        try
            % Calculate the new position for the text
            new_x = x_center + ratio_etiquette * (etiquettes_locs(etIdx, 1) - x_center);
            new_y = y_center + ratio_etiquette * (etiquettes_locs(etIdx, 2) - y_center);

            % Ensure the text is within the image bounds
            new_x = max(1, min(new_x, imgWidth)); % Clamp x to image width
            new_y = max(1, min(new_y, imgHeight)); % Clamp y to image height

            % Add the text
            t = text(new_x, new_y, sprintf(string(etiquettes_values(etIdx))), ...
                "FontWeight", "bold", ...
                "FontSize", opt.Fontsize, ...
                "Color", "white", ...
                "BackgroundColor", "black");

            t_matrix(etIdx) = t;
        catch
        end

    end

end

end
