function idx = nms_tlwh(boxes_tlwh, scores, thresh)

boxes_xyxy = [boxes_tlwh(:,1), ...
              boxes_tlwh(:,2), ...
              boxes_tlwh(:,1)+boxes_tlwh(:,3), ...
              boxes_tlwh(:,2)+boxes_tlwh(:,4)];

[~, order] = sort(scores, 'descend');
idx = [];

while ~isempty(order)
    i = order(1);
    idx(end+1) = i;

    if numel(order) == 1
        break
    end

    rest = order(2:end);

    xx1 = max(boxes_xyxy(i,1), boxes_xyxy(rest,1));
    yy1 = max(boxes_xyxy(i,2), boxes_xyxy(rest,2));
    xx2 = min(boxes_xyxy(i,3), boxes_xyxy(rest,3));
    yy2 = min(boxes_xyxy(i,4), boxes_xyxy(rest,4));

    w = max(0, xx2 - xx1);
    h = max(0, yy2 - yy1);
    inter = w .* h;

    area_i = (boxes_xyxy(i,3)-boxes_xyxy(i,1)) * ...
             (boxes_xyxy(i,4)-boxes_xyxy(i,2));
    area_rest = (boxes_xyxy(rest,3)-boxes_xyxy(rest,1)) .* ...
                (boxes_xyxy(rest,4)-boxes_xyxy(rest,2));

    union = area_i + area_rest - inter;
    iou = inter ./ union;

    % ALWAYS remove i → guaranteed termination
    order = rest(iou < thresh);
end
