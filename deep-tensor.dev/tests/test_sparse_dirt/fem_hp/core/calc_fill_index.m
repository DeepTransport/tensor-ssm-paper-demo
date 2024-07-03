function [subs, vec] = calc_fill_index(node_map)

[num_elem, num_node] = size(node_map);
is = zeros(num_node^2, num_elem);
js = zeros(num_node^2, num_elem);
vec = zeros(num_node, num_elem);
for i = 1:num_elem
    local_ind = node_map(i,:);
    vec(:,i) = local_ind(:);
    local_is = repmat(local_ind', 1, num_node);
    local_js = repmat(local_ind , num_node, 1);
    is(:,i) = local_is(:);
    js(:,i) = local_js(:);
end
subs = [is(:),js(:)];

end
