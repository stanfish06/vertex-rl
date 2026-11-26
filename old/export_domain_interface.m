function export_domain_interface(domain_coords, interface_coords, interface_nv, path2initial)
% Compute domain normals and export to files
%   <domain_coords> should be n by 2 array, where the first column is xs and
%   the second column is ys of the domain coordinates.
%
%   <interface_coords> should be 2 by 2 * 3, where each 2 by 2 block
%   represents one interface line. The first and second columns of each 2
%   by 2 block are xs and ys of the two interface end points.
%
%   <interface_nv> should be 3 by 2, where ith row is the normal vector of
%   the ith interface than points into the region closer to the proximal
%   end of the bone (e.g. for interface between resting and proliferative
%   zones, the normal vector should point into the resting zone)
%
%   <path2initial> provide path to the /cartilage/initial folder
    [~, lastpart, ~] = fileparts(path2initial);
    if (lastpart ~= "initial")
        error("File path not ended with initial.");
        return;
    end
    coords = domain_coords;
    if polygon_winding(coords) == -1
        coords = coords(size(coords, 1):-1:1, :);
    end
    v_ns_D = domain_normal_vectors(coords);
    v_ts_D = domain_tangent_vectors(coords);
    writematrix(coords, fullfile(path2initial, 'domain_coords.txt'));
    writematrix(v_ns_D, fullfile(path2initial, 'domain_v_ns.txt'));
    writematrix(v_ts_D, fullfile(path2initial, 'domain_v_ts.txt'));
    writematrix(interface_coords, fullfile(path2initial, 'interface_coords.txt'));
    writematrix(interface_nv, fullfile(path2initial, 'interface_v_ns.txt'));
end

function v_ns = domain_normal_vectors(coords)
    v_ns = zeros(size(coords, 1), 2);
    for i = 1:size(coords, 1)
        p_current = coords(i, :);
        if i == size(coords, 1)
            p_next = coords(1, :);
        else
            p_next = coords(i + 1, :);
        end
        v_t = (p_next - p_current) / norm(p_next - p_current);
        % v_n points into the domain
        v_n = [-v_t(2) v_t(1)];
        v_ns(i, :) = v_n;
    end
end

function v_ts = domain_tangent_vectors(coords)
    v_ts = zeros(size(coords, 1), 2);
    for i = 1:size(coords, 1)
        p_current = coords(i, :);
        if i == size(coords, 1)
            p_next = coords(1, :);
        else
            p_next = coords(i + 1, :);
        end
        v_t = (p_next - p_current) / norm(p_next - p_current);
        v_ts(i, :) = v_t;
    end
end

function direction = polygon_winding(coords)
    area = 0;
    for i = 1:size(coords, 1)
        if i == size(coords, 1)
            coord_pre = coords(i, :);
            coord_next = coords(1, :);
        else
            coord_pre = coords(i, :);
            coord_next = coords(i + 1, :);
        end
        area = area + (coord_next(1) - coord_pre(1)) * (coord_pre(2) + coord_next(2));
    end
    if area < 0
        % CCW
        direction = 1;
    elseif area > 0
        % CW 
        direction = -1;
    else
        direction = 0;
    end
end
