
% test_rings
grid_size = 41;
center_array = zeros(grid_size);
center_array(ceil(grid_size/2), ceil(grid_size/2)) = 1;
h_in = ceil(fspecial('disk',2));
h_out = ceil(fspecial('disk',3));
array_in = conv2(center_array, h_in, 'same');
array_out = conv2(center_array, h_out, 'same');
array_ring = array_out - array_in;
for rr = 0:floor(grid_size/2)
    if rr == 0
        array_ring = center_array;
    else
        if rr == 1
            h_in = 1;
            h_out = ceil(fspecial('disk',rr));
            array_in = conv2(center_array, h_in, 'same');
            array_out = conv2(center_array, h_out, 'same');
            array_ring = array_out - array_in;
        else
            h_in = ceil(fspecial('disk',rr-1));
            h_out = ceil(fspecial('disk',rr));
            array_in = conv2(center_array, h_in, 'same');
            array_out = conv2(center_array, h_out, 'same');
            array_ring = array_out - array_in;
        end
    end
    figure(38)
    clf
    spy(array_ring)
    pause(0.5)
    drawnow
end


