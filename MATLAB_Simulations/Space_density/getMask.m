function normalized_mask = getMask(mask_size)

[mask_j, mask_i] = meshgrid(1:mask_size);
mask_center = ceil([mask_size/2, mask_size/2]);
mask_rs = sqrt((mask_i-mask_center(1)).^2 +...
    (mask_j-mask_center(2)).^2);
% mask_array = 0.5*(1 + 1./mask_rs.^2);
mask_array = 0.25./mask_rs.^2;
mask_array(mask_center, mask_center) = 1;
normalized_mask = mask_array/sum(mask_array(:));

end