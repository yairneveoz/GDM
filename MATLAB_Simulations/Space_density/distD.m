function dist_d = distD(d)

d_ones = repmat(d,size(d'));

dist_d = d_ones - d_ones';

end