function [diff_path, diff_length] = getDiffPath(path1)

diff_path = diff(path1,1,1);
diff_length = sqrt(sum(diff_path.^2,2));

end