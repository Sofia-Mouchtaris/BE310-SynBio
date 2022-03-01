function[r, c, ppm, avg_disk_radius] = find_center(img)
    sz = size(img); % get image dimensions
    % find row and col with the most white pixels
    thresh = 65;
    disk_pix = img >= thresh; % find the white pixels
    disk_center = disk_pix((round(sz(1) / 2) - 200) : (round(sz(1) / 2) + 200), (round(sz(2) / 2) - 200) : (round(sz(2) / 2) + 200));
    disk_rows_sum = sum(disk_center, 2); % sum the number of white pixels in each row
    disk_cols_sum = sum(disk_center, 1); % sum the number of white pixels in each column
    [m, r] = max(disk_rows_sum); % find the row with the most white pixels, restricted to the center of th eplate
    [m, c] = max(disk_cols_sum); % same for columns
    r = r + round(sz(1) / 2) - 200 - 1; % convert to absolute row
    c = c + round(sz(2) / 2) - 200 - 1; % convert to absolute column

    % get line segments going out from center
    left = flip(img(r, 1:c));
    right = img(r, c:end);
    up = flip(img(1:r, c));
    down = img(r:end, c);

    % determine disk radius (pixels) from center of front in each direction
    left_disk = left <= thresh;
    left_disk = find(left_disk, 1);
    right_disk = right <= thresh;
    right_disk = find(right_disk, 1);
    up_disk = up <= thresh;
    up_disk = find(up_disk, 1);
    down_disk = down <= thresh;
    down_disk = find(down_disk, 1);

    % average
    avg_disk_radius = mean([left_disk right_disk up_disk down_disk]);

    % convert pixels to mm
    ppm = avg_disk_radius / 3;
    
end