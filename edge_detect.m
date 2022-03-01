addpath(genpath('plate_images'))

% load in images

strain = 1;
rep = 1;

images = {dir(sprintf('plate_images/strain%d_rep%d/1.0s_Exposure/*.png', strain, rep)).name};
images = cellfun(@(x) sprintf('plate_images/strain%d_rep%d/1.0s_Exposure/%s', strain, rep, x),images,'UniformOutput',false);

% locate image center
last = images{end};
images = [last, images(1:end-1)];
edge = zeros(length(images), 2);

imgStart = rgb2gray(imread(images{1}));

[r, c, ppm, rad] = find_center(imgStart);

imshow(imgStart)
vline(c - rad, 'b')
vline(c + rad, 'b')
vline(c)

edge(1,:) = [round(c - rad), round(c + rad)];

%%

for i_img = 2:length(images)

img2 = rgb2gray(imread(images{i_img}));

imgdiff = img2 - imgStart;

% get middle row, representative of entire plate
middle = double(imgdiff(r,:));

% calculate baseline outside of previous step edge
meanVal = mean([middle(1:edge(i_img-1,1)), middle(edge(i_img-1,2):end)]);
stdVal = std([middle(1:edge(i_img-1,1)), middle(edge(i_img-1,2):end)]);

% smooth data
avgd = movmean(middle,30);

% find first time data dips below threshold --> new edge
right = c - 1 + find(middle(c:end) <= meanVal + .18*stdVal, 1, 'first');
left = find(middle(1:c) <= meanVal + .18*stdVal, 1, 'last');

if right < edge(i_img-1,2) && left > edge(i_img-1,1)
    edge(i_img,:) = edge(i_img-1, :);
else
    edge(i_img,:) = [right, left];
end

% plot plate
% imshow(img2)
% vline(edge(i_img,1))
% vline(edge(i_img,2))
% pause

end

edge(:,2) = abs(edge(:,2) - c); % left
edge(:,1) = abs(edge(:,1) - c); % right

figure(2)
plot(edge./ppm)
hline(40)

fin = [0:length(images)-1; (edge(:,1)./ppm)'; (edge(:,2)./ppm)']';