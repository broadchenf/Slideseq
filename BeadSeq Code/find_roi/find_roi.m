% Locates a region of interest within a map. Uses the scale invariant
% feature transform (SIFT) with random sample consensus (RANSAC) to fit
% either a geometric mapping with rotation, translation, and stretch. SIFT
% keypoints are computed using the VLFeat package, available at
% www.vlfeat.org. For more information about SIFT, see D. G. Lowe, Int. J.
% Comput. Vision 60(2), 91--110, 2004.

close all;
clearvars -except f1s d1s database;
run('helpers\vlfeat-0.9.20\toolbox\vl_setup.m');
addpath helpers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
tileSize = 3500; %Decrease if out of memory
peakThresh = 0; %SIFT peak threshold
edgeThresh = 10; %SIFT edge threshold
nRansacTrials = 2000; %Increase for more reliable matching
nPtsFit = 2; %For each RANSAC trial
nKeypointsThresh = 50; %Minimum number of keypoint matches
radius = 20; %Maximum tolerated RANSAC distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load image data
%[mapName,queryName] = get_names();
mapName='testimage1.tiff'
queryName='testimage2.tiff'
mapOrig = imread(mapName);
map = imadjust(im2single(mapOrig));
query = im2single(imadjust(imread(queryName)));
disp(['Loaded map image ',mapName]);
disp(['Loaded query image ',queryName]);
if exist('database','var')
    db_old = database; %Store copy of existing map tiles
else
    db_old = {};
end

% Generate database of tiles
[m,n] = size(map);
N = ceil(n/tileSize);
M = ceil(m/tileSize);
nTiles = M*N;
database = cell(1,nTiles);
for i = 1:M
    for j = 1:N
        img = im2single(map((i-1)*floor(m/M)+1:i*min(floor(m/M),m),...
            (j-1)*floor(n/N)+1:j*min(floor(n/N),n)));
        database{(i-1)*N+j} = img/max(img(:)); %Normalize intensity
    end
end

% Check if SIFT keypoints already exist in workspace
isPrecomputed = exist('f1s','var') && exist('d1s','var') ...
    && exist('database','var') && iscell(f1s) && iscell(d1s) ...
    && iscell(database) && (numel(f1s) == numel(database)) && ...
    (numel(d1s) == numel(database)) && isequal(database,db_old);

% Compute SIFT keypoints for map
if isPrecomputed
    disp('Existing SIFT keypoints found in workspace.');
else
    f1s = cell(1,nTiles);
    d1s = cell(1,nTiles);
    disp('Computing SIFT keypoints...');
    reverseStr = '';
    for i = 1:nTiles
        [f1s{i},d1s{i}] = vl_sift(database{i},...
            'PeakThresh',peakThresh,'EdgeThresh',edgeThresh);
        msg = sprintf('Completed tile %d of %d\n',i,nTiles);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
    end
end

% SIFT keypoints for query image
[f2,d2] = vl_sift(query,'PeakThresh',peakThresh,...
    'EdgeThresh',edgeThresh);

% Find matching image
nMatches = zeros(1,nTiles);
matchResult = cell(1,nTiles);
disp('Matching with isometry model...');
reverseStr = '';
for i = 1:nTiles
    % SIFT matches
    tile = database{i};
    f1 = f1s{i};
    d1 = d1s{i};
    [matches, ~] = vl_ubcmatch(d1,d2); % Distance ratio test
    
    % Remove many-to-one matches
    [uniqueRow2, IA, ~] = unique(matches(2,:));
    uniqueRow1 = matches(1,IA);
    matches = [uniqueRow1; uniqueRow2];
    numMatches = size(matches,2);
    
    X1 = f1(1:2,matches(1,:)); X1(3,:) = 1;
    X2 = f2(1:2,matches(2,:)); X2(3,:) = 1;
    
    % RANSAC with geometric model
    H = cell(1,nRansacTrials);
    isMatch = cell(1,nRansacTrials);
    score = zeros(1,nRansacTrials);
    matched_mse = zeros(1,nRansacTrials);
    for j = 1:nRansacTrials
        % Estimate isometry
        subset = vl_colsubset(1:numMatches,nPtsFit);
        X = X1(1:2,subset);
        Y = X2(1:2,subset);
        [Q,s,t] = fit_isometry(X,Y);
        H{j} = [cat(2,s*Q,t); 0 0 1];
        % Score isometry
        X2_ = H{j}*X1;
        du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:);
        dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:);
        isMatch{j} = (du.^2 + dv.^2 < radius^2);
        score(j) = sum(isMatch{j});
        matched_mse(j) = sum((du.^2 + dv.^2).*isMatch{j});
    end
    
    % Find best mapping for current tile
    [~,best] = find(score == max(score));
    [~,idx] = min(matched_mse(best));
    H = H{best(idx)};
    isMatch = isMatch{best(idx)};
    msg = sprintf('Completed tile %d of %d\n',i,nTiles);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
    
    % Plot
    figure(1); subplot(M,N,i); hold on;
    dh1 = max(size(query,1)-size(tile,1),0);
    dh2 = max(size(tile,1)-size(query,1),0);
    imshow([padarray(tile,dh1,'post') padarray(query,dh2,'post')],...
        'InitialMagnification','fit');
    title(sprintf('Tile %d of %d\n%d keypoints matched',...
        i,nTiles,sum(isMatch)));
    o = size(tile,2);
    hold on;
    plot([f1(1,matches(1,isMatch))+1;f2(1,matches(2,isMatch))+o+1], ...
        [f1(2,matches(1,isMatch))+1;f2(2,matches(2,isMatch))+1],...
        ['y' '-'],'LineWidth',1);
    plot(f1(1,matches(1,isMatch))+1,f1(2,matches(1,isMatch))+1,...
        'yo','MarkerSize',5,'MarkerFaceColor','y');
    plot(f2(1,matches(2,isMatch))+o+1,f2(2,matches(2,isMatch))+1,...
        'ko','MarkerSize',5,'MarkerFaceColor','y');
    hold off;
    
    % Update
    matchResult{i}.ratio_test = size(matches,2);
    matchResult{i}.ransac = sum(isMatch);
    matchResult{i}.mse = matched_mse(best(idx));
    matchResult{i}.model = H;
    matchResult{i}.matches = matches(:,isMatch);
    matchResult{i}.f1 = f1;
    matchResult{i}.f2 = f2;
    matchResult{i}.loc1 = f1(1:2,matches(1,isMatch));
    matchResult{i}.loc2 = f2(1:2,matches(2,isMatch));
    
    nMatches(i) = sum(isMatch);
    if nMatches(i) >= nKeypointsThresh
        disp('Match found!');
        break;
    end
end

% Find best-matched tile
[~,idx] = max(nMatches);
matchResult = matchResult{idx};
fprintf('Best match: Tile %d of %d with %d keypoints matched\n\n',...
    idx,nTiles,max(nMatches));

% Calculate image borders using homography
H = inv(matchResult.model);
[h,w] = size(query);
x = [(1:w)'; (1:w)'; ones(h-2,1); w*ones(h-2,1)];
y = [ones(w,1); h*ones(w,1); (2:h-1)'; (2:h-1)'];
z = H(3,1)*x + H(3,2)*y + H(3,3);
x_ = (H(1,1)*x + H(1,2)*y + H(1,3))./z;
y_ = (H(2,1)*x + H(2,2)*y + H(2,3))./z;
loc = [x_,y_]';

% Display
disp(matchResult);
nLoc = size(loc,2);
mask = zeros(size(map));
i = ceil(idx/M);
j = idx - (i-1)*N;
x_offset = (i-1)*floor(m/M);
y_offset = (j-1)*floor(n/M);
y = round(loc(1,:)) + y_offset*ones(1,nLoc);
x = round(loc(2,:)) + x_offset*ones(1,nLoc);
for k = 1:nLoc
    if 1 <= y(k) && y(k) <= n && 1 <= x(k) && x(k) <= m
        mask(x(k),y(k)) = 1;
    end
end
warning('off','Images:initSize:adjustingMag');
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(121); imshow(imadjust(im2double(mapOrig) + ...
    imdilate(mask,strel('disk',8))));
subplot(122); imshow(imtransform(query,maketform('affine',H')));
warning('on','Images:initSize:adjustingMag');

% Crop
ly = find(sum(mask) > 0,1,'first');
ry = find(sum(mask) > 0,1,'last');
ux = find(sum(mask,2) > 0,1,'first');
lx = find(sum(mask,2) > 0,1,'last');
cropped_roi = mapOrig(ux:lx,ly:ry);
mask = mask(ux:lx,ly:ry);

H = [cat(2,normc(matchResult.model(1:2,1:2))',[0;0]); 0 0 1];
tform = maketform('affine',H);
[cropped_roi] = imtransform(cropped_roi,tform);
mask = imtransform(mask,tform);
ly = find(sum(mask) > 0,1,'first');
ry = find(sum(mask) > 0,1,'last');
ux = find(sum(mask,2) > 0,1,'first');
lx = find(sum(mask,2) > 0,1,'last');
%cropped_roi = cropped_roi(ux:lx,ly:ry);
imwrite(cropped_roi,'result.tif');

% clearvars -except mapOrig map query f1s d1s f2 d2 database ...
%     isPrecomputed matchResult mapName queryName ...
%     peakThresh edgeThresh tileSize radius nKeypointsThresh ...
%     n m N M nTiles nMatches idx cropped_roi;