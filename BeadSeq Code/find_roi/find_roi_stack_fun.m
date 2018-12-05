function find_roi_stack_fun(BaseName,suffix,ImageSize,varargin)
% Locates a region of interest within a map. Uses the scale invariant
% feature transform (SIFT) with random sample consensus (RANSAC) to fit
% either a geometric mapping with rotation, translation, and stretch. SIFT
% keypoints are computed using the VLFeat package, available at
% www.vlfeat.org. For more information about SIFT, see D. G. Lowe, Int. J.
% Comput. Vision 60(2), 91--110, 2004.

%Compared to 170818, I am changing parameters to try to get a better fit.



delete(gcp('nocreate'));
clearvars -except f1s d1s database2 BaseName suffix ImageSize varargin;
run('C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\helpers\vlfeat-0.9.20\toolbox\vl_setup.m');
addpath('C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\helpers');

displayfigs=0;

RegisterColorChannels=0;
index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="RegisterColorChannels"), varargin, 'UniformOutput', 1));
if ~isempty(index)
    RegisterColorChannels=varargin{index+1};
end

XCorrBounds=[1,ImageSize,1,ImageSize];
index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="XCorrBounds"), varargin, 'UniformOutput', 1));
if ~isempty(index)
    XCorrBounds=varargin{index+1};
end

PixelCutoff=300;
index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="PixelCutoff"), varargin, 'UniformOutput', 1));
if ~isempty(index)
    PixelCutoff=varargin{index+1};
end

channelnum=4; %number of channels
index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="NumChannels"), varargin, 'UniformOutput', 1));
if ~isempty(index)
    channelnum=varargin{index+1};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
tileSize = 10500; %Decrease if out of memory. Default was 3500. One problem we're having possibly is that when the images are translate too far down in Y, then there aren't enough keypoints in the upper left quadrant to make the match.
peakThresh = 0; %SIFT peak threshold
edgeThresh = 10; %SIFT edge threshold. default 10
nRansacTrials = 100000; %2000 by default. Increase for more reliable matching
nPtsFit = 2; %For each RANSAC trial
nKeypointsThresh = 30; %Minimum number of keypoint matches. Default was 50. We use a very low threshold because the match is essentially always on tile 1.
radius = 5; %Maximum tolerated RANSAC distance. 20 by default.
MatchThresh=1.5; %default is 1.5. Two keypoints d1 and d2 are matched only if the distance between them times this number is not greater than the distance between d1 and all other keypoints.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTE: For a typical nimage on which the algorithm works, it will find 12000
%features in one tile of the master image of size 3500; and it will find 400000 features in the
%whole query image.


l=0;
while true
    if exist([BaseName,pad(num2str(l+1),2,'left','0'),suffix,'.tif'],'file')
        l=l+1;
    else
        break
    end
end


for k=1:channelnum
    mapOrigstack(:,:,k)=imread([BaseName,pad(num2str(1),2,'left','0'),suffix,'.tif'],'index',k,'PixelRegion',{[1,ImageSize],[1,ImageSize]});
end

if RegisterColorChannels
    mapOrigstack=uint16(FindTranslationXCorr(mapOrigstack,'XCorrBounds',XCorrBounds,'PixelCutoff',PixelCutoff));
end

for k=1:channelnum
    mapstack(:,:,k) = im2single(imadjust(mapOrigstack(:,:,k)));
end

map=max(mapstack,[],3);

for k=1:channelnum
    imwrite(mapOrigstack(:,:,k),[BaseName,'01',' channel ',num2str(k),suffix,' transform.tif'])
end

for mm=2:l
    display(strcat('Loading file for ligation ',num2str(mm)))
    for k=1:channelnum
        queryOrigstack(:,:,k,mm)=imread([BaseName,pad(num2str(mm),2,'left','0'),suffix,'.tif'],'index',k,'PixelRegion',{[1,ImageSize],[1,ImageSize]});
    end
%    queryOrigstack1(:,:,mm)=imread([BaseName,pad(num2str(mm),2,'left','0'),suffix,'.tif'],'index',1,'PixelRegion',{[1,ImageSize],[1,ImageSize]});    %this is clumsy, but makes it easier for matlab to figure out which parts of the array it needs to send to each core in the parfor loop.
%    queryOrigstack2(:,:,mm)=imread([BaseName,pad(num2str(mm),2,'left','0'),suffix,'.tif'],'index',2,'PixelRegion',{[1,ImageSize],[1,ImageSize]});    
%    queryOrigstack3(:,:,mm)=imread([BaseName,pad(num2str(mm),2,'left','0'),suffix,'.tif'],'index',3,'PixelRegion',{[1,ImageSize],[1,ImageSize]});
%    queryOrigstack4(:,:,mm)=imread([BaseName,pad(num2str(mm),2,'left','0'),suffix,'.tif'],'index',4,'PixelRegion',{[1,ImageSize],[1,ImageSize]});
end
if RegisterColorChannels
    for mm=2:l
        display(['Registering color channels for ligation ',num2str(mm)])
        queryOrigstack(:,:,:,mm)=uint16(FindTranslationXCorr(squeeze(queryOrigstack(:,:,:,mm)),'XCorrBounds',XCorrBounds,'PixelCutoff',PixelCutoff));
    end
end

for mm=2:l
    for k=1:channelnum
        querystack(:,:,k,mm) = im2single(imadjust(queryOrigstack(:,:,k,mm)));
    end
end
pool=parpool(20);
addAttachedFiles(pool,'C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\helpers\fit_isometry.m');

parfor mm=2:l
%for mm=12
    %queryOrig=max(queryOrigstack(:,:,:,mm),[],3);
    query=max(querystack(:,:,:,mm),[],3);

    
    %When isPrecomputed is not in use, db_old is not used either
    if 0
        if exist('database2','var')
            db_old = database2; %Store copy of existing map tiles
        else
            db_old = {};
        end
    end    
    % Generate database of tiles
    [m,n] = size(map);
    N = ceil(n/tileSize);
    M = ceil(m/tileSize);
    nTiles = M*N;
    database2 = cell(1,nTiles);
    for i = 1:M
        for j = 1:N
            img = im2single(map((i-1)*floor(m/M)+1:i*min(floor(m/M),m),...
                (j-1)*floor(n/N)+1:j*min(floor(n/N),n)));
            database2{(i-1)*N+j} = img/max(img(:)); %Normalize intensity
        end
    end

    % Check if SIFT keypoints already exist in workspace
    if 0 %The SIFT keypoints are never precomputed in the parfor
    isPrecomputed = exist('f1s','var') && exist('d1s','var') ...
        && exist('database2','var') && iscell(f1s) && iscell(d1s) ...
        && iscell(database2) && (numel(f1s) == numel(database2)) && ...
        (numel(d1s) == numel(database2)) && isequal(database2,db_old);
    else
        isPrecomputed=0;
    end

    % Compute SIFT keypoints for map
    if isPrecomputed
        disp('Existing SIFT keypoints found in workspace.');
    else
        f1s = cell(1,nTiles);
        d1s = cell(1,nTiles);
        disp('Computing SIFT keypoints...');
        reverseStr = '';
        for i = 1:nTiles
            [f1s{i},d1s{i}] = vl_sift(database2{i},...
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
        tile = database2{i};
        f1 = f1s{i};
        d1 = d1s{i};
        [matches, ~] = vl_ubcmatch(d1,d2,MatchThresh); % Distance ratio test. This takes most of the time.

        % Remove many-to-one matches
        [uniqueRow2, IA, ~] = unique(matches(2,:));
        uniqueRow1 = matches(1,IA);
        matches = [uniqueRow1; uniqueRow2]; %the number of matches we get should be close to the number of fixed points in the tile of the map image.
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
            dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:); %The model we use here for the RANSAC is that the images are related to each other by an affine transformation.
            isMatch{j} = (du.^2 + dv.^2 < radius^2);
            score(j) = sum(isMatch{j});
            matched_mse(j) = sum((du.^2 + dv.^2).*isMatch{j});
        end

        % Find best mapping for current tile
        [~,best] = find(score == max(score)); %A good score here is like 1574, for example. By contrast, in the failing images we get scores around 5.
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
    %disp(matchResult);
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

    if displayfigs
        figure('units','normalized','outerposition',[0 0 1 1]);

        subplot(121); imshow(imadjust(im2double(mapOrig) + ...
        imdilate(mask,strel('disk',8))));
        subplot(122); imshow(imtransform(query,maketform('affine',H')));
        warning('on','Images:initSize:adjustingMag');

        testmask1=mask;
        testmask2=imdilate(mask,strel('disk',8));

        % Crop
        ly = find(sum(mask) > 0,1,'first');
        ry = find(sum(mask) > 0,1,'last');
        ux = find(sum(mask,2) > 0,1,'first');
        lx = find(sum(mask,2) > 0,1,'last');
        cropped_roi = mapOrig(ux:lx,ly:ry);
        mask = mask(ux:lx,ly:ry);

    end
    

%In the original code, linyi was using the forward transformation here on cropped_roi, which is derived from map, whereas we are using the inverse transformation on query    
%original:
%    H = [cat(2,normc(matchResult.model(1:2,1:2))',[0;0]); 0 0 1];
%    H = [cat(2,normc(matchResult.model(1:2,1:2))',[0;0]); 0 0 1]';
%    H = [cat(2,normc(matchResult.model(1:2,1:2))',[0;0]); cat(2,matchResult.model(1:2,3)', [1])];
%    H = [normc(matchResult.model(1:2,1:3)); 0 0 1]';
    H=inv(matchResult.model);
    tform = maketform('affine',H');

    if displayfigs
        [cropped_roi] = imtransform(cropped_roi,tform);
        mask = imtransform(mask,tform);
        ly = find(sum(mask) > 0,1,'first');
        ry = find(sum(mask) > 0,1,'last');
        ux = find(sum(mask,2) > 0,1,'first');
        lx = find(sum(mask,2) > 0,1,'last');
        figure(6)
        clf
        subplot(221)
        imshow(imtransform(query,tform,'XData',[1, size(query,2)],'YData',[1, size(query,1)])+testmask2);
        subplot(222)
        imshow(map)
    end
    
    %just debugging code
%    translationmatrix=eye(3);
%    translationmatrix(3,1)=100;
%    translationtform=maketform('affine',translationmatrix);
%    figure(10)
%    subplot(2,2,1)
%    imshow(query);
%    subplot(2,2,2)
%    imshow(imtransform(query,translationtform));
%    subplot(2,2,3)
%    imshow(imtransform(query,translationtform,'XData',[1 querysize(2)],'YData',[1,querysize(1)]));
    

%    figure(6)
%    transformedquery=imtransform(query,tform);
%    Hbar=matchResult.model';
    %hbarinv=inv(Hbar)';
    %translationamount=hbarinv(3,1:2);
%    translationamount=size(query)-Hbar(3,1:2);
%    translatedquery=imtranslate(transformedquery,translationamount,'outputView','full');
%    subplot(2,2,3)
%    imshow(translatedquery)
    %    t=1;
%    for t=1:querysize(1)
%        if t+ytranslate<=querysize(1) && t+ytranslate>0
%            if xtranslate>0
%                translatedquery(t+ytranslate,:)=[zeros(1,xtranslate),transformedquery(t,max(1,1+xtranslate):min(querysize(2),querysize(2)+xtranslate))];
%            elseif xtranslate<0
%                translatedquery(t+ytranslate,:)=[transformedquery(t,max(1,1+xtranslate):min(querysize(2),querysize(2)+xtranslate)),zeros(1,-xtranslate)];
%            end            
%        end
%    end
    %we do this rather than a for loop, because it avoids matlab trying to
    %end the entire queryorigstack variable to all the cores on the parfor.
    try
        for k=1:4
            imwrite(imtransform(squeeze(queryOrigstack(:,:,k,mm)),tform,'XData',[1, size(query,2)],'YData',[1, size(query,1)]),[BaseName,pad(num2str(mm),2,'left','0'),' channel ',int2str(k),suffix,' transform.tif'])
        end
%       k=1
%        imwrite(imtransform(queryOrigstack1(:,:,mm),tform,'XData',[1, size(query,2)],'YData',[1, size(query,1)]),[BaseName,pad(num2str(mm),2,'left','0'),' channel ',int2str(k),suffix,' transform.tif'])
%        k=2;
%        imwrite(imtransform(queryOrigstack2(:,:,mm),tform,'XData',[1, size(query,2)],'YData',[1, size(query,1)]),[BaseName,pad(num2str(mm),2,'left','0'),' channel ',int2str(k),suffix,' transform.tif'])
%        k=3;
%        imwrite(imtransform(queryOrigstack3(:,:,mm),tform,'XData',[1, size(query,2)],'YData',[1, size(query,1)]),[BaseName,pad(num2str(mm),2,'left','0'),' channel ',int2str(k),suffix,' transform.tif'])
%        k=4;
%        imwrite(imtransform(queryOrigstack4(:,:,mm),tform,'XData',[1, size(query,2)],'YData',[1, size(query,1)]),[BaseName,pad(num2str(mm),2,'left','0'),' channel ',int2str(k),suffix,' transform.tif'])
    catch ME
        disp(['The following error was caught while writing out the transformed images for mm=',num2str(mm),' and BaseName=',BaseName,':'])
        disp([ME.identifier,': ',ME.message])
        disp('The transform is:')
        disp(tform.tdata)
    end

%    [transformedquery,xmove,ymove]=imtransform(query, tform);
    %NOTE: We are going to try only using the first entry of xmove and
    %ymove, which is the location of the first column in each dimension

%     outputstack((dims(1)+1+xmove(1)):(2*dims(1)+xmove(1)),(dims(2)+1+ymove(1)):(2*dims(2)+ymove(1)),:,mm)=transformedquery
    
%    for k=1:channelnum
%        outputstack((dims(1)+1+xmove(1)):(2*dims(1)+xmove(1)),(dims(2)+1+ymove(1)):(2*dims(2)+ymove(1)),:,mm)=transformedquery;
%    end



%cropped_roi = cropped_roi(ux:lx,ly:ry);
%    imwrite(cropped_roi,'result.tif');
%    for k=1:channelnum
%        alignedimage(:,:,k,mm)=
%    end
    % clearvars -except mapOrig map query f1s d1s f2 d2 database ...
    %     isPrecomputed matchResult mapName queryName ...
    %     peakThresh edgeThresh tileSize radius nKeypointsThresh ...
    %     n m N M nTiles nMatches idx cropped_roi;
end
delete(pool);
end

