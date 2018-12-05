function AngleVals=image_registr_xcorr(im1, im2,listofangles)
%To register J to im2, use the following code:
%    xoffset=xpeak-size(im1,2);
%    yoffset=ypeak-size(im1,1);
%    TransformOffset=[xoffset,yoffset];
%    transformedim1=imtranslate(imrotate(im1,TransformAngle,'bilinear','crop'),TransformOffset,'OutputView','full');        
%
%Note that this will depend on whether the translation is acting on im1 or
%im2. See FindTranslationXCorr for an example.

%NOTE: The built in function imregcorr might be better.

delete(gcp('nocreate'))
pool=parpool(12);
AngleVals={};
%listofangles=0:0.5:360;
%we could do gradient descent actually, to maximize the maximum normxcorr.
%Especially if you look for the maximum normxcorr close to where you
%previously found the maximum.
parfor angleindex=1:length(listofangles)
    display(['Beginning on angleindex ',num2str(angleindex)])
    angle=listofangles(angleindex);
    J = imrotate(im1, angle,'bilinear','crop'); %rotated cropped IMAGE2
    tmp=normxcorr2(J,im2);
    
    [max_c, imax] = max(abs(tmp(:)));
    [ypeak, xpeak] = ind2sub(size(tmp),imax(1));
%    corr_offset = [(xpeak-size(J,2)) 
%                   (ypeak-size(J,1))];
    
    AngleVals{angleindex}={max_c,[xpeak,ypeak]};
end
    
for k=1:length(AngleVals)
    AngleValValues(k)=cell2mat(AngleVals{k}(1));
end