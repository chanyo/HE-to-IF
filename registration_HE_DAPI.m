
clear all;
close all;
clc;

% ORG for H&E


% without ORG -> IF

N_SMPL = 5000;
output_dir = './registered_HE/';
mkdir(output_dir);


dir_root = './4373_D09_Subset_2_TIFF_Images/';
HE_name = sprintf('%s4373_D09_Post-IFH&E_Subset_2_ORG.tif', dir_root);
IF_name = sprintf('%s4373_D09_IF_Subset_2/4373_D09_IF_Subset_2_DAPI_ORG.tif', dir_root);

%dir_root = './4373_D06_TIFF_Images/';
%HE_name = sprintf('%sPost_IF_HE_TIFF_Subsets/4373_D06_PostIFH&E_Subset_3_ORG.tif', dir_root);
%IF_name = sprintf('%sIF_TIFF_Subsets/4373_D06_IF_Subset_3_DAPI_ORG.tif', dir_root);



HE = imread(HE_name);
HE_info = imfinfo(HE_name);
I_DAPI = imread(IF_name);
%I_DAPI = imread(IF_name, 'PixelRegion', {[1 W], [1 H]});
DAPI_info = imfinfo(IF_name);

%%
HE = imresize(HE, DAPI_info.XResolution/HE_info.XResolution);

I_mask = HE(:,:,1);
bw = imbinarize(I_mask, graythresh(I_mask));
I_mask = (bw==0);


I_DAPI = I_DAPI(:,:,1);


%%

if class(I_DAPI) == 'uint16'
    I_ref = uint8((double(I_DAPI)*2^8/2^16)); %imresize(I_DAPI,0.1);
else
    I_ref = I_DAPI;
end

I_ref = imbinarize(I_ref, graythresh(I_ref));
I_obj = I_mask; %imresize(I_mask,0.1);





ptsRef = detectSURFFeatures(I_ref); %SURFFeatures(I_ref);
ptsObj = detectSURFFeatures(I_obj); %SURFFeatures(I_obj);ptsRef


ptsRef = ptsRef.selectStrongest(min(N_SMPL, length(ptsRef)));
ptsObj = ptsObj.selectStrongest(min(N_SMPL, length(ptsObj))); 

[featuresRef, validPtsRef] = extractFeatures( I_ref, ptsRef);
[featuresObj, validPtsObj] = extractFeatures( I_obj, ptsObj);

indxPairs = matchFeatures(featuresRef, featuresObj);

matchedRef = validPtsRef(indxPairs(:,1));
matchedObj = validPtsObj(indxPairs(:,2));



% [tform, inlierDistorted, inlierOriginal,status] = estimateGeometricTransform(...
%     matchedObj, matchedRef, 'affine');

[tform, inlierDistorted, inlierOriginal,status] = estimateGeometricTransform(...
                matchedObj, matchedRef,  'similarity', 'MaxNumTrials',1000, 'Confidence',99.9, 'MaxDistance', 2.5);
               % matchedObj, matchedRef,  'similarity', 'MaxNumTrials',1000, 'Confidence',99.9, 'MaxDistance', 2.5);
                %matchedObj, matchedRef, 'affine');




%
% [optimizer, metric] = imregconfig('multimodal')
%   optimizer.InitialRadius = 0.009;
%     optimizer.Epsilon = 1.5e-4;
%     optimizer.GrowthFactor = 1.01;
%     optimizer.MaximumIterations = 300;
%   tform = imregtform(I_obj, I_ref, 'affine', optimizer, metric);
%    
%    movingRegistered = imwarp(I_obj,tform,'OutputView',imref2d(size(I_ref)));
%    

Tinv  = tform.invert.T;

ss = Tinv(2,1);
sc = Tinv(1,1);
scale_recovered = sqrt(ss*ss + sc*sc)
theta_recovered = atan2(ss,sc)*180/pi



outputView = imref2d( [size(I_ref,1) size(I_ref,2)]);

recovered_obj = imwarp( I_obj, tform,'OutputView', outputView);
%
recovered_HE = uint8(zeros(size(I_ref,1), size(I_ref,2), 3));
for j=1:3
    recovered_HE(:,:,j) = imwarp( HE(:,:,j), tform,'OutputView', outputView);
end
%overlay = imoverlay( I_ref, bwperim(recovered_obj), [1 0 0]);
% imwrite(overlay, 'overlay_registered.png','png');
%
figure
ax(1) = subplot(121); imagesc(I_ref);
ax(2) = subplot(122); imagesc(recovered_HE); %imagesc(recovered_obj);
linkaxes(ax,'xy');
%%

if size(I_ref, 1)*size(I_ref,2)*3 > 2^32-1 
    imwrite(recovered_HE(1:floor(size(I_ref,1)/2),:,:), sprintf('%sregistered_HE_TILE1.tif', output_dir), 'tif')
    imwrite(recovered_HE(floor(size(I_ref,1)/2)+1:end,:,:), sprintf('%sregistered_HE_TILE2.tif', output_dir), 'tif');
else
    imwrite(recovered_HE, sprintf('%sregistered_HE.tif', output_dir),  'tif');
end






