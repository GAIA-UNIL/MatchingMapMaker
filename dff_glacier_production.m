% clear all; clc;
% update teh list of files, that need to have the same resolution
imageName={'file1.tif',  
'file2.tif'
}

% directory with images, need to finish with /
dataDirectory='imageGlacier/'

% patteren size, dosen't change the computation speed
patternSize=51;

% change number of thread to use
maxNumCompThreads(20);

% change the explotation lags
minX=-200;%exclude
maxX=50;%incldue
minY=-200;%exclude
maxY=50;%include

% automatic naming
im1Date=imageName{i}(1:10)
im2Date=imageName{i+1}(1:10)

tic;
[im1,R1]=geotiffread(strcat(dataDirectory,imageName{i}));
[im2,R2]=geotiffread(strcat(dataDirectory,imageName{i+1}));
%%
R1map=[R1.XWorldLimits;R1.YWorldLimits];
R2map=[R2.XWorldLimits;R2.YWorldLimits];
RFinal=[max(R1map(:,1),R2map(:,1)),min(R1map(:,2),R2map(:,2))];

deltaR1=int32(round((RFinal-R1map)./[R1.CellExtentInWorldX; R1.CellExtentInWorldY]));
deltaR2=int32(round((RFinal-R2map)./[R2.CellExtentInWorldX; R2.CellExtentInWorldY]));

stepSize=1;
aoi1=im1(1-deltaR1(2,2):stepSize:end-deltaR1(2,1),1+deltaR1(1,1):stepSize:end+deltaR1(1,2),:);
aoi2=im2(1-deltaR2(2,2):stepSize:end-deltaR2(2,1),1+deltaR2(1,1):stepSize:end+deltaR2(1,2),:);
clear('im1')
clear('im2')

result=mean(aoi1(:,:,1:3),3)-mean(aoi2(:,:,1:3),3);
imagesc(result(1:100:end,1:100:end))

size(aoi1);
size(aoi2);

aoi1f=rgb2gray(single(aoi1(:,:,1:3))/255.);
aoi1f(aoi1(:,:,4)==0)=nan;

aoi2f=rgb2gray(single(aoi2(:,:,1:3))/255.);
aoi2f(aoi2(:,:,4)==0)=nan;

clear('aoi1')
clear('aoi2')

sourceSize=size(aoi1f);
finalSize=sourceSize;

while (max(factor(finalSize(1)))>5)
	finalSize(1)=finalSize(1)+1;
end

while (max(factor(finalSize(2)))>5)
	finalSize(2)=finalSize(2)+1;
end

aoi1f=padarray(aoi1f,finalSize-sourceSize,nan,'post');
aoi2f=padarray(aoi2f,finalSize-sourceSize,nan,'post');


pattern=ones(patternSize);
finalSize

[a,b]=ind2sub([maxX-minX,maxY-minY],1:(maxX-minX)*(maxY-minY));
possibleVector=cat(2,a(:)+minX,b(:)+minY);
possibleVector=possibleVector(randperm(size(possibleVector,1)),:);

tic;
[lagIndex,quality]=mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'SAE');
toc

save(strcat(im1Date,'-',im2Date,'-',int2str(patternSize),'.mat'),'-v7.3')
