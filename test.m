minX=-2;%exclude
maxX=3;%incldue
minY=-2;%eclude
maxY=3;%include
pattern=ones(11);

[a,b]=ind2sub([maxX-minX,maxY-minY],1:(maxX-minX)*(maxY-minY));
posibleVector=cat(2,a(:)+minX,b(:)+minY);

aoi1f=single(rand(5000,3000));
aoi2f=circshift(aoi1f,[2,3]);


disp("start")
tic;
[lagIndex,quality]=movsae2(single(pattern),aoi1f,aoi2f,int32(posibleVector));
toc