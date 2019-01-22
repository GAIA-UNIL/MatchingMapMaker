minX=-2;%exclude
maxX=3;%incldue
minY=-2;%eclude
maxY=3;%include
pattern=ones(11);

[a,b]=ind2sub([maxX-minX,maxY-minY],1:(maxX-minX)*(maxY-minY));
possibleVector=cat(2,a(:)+minX,b(:)+minY);

aoi1f=single(rand(5000,3000));
%single(filter2(exp(-0.2*bwdist(padarray(1,[50 50]))),randn(750,500)));%
aoi2f=circshift(aoi1f,[2,3]);


disp("start")
tic;
[ lagIndex_SAE, quality_SAE ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'SAE'); %SAE
[ lagIndex_MAE, quality_MAE ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'MAE'); %MAE
[ lagIndex_SSE, quality_SSE ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'SSE'); %SSE
[ lagIndex_MSE, quality_MSE ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'MSE'); %MSE
[ lagIndex_CC, quality_CC ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'CC'); %CC
[ lagIndex_NMAE, quality_NMAE ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'NMAE'); %NMAE
[ lagIndex_NMSE, quality_NMSE ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'NMSE'); %NMSE
[ lagIndex_NCC, quality_NCC ] = mmm(single(pattern),aoi1f,aoi2f,int32(possibleVector),'NCC'); %NCC
toc

fprintf("SAE %f %%\n",sum(lagIndex_SAE(:)==24)/numel(lagIndex_SAE)*100);
fprintf("MAE %f %%\n",sum(lagIndex_MAE(:)==24)/numel(lagIndex_MAE)*100);
fprintf("SSE %f %%\n",sum(lagIndex_SSE(:)==24)/numel(lagIndex_SSE)*100);
fprintf("MSE %f %%\n",sum(lagIndex_MSE(:)==24)/numel(lagIndex_MSE)*100);
fprintf("CC %f %%\n",sum(lagIndex_CC(:)==24)/numel(lagIndex_CC)*100);
fprintf("NMAE %f %%\n",sum(lagIndex_NMAE(:)==24)/numel(lagIndex_NMAE)*100);
fprintf("NMSE %f %%\n",sum(lagIndex_NMSE(:)==24)/numel(lagIndex_NMSE)*100);
fprintf("NCC %f %%\n",sum(lagIndex_NCC(:)==24)/numel(lagIndex_NCC)*100);
