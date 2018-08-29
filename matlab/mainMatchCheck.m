more off
% workOutDir = '../images/GRAF/OutMatch/';
% imgDir      = '../images/GRAF/org/';
workDir      = '../images/Boat/';
workOutDir   = '../images/Boat/OutMatch/';

tmplImgName = 'img1'; % Not dat but image for display
inpImgName  = 'img6';
stOrd       = 1;      % Start to display match point in order.
enOrd       = 4;      % End to display match point in order.

datFileTmpl = strcat(workOutDir, tmplImgName, '.template.dat');
datFile     = strcat(workOutDir, inpImgName, '.dat');

[featureLTmpl infVecLTmpl nxTmpl nyTmpl] = readFeatures(datFileTmpl);
[featureL     infVecL     nx     ny   ]  = readFeatures(datFile);

nDataTmpl = size(infVecLTmpl, 2);
nData     = size(infVecL, 2);

matchDegInf  = zeros(size(infVecL, 1) + 2, nDataTmpl * nData);
matchDegInf(1, :) = reshape(featureLTmpl' * featureL, 1, nDataTmpl * nData);
matchDegInf(2, :) = reshape((1:nDataTmpl)' * ones(1, nData), 1, nDataTmpl * nData);

for pos = 1:nData
  matchDegInf(3:size(matchDegInf, 1), ((pos - 1) * nDataTmpl + 1):(pos * nDataTmpl)) = infVecL(:,pos) * ones(1, nDataTmpl);
end
[val sortedIndex] = sort(matchDegInf(1, :)', 'descend');
sortedMDInf       = matchDegInf(:, sortedIndex);
% sortedFeatureL  = featureL(:, sortedIndex);

pos = 1;
while sortedMDInf(1, pos) >= 0.800
  txt = sprintf('cor = %6.3f (%3d, %3d) trans = %d angle = %2d  scale = %6.1f order = %1d\n', sortedMDInf(1, pos), sortedMDInf(3, pos), sortedMDInf(4, pos), sortedMDInf(2, pos), sortedMDInf(6, pos), sortedMDInf(7, pos), sortedMDInf(5, pos));
  pos = pos + 1;
end

% No transformation
imgFileTmpl = strcat(workDir,     tmplImgName, '.pgm');
tmplImg     = imread(imgFileTmpl);
tmplImgM    = drawMatch(tmplImg, infVecLTmpl(1, 1),  infVecLTmpl(2, 1), infVecLTmpl(5, 1), 1, pi * infVecLTmpl(4, 1) / 8);

figure(1);
imshow(tmplImgM);
imgFile  = sprintf('%s%sM.tif', workOutDir, inpImgName);
imwrite(tmplImgM, imgFile);

% Output the matched template and input images
imgFile  = strcat(workDir,     inpImgName, '.pgm');
inpImg   = imread(imgFile);
for sel = stOrd:enOrd
  % Transformed template image
  tmplN           = sortedMDInf(2, sel);
  trans           = infVecLTmpl(6,tmplN);
  imgFileTmplTran = sprintf('%stranImg%02d.pgm', workOutDir, trans);
  tmplImgTran     = imread(imgFileTmplTran);
  tmplImgTranM    = drawMatch(tmplImgTran, infVecLTmpl(1, tmplN), infVecLTmpl(2, tmplN), infVecLTmpl(5, tmplN), 1, pi * infVecLTmpl(4, tmplN) / 8);

  % input image
  inpImgM = drawMatch(inpImg, sortedMDInf(3, sel), sortedMDInf(4, sel), sortedMDInf(7, sel), 1, pi * sortedMDInf(6, sel) / 8);

  compImgM = [tmplImgTranM inpImgM];
  figure(sel - stOrd + 2);
  imshow(compImgM);
  imgFile  = sprintf('%s%s_%dM.tif', workOutDir, inpImgName, sel - stOrd + 1);
  imwrite(compImgM, imgFile);
end
