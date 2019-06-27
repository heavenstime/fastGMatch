more off

workHome    = '../images/';
tmplImgName = 'img1'; % Not dat but image for display
degradL     = {'', '/blur3', '/blur6', '/mblur3', '/mblur6', '/mblur12', '/gnoisy10', '/gnoisy20', '/gnoisy40', '/gnoisy80'};
inpImgNameL = {'img2', 'img3', 'img4', 'img5', 'img6'};
dNameL      = {'Graffiti', 'Boat', 'Bark', 'UBC', 'mnist2'};
charImgNameL = {'0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008', '0009', '0010'};

%for dgN = [1 3 6 10] % Degradation No.
for dgN = 1:1 % Degradation No.
  degrad = degradL{dgN};
  for dataType = 1:1  % Data Base No. 1:Graffiti, 2:Boat, 3:Bar', 4:UBC, 5:mnist2};
    ordHist = zeros(1, 101); 
    transData;
    if dataType == 5
      cTypeEnd = 9;
      charNEnd = 10;
    else
      cTypeEnd = 0;
      charNEnd = 1;
    end

    for cType = 0:cTypeEnd % For MNIST character type
      nMaru    = 0;
      nSankaku = 0;
      nBatsu   = 0;
      for charN = 1:charNEnd % For MNIST character type
        if dataType == 5
          dName = strcat(dNameL{dataType}, '/', num2str(cType), '/', charImgNameL{charN});
        else
          dName = dNameL{dataType};
        end

        workTmplDataDir = strcat(workHome, dName ,'/OutMatch/');
        workTmplDir     = strcat(workHome, dName,'/');
        workTmplOutDir  = strcat(workHome, dName ,'/OutMatch/');

        workInDir       = strcat(workHome, dName, degrad,'/');
        workOutDir      = strcat(workHome, dName, degrad, '/OutMatch/');

        for tlN = 1:5
          inpImgName  = inpImgNameL{tlN};
          H = Hs((3*tlN - 2):(3*tlN),:);
          stOrd       = 1;      % Start to display match point in order.
          enOrd       = 1;      % End to display match point in order.

          datFileTmpl = strcat(workTmplDataDir, tmplImgName, '.template.dat');
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
            % disp(sprintf('cor = %6.3f (%3d, %3d) trans = %d angle = %2d  scale = %6.1f order = %1d\n', sortedMDInf(1, pos), sortedMDInf(3, pos), sortedMDInf(4, pos), sortedMDInf(2, pos), sortedMDInf(6, pos), sortedMDInf(7, pos), sortedMDInf(5, pos)));
            pos = pos + 1;
	    if pos > nData
	      break
	    end
          end

          % No transformation
          imgFileTmpl = strcat(workTmplDir,     tmplImgName, '.pgm');
          tmplImg     = imread(imgFileTmpl);
          tmplImgM    = drawMatch(tmplImg, infVecLTmpl(1, 1),  infVecLTmpl(2, 1), infVecLTmpl(5, 1), 1, pi * infVecLTmpl(4, 1) / 8);

          figure(1);
          imshow(tmplImgM);
          imgFile  = sprintf('%s%sM.tif', workTmplOutDir, inpImgName);
          imwrite(tmplImgM, imgFile);

          % Output the matched template and input images
          imgFile  = strcat(workInDir,     inpImgName, '.pgm');
          inpImg   = imread(imgFile);
          for sel = stOrd:enOrd
            % Transformed template image
            tmplN           = sortedMDInf(2, sel);
            trans           = infVecLTmpl(6,tmplN);
            imgFileTmplTran = sprintf('%stranImg%02d.pgm', workTmplOutDir, trans);
            tmplImgTran     = imread(imgFileTmplTran);
            tmplImgTranM    = drawMatch(tmplImgTran, infVecLTmpl(1, tmplN), infVecLTmpl(2, tmplN), infVecLTmpl(5, tmplN), 1, pi * infVecLTmpl(4, tmplN) / 8);

            % input image
            inpImgM = drawMatch(inpImg, sortedMDInf(3, sel), sortedMDInf(4, sel), sortedMDInf(7, sel), 1, pi * sortedMDInf(6, sel) / 8);

            compImgM = [tmplImgTranM inpImgM];
            figure(sel - stOrd + 2);
            imshow(compImgM);
            imgFile  = sprintf('%s%s_%dTM.tif', workOutDir, inpImgName, sel - stOrd + 1);
            imwrite(compImgM, imgFile);
            imgFile  = sprintf('%s%s_%dM.tif', workOutDir, inpImgName, sel - stOrd + 1);
            imwrite(inpImgM, imgFile);
          end
      
%          [sortedMDInf(3, 1)  sortedMDInf(4, 1)]
          ppInvTran = inv(H) * [sortedMDInf(3, 1)  sortedMDInf(4, 1)  1]';
          pInvDiff = ppInvTran(1:2, :)/ppInvTran(3) - pCenter;
          disDiff = norm(pInvDiff);
          if (disDiff <= tmplScale * 1.4142136 * 0.2)
            nMaru = nMaru + 1;
	    res = 'Maru';
          elseif (disDiff <= tmplScale * 1.4142136)
            nSankaku = nSankaku + 1;
	    res = 'Triangle';
          else
            nBatsu = nBatsu + 1;
	    res = 'Times';
          end
	  %disp(sprintf("Degradation %d  DN = %d  IN = %d %s \n", dgN, dataType, tlN, res));
          for sel = 1:101
            if sel == 101
              ordHist(sel) = ordHist(sel) + 1;
            end
	    
	    % Ranking serach. The order first satisfied
            ppInvTran = inv(H) * [sortedMDInf(3, sel)  sortedMDInf(4, sel)  1]'; %Invers transform
            pInvDiff  = ppInvTran(1:2, :)/ppInvTran(3) - pCenter;
            disDiff   = norm(pInvDiff);
%            if (disDiff <= tmplScale * 1.4142136 * 0.2)
            if (disDiff <= tmplScale * 1.4142136 )
              ordHist(sel) = ordHist(sel) + 1;
              break
            end
          end
          clear matchDegInf val sortedIndex sortedMDInf
        end
      end
      disp(sprintf('%d  %d  %d  %d\n', cType, nMaru, nSankaku, nBatsu))
    end
  end
end
