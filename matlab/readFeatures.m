function [featureL infVecL nx ny] = readFeatures(fileName)
  fp = fopen(fileName, "r");
  if fp == 0
    printf("     Cannot open! %s \n\n", fileName);
    exit(1);
  end

  nx         = fread(fp, 1, "int32");
  ny         = fread(fp, 1, "int32");
  nFeature   = fread(fp, 1, "int32");
  nAngleCoef = fread(fp, 1, "int32");

  infVecL    = zeros(5, nFeature);
  featureL   = zeros(nAngleCoef * 4 * 4, nFeature);

  for pos = 1:nFeature
    infVecL(1, pos)  = fread(fp, 1, "int32");   % ix
    infVecL(2, pos)  = fread(fp, 1, "int32");   % iy
    infVecL(3, pos)  = fread(fp, 1, "int32");   % order in histrgram 0/1
    infVecL(4, pos)  = fread(fp, 1, "int32");   % iTheta
    infVecL(5, pos)  = fread(fp, 1, "float32"); % scale
    infVecL(6, pos)  = fread(fp, 1, "int32");   % transType
    featureL(:, pos) = fread(fp, nAngleCoef * 4 * 4, "float32"); % scale
  end
fclose(fp);

end
