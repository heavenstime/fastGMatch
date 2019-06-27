/*
 * fastGMatch.h
 *
 *  Created on: 2019/03/13
 *      Author: yamasita
 */

#ifndef FASTGMATCH_H_
#define FASTGMATCH_H_



class FastGMatch {
public:
	int nx, ny, nxy;
	int nAngleCoef;
	int nFeature = 0;
	int K2 = 0, K4 = 0;
	int step = 0, order = 0, pos = 0, stPos = 0, stX = 0, stY = 0, endX = 0, endY = 0;
	int ix = 0, iy = 0, iFeature = 0, transType = 0;
	int nMaxTheta = 0;
	int maxITheta[2];
	int mkTemplate = 0, nTrans = 0;
	double scale = 0.0, scaleRatio = 0.0, scaleMax = 0.0;

	std::string fileName;
	std::string dFileName;      /* Data file name*/
	std::string tmplDFileName;  /* template data file name */
	int printTime = 1;

	clock_t start = 0.0, now = 0.0;

	int       *inImgI;
	FPTYPE    *inImgOrg;
	FPTYPE    *inImg;
	Feature   **features;
	Feature   **featuresTmpl;
	Feature   **featuresP = NULL;
	ExtFeature *extFeature;
	ImageIO   *imageIO;
	GptTbl    *gptTbl;
	Utilities *utilities;

	FastGMatch(int nx, int ny, int mkTemplate);
	int fastGMatch();
	int histSave(FPTYPE *dirHist);
};

#endif /* FASTGMATCH_H_ */
