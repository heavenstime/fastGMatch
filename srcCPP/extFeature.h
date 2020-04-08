/*
 * extFeature.h
 *
 *  Created on: 2019/05/09
 *      Author: yamasita
 */

#ifndef EXTFEATURE_H_t

#define EXTFEATURE_H_

class WorkImg {
public:
	WorkImg(int nx, int ny);
	int nxy;
	FPTYPE *xBlurImg;
	FPTYPE *xDiffXImg;
};

class WorkMaxT {
public:
	WorkMaxT(int nx, int ny);
	int maxDirection(FPTYPE *blurLargeDirHist, int *maxITheta);
	int nx, ny;
	int nAngleCoef = 2 * PANGLE + 1;;
	FPTYPE *invFourTbl;
	FPTYPE *rotTbl;
	FPTYPE *largeScaleRelTbl;
	int    *largeScaleTbl;
	FPTYPE *largeScaleATbl;
	FPTYPE *tmpVector;
}; /* Work max for theta*/

class WorkIIR {
public:
	WorkIIR(int nOrd0, int maxNxy, FPTYPE *coefG, FPTYPE *coefDG);
	int K = 0;
	int nInt = 0;
	int P;
	int extType = 0; /* 0 : zero extension, 1 : extension */
	FPTYPE *cosL;
	FPTYPE *sinL;
	FPTYPE *lineExt;
	FPTYPE *blurCoef;
	FPTYPE *diffCoef;
	FPTYPE *inteCos;
	FPTYPE *inteSin;
	FPTYPE *cosInitCoef;
	FPTYPE *sinInitCoef;
	int initByK(int K);
};

class ImageProc {
public:
	int nx, ny, nxy;
	WorkImg  *workImg;
	WorkMaxT  *workMaxT;
	ImageProc(int nx0, int ny0, WorkImg *workImg0, WorkMaxT *workMT0);
	int calDirHistPoint(FPTYPE *diffXImg, FPTYPE *diffYImg, FPTYPE *dirHist);
	int rotateDirHist(FPTYPE *vectorIn, FPTYPE *vectorOut, int stRot);
	int approxLargeScale(FPTYPE *blurSmallDirHist, FPTYPE *outHist, int ix, int iy);

	int  gaussDiff(FPTYPE *inImg, FPTYPE *diffXImg,  FPTYPE *diffYImg, WorkIIR *workIIR) ;
	int  gaussSmooth(FPTYPE *inImg, FPTYPE *blurImg, WorkIIR *workIIR);
	int  gaussSmoothSubTr(int type, FPTYPE *inImg, FPTYPE *bImg, FPTYPE *dImg, WorkIIR *workIIR);
};

class Feature {
public:
	int ix = 0, iy = 0;    /* Position */
	int ordHist = 0;   /* 0 or 1 */
	int iTheta = 0;
	FPTYPE scale = 0.0;  /* Scale */
	int transType = 0; /* Transformation type */
	FPTYPE *vector;
	Feature(int ix, int iy, int ordHist, int *maxIThetaL, double scale, FPTYPE *histFourS4, int *relativePosScL, WorkMaxT *workMT, ImageProc *imageProc);
	Feature();
	void prFeatureInf();
	void prFeatureVal();
};

class ExtFeature {
public:
	int nx, ny, nxy;
	int nAngleCoef = 2 * PANGLE + 1;
	WorkIIR *workIIR1;
	WorkIIR *workIIR2;
	WorkIIR *workIIR4;
	WorkImg *workImg;
	ImageProc *imageProc;
	WorkMaxT *workMaxT;
	ImageIO *imageIO = NULL; /* For debug */

	FPTYPE *diffXImg, *diffYImg;
	FPTYPE *dirHist, *histFourSL, *histFourSS;
	FPTYPE *relativePosL;
	int    *relativePosScL;

	ExtFeature(int nx, int ny);
	int extract(int cont, FPTYPE *inImg, int ix, int iy, double scale, Feature **feature1, Feature **feature2);
	void calHistgram();
	int  histSave(FPTYPE *dirHist, ImageIO *imageIO);
	FPTYPE *getImageFp(int);
	void loadFeatures(const char *fileName, int *nFeatureP, int *nAngleCoefP, Feature **featureL);
	void saveFeatures(const char *fileName, int nFeature, int nAngleCoef, Feature **featureL);
};

#endif /* EXTFEATURE_H_ */
