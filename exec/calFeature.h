/*
 * fastFeature.h
 *
 *  Created on: 2018/03/29
 *      Author: yamasita
 */

#ifndef CALFEATURE_H_
#define CALFEATURE_H_

typedef struct {
	int ix, iy;    /* Position */
	int ordHist;   /* 0 or 1 */
	int iTheta;
	FPTYPE scale;  /* Scale */
	int transType; /* Transformation type */
	FPTYPE *vector;
} Feature;

typedef struct {
	int nx, ny, nxy;
	int nFeature;
	int nAngleCoef;
	int nEstFeature;
	FPTYPE  *relativePosL;
	int     *relativePosScL;
	Feature *featureL;
	FPTYPE  *tmpVector;
} Features;

typedef struct {
	int K;
	int nInt;
	int nOrd;
	int maxNInt;
	int extType; /* 0 : zero extension, 1 : extension */
	FPTYPE *cosL;
	FPTYPE *sinL;
	FPTYPE *lineExt;
	FPTYPE *intSin;
	FPTYPE *intCos;
	FPTYPE *blurCoef;
	FPTYPE *diffCoef;
} WorkIIR;

typedef struct {
	int nx, ny;
	FPTYPE *xBlurImg;
	FPTYPE *xDiffXImg;
	FPTYPE *trans0;
	FPTYPE *trans1;
	FPTYPE *trans2;
} WorkImg;

typedef struct {
	int nx, ny;
	FPTYPE *invFourTbl;
	FPTYPE *rotTbl;
	FPTYPE *largeScaleRelTbl;
	int    *largeScaleTbl;
	FPTYPE *largeScaleATbl;
} WorkMT; /* Work max for theta*/

typedef struct {
	int nTran;
	FPTYPE *gptTbl;
} GptTbl; /* Work max for theta*/

int calDirHistPoint(FPTYPE *diffXImg, FPTYPE *diffYImg, FPTYPE *dirHist, int nx, int ny);
int maxDirection(FPTYPE *blurLargeDirHist, WorkMT *workMT, int *maxITheta);
int rotateDirHist(FPTYPE *vectorIn, FPTYPE *vectorOut, int stRot, WorkMT *workMT);
int mkFeature(FPTYPE *blurSmallDirHist, WorkMT *workMT, Features *features, Feature *feature);
int approxLargeScale(FPTYPE *blurSmallDirHist, FPTYPE *outHist, int ix, int iy, WorkMT *workMt);

int gaussDiff(FPTYPE *inImg, FPTYPE *diffXImg,  FPTYPE *diffYImg, WorkImg *workImg, WorkIIR *workIIR) ;
int gaussSmooth(FPTYPE *inImg, FPTYPE *blurImg, WorkImg *workImg, WorkIIR *workIIR);
int gaussSmooth1(int type, FPTYPE *inImg, FPTYPE *bImg, FPTYPE *dImg, int nx, int ny,  WorkIIR *workIIR);
int sincosIIRF(int nData, WorkIIR *workFP);
int transposition(FPTYPE *inImg, FPTYPE *outImg, int nx, int ny);

#endif /* CALFEATURE_H_ */

