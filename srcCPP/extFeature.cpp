/*
 * extFeature.cpp
 *
 *  Created on: 2019/05/09
 *      Author: yamasita
 */

#include <cmath>
#include <string>
#include <iostream>

#include "parameter.h"
#include "utilities.h"
#include "extFeature.h"
#include "fastGMatch.h"

ExtFeature::	ExtFeature(int nx, int ny) {
	this->nx = nx;
	this->ny = ny;
	nxy = nx * ny;

	/* Initialize work image regions */
	workImg   = new WorkImg(nx, ny);

	/* Initialize workMT */
	workMaxT    = new WorkMaxT(nx, ny);
	imageProc   = new ImageProc(nx, ny, workImg, workMaxT);

	/* Coefficients P = 2, sigma = 1.1  created by totalError201803.m */
	FPTYPE coefG2[]   = {1.5846315202e-01, 1.7508079293e-01, 2.7323762984e-02};
	FPTYPE coefDG2[]  = {0.0000000000e+00, -1.7506270969e-01, -5.4683692457e-02};

	/* Coefficients P = 4, sigma = 0.8 */
	FPTYPE coefG4[]   = {1.5914081930e-01,  2.3116780826e-01,  8.8476942331e-02,  1.7890179087e-02,  1.8836364529e-03};
	FPTYPE coefDG4[]  = {0.0000000000e+00, -2.3116693143e-01, -1.7695563833e-01, -5.3667906760e-02, -7.5380531476e-03};

	/* Work region for IIR filter */
	int maxNxy = (nx > ny) ? nx : ny;
	workIIR2  = new WorkIIR(2, maxNxy, coefG2, coefDG2);
	workIIR4  = new WorkIIR(4, maxNxy, coefG4, coefDG4);

	/* Differential Image */
	diffXImg       = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	diffYImg       = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	/* Histogram  region */
	dirHist        = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy * nAngleCoef);
	histFourSL     = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy * nAngleCoef);
	histFourSS     = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy * nAngleCoef);

	/* Work region for mesh of features */
	relativePosL   = (FPTYPE *) malloc(sizeof(FPTYPE) * NMAXTHETA * GRIDSIZE * GRIDSIZE * 2);
	relativePosScL = (int *) malloc(sizeof(int) * NMAXTHETA * GRIDSIZE * GRIDSIZE * 2);
	int pos = 0;
	for (int iTheta = 0 ; iTheta < NMAXTHETA ; ++iTheta) {
		FPTYPE cosTheta = cos(2 * PI * iTheta / NMAXTHETA);
		FPTYPE sinTheta = sin(2 * PI * iTheta / NMAXTHETA);
		FPTYPE y = - (GRIDSIZE - 1.0) / GRIDSIZE;
		for (int iy = 0 ; iy < GRIDSIZE ; ++iy) {
			FPTYPE x = - (GRIDSIZE - 1.0) / GRIDSIZE;
			for (int ix = 0 ; ix < GRIDSIZE ; ++ix) {
				relativePosL[pos] = x * cosTheta - y * sinTheta; ++pos;// x
				relativePosL[pos] = x * sinTheta + y * cosTheta; ++pos;// y
				x = x + 2.0 / GRIDSIZE;
			}
			y = y + 2.0 / GRIDSIZE;
		}
	}
}

/* cont == 0: first time for the scale, cont == 1: the previous calculated data can be used. */
int ExtFeature::extract(int cont, FPTYPE *inImg, int ix, int iy, double scale, Feature **feature1, Feature **feature2) {
	if (cont == 0) {

		/* Calculate differential images */
		int 	K4 = (int) (PI * scale * GRADRATIO / (GRIDSIZE * SIGMA4) + ROUNDFRAC);
		workIIR4->extType = 1;
		workIIR4->initByK(K4);

		imageProc->gaussDiff(inImg, diffXImg, diffYImg, workIIR4); //Calculate differential  image;
		//imageProc->gaussSmooth(inImg, diffXImg, workIIR4); //Calculate differential  image;
		/* Output differential image */
		std::string fileNameDX(WORKBASE); fileNameDX = fileNameDX + IMGPRE + BLURPRE + OUTPRE + "diffX.pgm";
		std::string fileNameDY(WORKBASE); fileNameDY = fileNameDY + IMGPRE + BLURPRE + OUTPRE + "diffY.pgm";
		imageIO->saveImageFileFp(fileNameDX.c_str(), diffXImg);
		imageIO->saveImageFileFp(fileNameDY.c_str(), diffYImg);
//exit(0);
		/* Fourier series expression of directions */
		imageProc->calDirHistPoint(diffXImg, diffYImg, dirHist);

#ifdef USELARGESCALE
		/* Blurred directions for large scale */
		int K2L = (int) (PI * scale / SIGMA2 + ROUNDFRAC);
		workIIR2->extType = 0;
		workIIR2->initByK(K2L);
		int stPos = 0;
		for (int order = 0 ; order < 2 * PANGLE + 1 ; ++order) {
			imageProc->gaussSmooth(&(dirHist[stPos]), &(histFourSL[stPos]), workIIR2);
			//printf("order %d  %f   %f\n", order, dirHist[stPos + nx / 2 + nx * ny /2], blurLargeDirHist[stPos + nx / 2 + nx * ny /2] );
			stPos += nxy;
		}
#else
		for (pos = 0 ; pos < 2 * NAPPROPOINT ; ++pos) {
			workMT->largeScaleTbl[pos] = (int) (scale * workMT->largeScaleRelTbl[pos] + ROUNDFRAC);
		}
#endif

		/* Smooth directions for small scale */
		int K2S = (int) (PI * scale / (SIGMA2 * GRIDSIZE) + ROUNDFRAC);
		workIIR2->extType = 0;
		workIIR2->initByK(K2S);
		stPos = 0;
		for (int order = 0 ; order < 2 * PANGLE + 1; ++order) {
			imageProc->gaussSmooth(&(dirHist[stPos]), &(histFourSS[stPos]), workIIR2);
			stPos += nxy;
		}

		/* Calculate grid points for feature vector. */
		for (int pos = 0 ; pos < NMAXTHETA * GRIDSIZE * GRIDSIZE * 2 ; ++pos) {
			relativePosScL[pos] = (int) (scale * relativePosL[pos] + ROUNDFRAC);
			//printf("relativePosScL[%d] = %d    %lf \n", pos, relativePosScL[pos], relativePosL[pos]);
		}
	}
#ifndef USELARGESCALE
	imageProc->approxLargeScale(histFourS4, & (histFourS1[ix + nx * iy]), ix, iy, workMT);
#endif
	int maxIThetaL[2];
	int nMaxTheta = workMaxT->maxDirection(& (histFourSL[ix + nx * iy]), maxIThetaL);

	if (nMaxTheta >= 1) *feature1 = new Feature(ix, iy, 0, maxIThetaL, scale, histFourSS, relativePosScL, workMaxT, imageProc);
	if (nMaxTheta >= 2) *feature2 = new Feature(ix, iy, 1, maxIThetaL, scale, histFourSS, relativePosScL, workMaxT, imageProc);
	return nMaxTheta;
}

FPTYPE *ExtFeature::getImageFp(int type) {
	switch(type) {
	case 0:
		return diffXImg;
	case 1:
		return diffYImg;
	case 2:
		return dirHist;
	case 3:
		return histFourSL;
	case 4:
		return histFourSS;
	}
	return NULL;
}

Feature::Feature(int ix, int iy, int ordHist, int *maxIThetaL, double scale, FPTYPE *histFourSS, int *relativePosScL, WorkMaxT *workMaxT, ImageProc *imageProc) {

	double sqSum, sqr;
	int *rePosSc;
	int iTheta = maxIThetaL[ordHist];
	this->ix = ix; this->iy = iy; this->iTheta = iTheta; this->ordHist = ordHist; this->scale = scale;
	vector =  (FPTYPE *) malloc(sizeof(FPTYPE) * (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE);

	rePosSc = &(relativePosScL[iTheta * GRIDSIZE * GRIDSIZE * 2]);

	int posXY = 0, posXY2 = 0;
	for (int ify = 0 ; ify < GRIDSIZE ; ++ify) {
		for (int ifx = 0 ; ifx < GRIDSIZE ; ++ifx) {
			//printf("rePosSc (%d, %d)  \n", rePosSc[posXY2], rePosSc[posXY2 + 1]);
			int ixg = ix + rePosSc[posXY2++];
			int iyg = iy + rePosSc[posXY2++];
			// printf("(%d, %d)\n", ix, iy);
			int posHist = ixg + imageProc->nx * iyg;
			//printf("%d %d \n", histFourSS, workMaxT->tmpVector);
			workMaxT->tmpVector[0] = histFourSS[posHist];
			for (int order = 1 ; order < workMaxT->nAngleCoef ; ++order) {
				posHist += imageProc->nxy;
				workMaxT->tmpVector[order] = histFourSS[posHist];
			}
			imageProc->rotateDirHist(workMaxT->tmpVector, &(vector[workMaxT->nAngleCoef * (posXY++)]), iTheta * workMaxT->nAngleCoef);
			// for (ix = 0 ; ix < nAngleCoef ; ++ix) printf("iCoef %d  %f, ", ix, feature->vector[nAngleCoef * (posXY - 1) + ix]);
			// printf("nAng  %d \n", nAngleCoef);
		}
	}
	/* Normalization */
	sqSum = 0.0;
	for (int pos = 0 ; pos < workMaxT->nAngleCoef * GRIDSIZE * GRIDSIZE ; ++pos) sqSum += vector[pos] * vector[pos];
	sqr = sqrt(sqSum);
	for (int pos = 0 ; pos < workMaxT->nAngleCoef * GRIDSIZE * GRIDSIZE ; ++pos) vector[pos] /= sqr;
}

/* Print information of features */
void Feature::prFeatureInf() {
	printf("Feature Info. (%3d, %3d) order = %2d, angle = %2d, scale = %7.2f\n", ix, iy, ordHist, iTheta, scale);
}

/* For load feature */
Feature::Feature() {
	vector =  (FPTYPE *) malloc(sizeof(FPTYPE) * (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE);
}

/* Load features */
void ExtFeature::loadFeatures(const char *fileName, int *nFeatureP, int *nAngleCoefP, Feature **featureL) {
	int  pos, ret = 0;
	Feature *feature;
	FILE *fp; /* File pointer */
	fp = fopen(fileName, "rb");
	if (fp == NULL) {
		printf("     Cannot open! %s \n\n", fileName);
		exit(1);
	}

	ret = fread(&(nx), sizeof(int), 1, fp);
	if (ret == 0) 	printf("Error occurs ! nx\n");
	ret = fread(&(ny), sizeof(int), 1, fp);
	if (ret == 0) 	printf("Error occurs ! ny \n");
	ret = fread(nFeatureP, sizeof(int), 1, fp);
	if (ret == 0) 	printf("Error occurs ! nFeature \n");
	ret = fread(nAngleCoefP, sizeof(int), 1, fp);
	if (ret == 0) 	printf("Error occurs ! nAngleCoef \n");
	for (pos = 0 ; pos < *nFeatureP ; ++pos) {
		featureL[pos] = feature = new Feature();
		ret = fread(&(feature->ix), sizeof(int), 1, fp);
		if (ret == 0) 	printf("Error occurs ! ix\n");
		ret = fread(&(feature->iy), sizeof(int), 1, fp);
		if (ret == 0) 	printf("Error occurs ! iy\n");
		ret = fread(&(feature->ordHist), sizeof(int), 1, fp);
		if (ret == 0) 	printf("Error occurs ! ordHist\n");
		ret = fread(&(feature->iTheta), sizeof(int), 1, fp);
		if (ret == 0) 	printf("Error occurs ! iTheta\n");
		ret = fread(&(feature->scale), sizeof(FPTYPE), 1, fp);
		if (ret == 0) 	printf("Error occurs ! scale\n");
		ret = fread(&(feature->transType), sizeof(int), 1, fp);
		if (ret == 0) 	printf("Error occurs ! transType\n");
		ret = fread(feature->vector, sizeof(FPTYPE), (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE, fp);
		if (ret == 0) 	printf("Error occurs ! vectors \n");
	}
}

/* Save features */
void ExtFeature::saveFeatures(const char *fileName, int nFeature, int nAngleCoef, Feature **featureL) {
	int  pos;
	Feature *feature;
	FILE *fp; /* File pointer */
	fp = fopen(fileName, "wb");
	if (fp == NULL) {
		printf("     Cannot open! %s \n\n", fileName);
		exit(1);
	}

	fwrite(&(nx), sizeof(int), 1, fp);
	fwrite(&(ny), sizeof(int), 1, fp);
	fwrite(&(nFeature), sizeof(int), 1, fp);
	fwrite(&(nAngleCoef), sizeof(int), 1, fp);

	for (pos = 0 ; pos < nFeature ; ++pos) {
		// printf("nFeature %d, pos %d \n", nFeature, pos);
		feature = featureL[pos];
		fwrite(&(feature->ix), sizeof(int), 1, fp);
		fwrite(&(feature->iy), sizeof(int), 1, fp);
		fwrite(&(feature->ordHist), sizeof(int), 1, fp);
		fwrite(&(feature->iTheta), sizeof(int), 1, fp);
		fwrite(&(feature->scale), sizeof(FPTYPE), 1, fp);
		fwrite(&(feature->transType), sizeof(int), 1, fp);
		fwrite(feature->vector, sizeof(FPTYPE), (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE, fp);
	}
}

ImageProc::ImageProc(int nx, int ny, WorkImg *workImg, WorkMaxT *workMaxT) {
		this->nx      = nx;
		this->ny      = ny;
		nxy           = nx * ny;
		this->workImg  = workImg;
		this->workMaxT = workMaxT;
}
/* Calculate Fourier expression of directional histogram */
int ImageProc::calDirHistPoint(FPTYPE *diffXImg, FPTYPE *diffYImg, FPTYPE *dirHist) {
	int ixy, nxy = nx * ny, nxy2 = 2 * nxy, stCos, stSin, order;
	FPTYPE diffX, diffY, diffPower, diffSqrt, ansCos, ansSin, diffCos, diffSin;

	for (ixy = 0 ; ixy < nxy ; ++ixy) {
		diffX     = diffXImg[ixy];
		diffY     = diffYImg[ixy];
		diffPower = diffX * diffX + diffY * diffY;

		/* Histogram with regularization */
		dirHist[ixy] = diffSqrt = sqrt(diffPower + EDGEPOWER);
		stCos = nxy; stSin = nxy2; /* 0 is for DC component */
		dirHist[ixy + stCos] = ansCos =  diffX;
		dirHist[ixy + stSin] = ansSin = - diffY;
		diffCos = diffX / diffSqrt;
		diffSin = diffY / diffSqrt;
		for (order = 2 ; order <= PANGLE ; ++order) {
			stCos = stCos + nxy2, stSin = stSin + nxy2;
			dirHist[ixy + stCos]          =  ansCos * diffCos + ansSin * diffSin;
			dirHist[ixy + stSin] = ansSin = -ansCos * diffSin + ansSin * diffCos;
			ansCos = dirHist[ixy + stCos];
		}
	}
	return 0;
}

/* Rotation of directional histogram */
int ImageProc::rotateDirHist(FPTYPE *vectorIn, FPTYPE *vectorOut, int stRotTbl) {
	int order, posRotTbl, posVect;
	FPTYPE vCos, vSin, rCos, rSin;

	posRotTbl  = stRotTbl + 1;
	posVect    = 0;
	vectorOut[posVect] = vectorIn[posVect]; ++posVect;
	for (order = 1 ; order <= PANGLE ; ++order) {
		vCos = vectorIn[posVect];
		vSin = vectorIn[posVect + 1];
		rCos = workMaxT->rotTbl[posRotTbl++];
		rSin = workMaxT->rotTbl[posRotTbl++];
		vectorOut[posVect++] =  vCos * rCos - vSin * rSin;
		vectorOut[posVect++] =  vCos * rSin + vSin * rCos;
	}
	return 0;
}

/* Large scale smoothing is approximated by small scale smoothing */
int ImageProc::approxLargeScale(FPTYPE *histFourS4, FPTYPE *histFourAS1, int ix, int iy) {
	int pInd, pos, tPos, tPos2, order, orderP;
	int nxy = nx * ny;
	FPTYPE a;

	/* Initialize the output */
	for (order = 0 ; order < 2 * PANGLE + 1 ; ++order) 	histFourAS1[order] = 0.0;

	tPos = tPos2 = 0;
	for (pInd = 0 ; pInd < NAPPROPOINT ; ++pInd) {
		//printf("(%d, %d) (%d, %d) a = %f \n", ix, iy, ix + workMT->largeScaleTbl[tPos2], iy + workMT->largeScaleTbl[tPos2+1], workMT->largeScaleATbl[tPos]);
		pos  = (ix + workMaxT->largeScaleTbl[tPos2]);
		pos += nx * (iy + workMaxT->largeScaleTbl[tPos2++]);
		a    = workMaxT->largeScaleATbl[tPos++];
		orderP = 0;
		for (order = 0 ; order < 2 * PANGLE + 1 ; ++order) {
			//printf("orderP = %d, pos = %d \n", orderP, pos);
			histFourAS1[orderP] += a * histFourS4[pos];
			pos += nxy; orderP += nxy;
		}
	}
	//for (order = 0 ; order < 2 * PANGLE + 1 ; ++order) printf("order %d %lf\n", order,	outHist[order]);
	return 0;
}

/* Differential with Gauss smoothing */
int ImageProc::gaussDiff(FPTYPE *inImg, FPTYPE *diffXImg,  FPTYPE *diffYImg, WorkIIR *workIIR) {
	ImageProc *imageProcTrans = new ImageProc(ny, nx, workImg, workMaxT);

	/* X direction */;
	gaussSmooth1Tr(2, inImg, workImg->xBlurImg, workImg->xDiffXImg,  workIIR);

	/* Y direction */
	imageProcTrans->gaussSmooth1Tr(1, workImg->xBlurImg, NULL, diffYImg,  workIIR);
	imageProcTrans->gaussSmooth1Tr(0, workImg->xDiffXImg, diffXImg, NULL,  workIIR);

	return 0;
}

/* Gauss smoothing */
int ImageProc::gaussSmooth(FPTYPE *inImg, FPTYPE *blurImg, WorkIIR *workIIR) {
	ImageProc *imageProcTrans = new ImageProc(ny, nx, workImg, workMaxT);

	/* X direction */;
	gaussSmooth1Tr(0, inImg, workImg->xBlurImg, NULL, workIIR);
	/* Y direction */
	imageProcTrans->gaussSmooth1Tr(0, workImg->xBlurImg, blurImg, NULL,  workIIR);
	return 0;
}

/* One dimensional transform (by interval integral) with Transposition */
int ImageProc::gaussSmooth1Tr(int type, FPTYPE *inImg, FPTYPE *bImg, FPTYPE *dImg, WorkIIR *workIIR) {
	int K = workIIR->K, K2 = 2 * K;
	int nInte = nx + K;
	int P = workIIR->P;
  FPTYPE extS, inteCosTmp, inteSinTmp, add, sub, val1, val2;
  FPTYPE *lineAdd = workIIR->lineAdd;
  FPTYPE *lineSub = workIIR->lineSub;
	FPTYPE *cosL    = workIIR->cosL;
	FPTYPE *sinL    = workIIR->sinL;
	FPTYPE *inteCos = workIIR->inteCos;
	FPTYPE *inteSin = workIIR->inteSin;
	FPTYPE *cosInitCoef = workIIR->cosInitCoef;
	FPTYPE *sinInitCoef = workIIR->sinInitCoef;

	FPTYPE *blurCoef = workIIR->blurCoef;
	FPTYPE *diffCoef = workIIR->diffCoef;


	/* Fixed Extension region */
	int H = (nx > K)? K2 : (nx + K);
	//printf("%d %d %d  %d \n", workIIR->lineSub, lineSub,  workIIR->lineAdd, lineAdd);
	for (int pos = 0  ; pos < K ; ++pos) lineSub[pos] = 0.0;
	if (workIIR->extType == 0) {
		for (int pos = nx ; pos < nx + K ; ++pos) lineAdd[pos] = 0.0;
		for (int pos = K  ; pos < H      ; ++pos) lineSub[pos] = 0.0;
	}

	for (int iy = 0 ; iy < ny ; ++iy) {
		int outPos = iy;
		FPTYPE *vect = &(inImg[nx * iy]);
    for (int pos = 0     ; pos < nx    ; ++pos) lineAdd[pos] = vect[pos];
    for (int pos = K2    ; pos < nInte ; ++pos) lineSub[pos] = vect[pos - K2];
    if (workIIR->extType == 1) {
	    for (int pos = nx    ; pos < nInte ; ++pos) lineAdd[pos] = vect[nx - 1];
	    for (int pos = K     ; pos < H     ; ++pos) lineSub[pos] = vect[0];
	  }
	   /* Initial value of integral signal */
	   extS = (workIIR->extType == 0) ? 0.0 : vect[0];
	   inteCos[0] = K * extS;
	   for (int p = 1 ; p <= P ; ++p) {
	     inteCos[p] = cosInitCoef[p] * extS;
	     inteSin[p] = sinInitCoef[p] * extS;
	   }
	   switch(type) {
	   case 0: /* Only blurring */
	  	 for (int pos = 0 ; pos < K ; ++pos) {
	  		 inteCos[0] += add = lineAdd[pos];
	  		 for (int p = 1 ; p <= P ; ++p) {
	  			 inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p];
	  			 inteCos[p] = cosL[p] * inteCosTmp - sinL[p] * inteSinTmp + add;
	  			 inteSin[p] = sinL[p] * inteCosTmp + cosL[p] * inteSinTmp;
	  		 }
	  	 }
	  	 for (int pos = K ; pos < nInte ; ++pos) {
	  		 inteCos[0] += add = lineAdd[pos];
	  		 val1 =  blurCoef[0] * inteCos[0];
	  		 inteCos[0] -= sub = lineSub[pos];
	  		 for (int p = 1 ; p <= P ; ++p) {
	  			 inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p];
	  			 inteCos[p] = cosL[p] * inteCosTmp - sinL[p] * inteSinTmp + add;
	  			 inteSin[p] = sinL[p] * inteCosTmp + cosL[p] * inteSinTmp;
	  			 val1 += blurCoef[p] * inteCos[p];
	  			 inteCos[p] -= sub;
	  		 }
	  		 //printf("val %lf \n", val1); printf("outPos = %d \n", outPos);
	  		 bImg[outPos] = val1;
	  		 outPos += ny;
	  	 }
	  	 break;
	   case 1: /* Only differential */
	  	 for (int pos = 0 ; pos < K ; ++pos) {
	  		 add = lineAdd[pos];
		  	 sub = lineSub[pos];
	  		 for (int p = 1 ; p <= P ; ++p) {
	  			 inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p];
	  			 inteCos[p] = cosL[p] * inteCosTmp - sinL[p] * inteSinTmp + add;
	  			 inteSin[p] = sinL[p] * inteCosTmp + cosL[p] * inteSinTmp;
	  		 }
	  	 }
	  	 for (int pos = K ; pos < nInte ; ++pos) {
	  		 add = lineAdd[pos];
		  	 val2 = 0.0;
		  	 sub = lineSub[pos];
	  		 for (int p = 1 ; p <= P ; ++p) {
	  			 inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p];
	  			 inteCos[p] = cosL[p] * inteCosTmp - sinL[p] * inteSinTmp + add;
	  			 inteSin[p] = sinL[p] * inteCosTmp + cosL[p] * inteSinTmp;
	  			 val2 += diffCoef[p] * inteSin[p];
	  			 inteCos[p] -= sub;
	  		 }
	  		 dImg[outPos] = val2;
	  		 outPos += ny;
	  	 }
	  	 break;
	   case 2: /* Blur and differential */
	  	 for (int pos = 0 ; pos < K ; ++pos) {
	  		 inteCos[0] += add = lineAdd[pos];
	  		 for (int p = 1 ; p <= P ; ++p) {
	  			 inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p];
	  			 inteCos[p] = cosL[p] * inteCosTmp - sinL[p] * inteSinTmp + add;
	  			 inteSin[p] = sinL[p] * inteCosTmp + cosL[p] * inteSinTmp;
	  		 }
	  	 }
	  	 for (int pos = K ; pos < nInte ; ++pos) {
	  		 inteCos[0] += add = lineAdd[pos];
	  		 val1 = blurCoef[0] * inteCos[0];
	  		 val2 = 0.0;
	  		 inteCos[0] -= sub = lineSub[pos];
	  		 for (int p = 1 ; p <= P ; ++p) {
	  			 inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p];
	  			 inteCos[p] = cosL[p] * inteCosTmp - sinL[p] * inteSinTmp + add;
	  			 inteSin[p] = sinL[p] * inteCosTmp + cosL[p] * inteSinTmp;
	  			 val1 += blurCoef[p] * inteCos[p];
	  			 val2 += diffCoef[p] * inteSin[p];
	  			 inteCos[p] -= sub;
	  		 }
	  		 bImg[outPos] = val1;
	  		 dImg[outPos] = val2;
	  		 outPos += ny;
	  	 }
	  	 break;
	   }
	}
	return 0;
}

WorkImg::WorkImg(int nx, int ny) {
	this->nxy = nx * ny;
	xBlurImg  = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	xDiffXImg = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
}

/* Constructor */
WorkMaxT::WorkMaxT(int nx, int ny) {
	this->nx         = nx;
	this->ny         = ny;
	invFourTbl       = (FPTYPE *) malloc(sizeof(FPTYPE) * 2 * NMAXTHETA * (2 * PANGLE + 1));
	rotTbl           = (FPTYPE *) malloc(sizeof(FPTYPE) * 2 * NMAXTHETA * (2 * PANGLE));
	largeScaleRelTbl = (FPTYPE *) malloc(sizeof(FPTYPE) * 2 * NAPPROPOINT);
	largeScaleTbl    = (int *)    malloc(sizeof(int)    * 2 * NAPPROPOINT);
	largeScaleATbl   = (FPTYPE *) malloc(sizeof(FPTYPE) * NAPPROPOINT);
	tmpVector         = (FPTYPE *) malloc(sizeof(FPTYPE) * NMAXTHETA * (2 * PANGLE + 1));

	/* For inverse Fourier transform of max angle*/
	double dTheta    = 2.0 * PI / NMAXTHETA;
	int pos = 0, order;
	for (int posTheta = 0 ; posTheta < NMAXTHETA ; ++posTheta) {
		invFourTbl[pos++] = 0.5 / PI;
		for (order = 1 ; order <= PANGLE ; ++order) {
			invFourTbl[pos++] =   cos(dTheta * order * posTheta) / PI;
			invFourTbl[pos++] = - sin(dTheta * order * posTheta) / PI;
		}
	}

	/* For rotation of direction histogram */
	pos = 0;
	for (int posTheta = 0 ; posTheta <= NMAXTHETA ; ++posTheta) {
		rotTbl[pos++] = 1.0;
		for (order = 1 ; order <= PANGLE ; ++order) {
			rotTbl[pos++] = cos(dTheta * order * posTheta);
			rotTbl[pos++] = sin(dTheta * order * posTheta);
		}
	}
	/* For approximation of Gaussian function of large variance by small variance */
#ifndef USELARGESCALE
	pos = order = 0;
	largeScaleRelTbl[pos++] = 0;
	largeScaleRelTbl[pos++] = 0;
	largeScaleATbl[order++] = A0AP;
	dTheta = 2 * PI / N1AP;
	for (int posTheta = 0 ; posTheta < N1AP ; ++posTheta) {
		largeScaleRelTbl[pos++] = R1AP * cos(dTheta * posTheta);
		largeScaleRelTbl[pos++] = R1AP * sin(dTheta * posTheta);
		largeScaleATbl[order++] = A1AP;
	}
	dTheta = 2 * PI / N2AP;
	for (int posTheta = 0 ; posTheta < N2AP ; ++posTheta) {
		largeScaleRelTbl[pos++] = R2AP * cos(dTheta * posTheta);
		largeScaleRelTbl[pos++] = R2AP * sin(dTheta * posTheta);
		largeScaleATbl[order++] = A2AP;
	}
	dTheta = 2 * PI / N3AP;
	for (int posTheta = 0 ; posTheta < N3AP ; ++posTheta) {
		largeScaleRelTbl[pos++] = R3AP * cos(dTheta * posTheta);
		largeScaleRelTbl[pos++] = R3AP * sin(dTheta * posTheta);
		largeScaleATbl[order++] = A3AP;
	}
#endif
}

/* Search direction with the max value in histogram */
int WorkMaxT::maxDirection(FPTYPE *histFourS1, int *maxITheta) {
	int nMaxTheta = 0;
	int nData = nx * ny;
	FPTYPE vals[NMAXTHETA + 2], maxVal1, maxVal2;
	int pos, maxPos1, maxPos2;
	pos = 0;
	for (int iTheta = 0 ; iTheta < NMAXTHETA ; ++iTheta) {
		FPTYPE val     = 0.0;
		int posCoef = 0;
		for (int order = 0 ; order < 2 * PANGLE + 1 ; ++order) {
			// printf("%d  %f\n", order, histFourS1[posCoef]);
			val += invFourTbl[pos++] * histFourS1[posCoef];
			posCoef += nData;
		}
		vals[iTheta + 1] = val;
		// printf("iTheta = %d val = %f \n", iTheta, vals[iTheta + 1]);
	}
	vals[0]             = vals[NMAXTHETA];
	vals[NMAXTHETA + 1] = vals[1];
	maxPos1 = maxPos2 = 0;
	for (int iThetaP = 1 ; iThetaP <= NMAXTHETA ; ++iThetaP) {
		if (vals[iThetaP] > vals[iThetaP - 1] && vals[iThetaP] > vals[iThetaP + 1]) { /* Local maximum */
			switch(nMaxTheta) {
			case 0:
				maxVal1 = vals[iThetaP];
				maxPos1 = iThetaP;
				nMaxTheta = 1;
				break;
			case 1:
				if (vals[iThetaP] > maxVal1) {
					maxVal2 = maxVal1;
					maxPos2 = maxPos1;
					maxVal1 = vals[iThetaP];
					maxPos1 = iThetaP;
				} else {
					maxVal2 = vals[iThetaP];
					maxPos2 = iThetaP;
				}
				nMaxTheta = 2;
				break;
			case 2:
				if (vals[iThetaP] > maxVal1) {
					maxVal2 = maxVal1;
					maxPos2 = maxPos1;
					maxVal1 = vals[iThetaP];
					maxPos1 = iThetaP;
				} else if (vals[iThetaP] > maxVal2) {
					maxVal2 = vals[iThetaP];
					maxPos2 = iThetaP;
				}
				break;
			}
		}
	}
	if (nMaxTheta == 2 && maxVal2 < maxVal1 * 0.7) nMaxTheta = 1;
	maxITheta[0] = maxPos1 - 1;
	maxITheta[1] = maxPos2 - 1;
	//printf("nMaxTheta = %d, maxITheta1 = %d maxTheta2 = %d  maxVal1 = %.1f maxVal2 = %.1f \n", nMaxTheta, maxITheta[0], maxITheta[1], maxVal1,  maxVal2 );
	return nMaxTheta;
}

WorkIIR::WorkIIR(int P, int maxNxy, FPTYPE *coefG, FPTYPE *coefDG) {
	/* Initialize image regions P = 4 */
	this->P = P;
  blurCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 2));
	diffCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 1));
	blurCoef[P + 1] = 0.0;
	FPTYPE sign = 1.0;
	for (int pos = 0; pos <= P ; ++pos) {
		blurCoef[pos]    = sign * coefG[pos];
		diffCoef[pos]    = sign * coefDG[pos];
		blurCoef[P + 1] += blurCoef[pos];
		sign *= -1.0;
	}

	cosL    = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 1));
	sinL    = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 1));
	lineAdd = (FPTYPE *) malloc(sizeof(FPTYPE) * (maxNxy + MAXK));
	lineSub = (FPTYPE *) malloc(sizeof(FPTYPE) * (maxNxy + MAXK));
	inteCos = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 1));
	inteSin = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 1));
	cosInitCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 1));
	sinInitCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (P + 1));
	}

int WorkIIR::initByK(int K0) {
	// printf("Init by K %d\n", K);
	if (K < 2) K = 2;
	K = K0;
	for (int p = 1 ; p <= P ; ++p) {
		cosL[p] = cos(PI * p / K);
		sinL[p] = sin(PI * p / K);
	}
  /* Initial coefficient of integral signal */
  for(int p = 1 ; p <= P ; ++p) {
    if (p % (K + K) == 0) {
      cosInitCoef[p] = K;
      sinInitCoef[p] = 0.0;
    } else if (p % 2 == 0) {
      cosInitCoef[p] = 0.0;
      sinInitCoef[p] = 0.0;
    } else {
      cosInitCoef[p] = 1.0;
      sinInitCoef[p] = sinL[p] / (1 - cosL[p]);
    }
  }
	return(0);
}

