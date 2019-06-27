/*
 * calFeature.c
 *
 *  Created on: 2018/03/29
 *      Author: yamasita
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "fastGMatch.h"
#include "calFeature.h"
#include "utilities.h"

/* Make feature vector */
int mkFeature(FPTYPE *histFourS4, WorkMT *workMT, Features *features, Feature *feature) {
	int ix, iy, ifx, ify, pos, order, posXY, posXY2, posHist;
	double sqSum, sqr;
	int *rePosSc;

  // printf("iTheta = %d\n", feature->iTheta);
	rePosSc = &(features->relativePosScL[feature->iTheta * GRIDSIZE * GRIDSIZE * 2]);

	posXY = posXY2 = 0;
	for (ify = 0 ; ify < GRIDSIZE ; ++ify) {
		for (ifx = 0 ; ifx < GRIDSIZE ; ++ifx) {
			//printf("(%d, %d)  ", rePosSc[posXY2], rePosSc[posXY2 + 1]);
			ix = feature->ix + rePosSc[posXY2++];
			iy = feature->iy + rePosSc[posXY2++];
			// printf("(%d, %d)\n", ix, iy);
			posHist = ix + features->nx * iy;
			features->tmpVector[0] = histFourS4[posHist];
			for (order = 1 ; order < features->nAngleCoef ; ++order) {
				posHist += features->nxy;
				features->tmpVector[order] = histFourS4[posHist];
			}
			rotateDirHist(features->tmpVector, &(feature->vector[features->nAngleCoef * (posXY++)]), feature->iTheta * features->nAngleCoef, workMT);
			//for (ix = 0 ; ix < features->nAngleCoef ; ++ix) printf("%f, ", feature->vector[features->nAngleCoef * (posXY - 1) + ix]);
			//printf("\n");
		}
	}
	/* Normalization */
	sqSum = 0.0;
	for (pos = 0 ; pos < features->nAngleCoef * GRIDSIZE * GRIDSIZE ; ++pos) 	sqSum += feature->vector[pos] * feature->vector[pos];
	sqr = sqrt(sqSum);
	for (pos = 0 ; pos < features->nAngleCoef * GRIDSIZE * GRIDSIZE ; ++pos) feature->vector[pos] /= sqr;
	return 0;
}

/* Calculate Fourier expression of directional histogram */
int calDirHistPoint(FPTYPE *diffXImg, FPTYPE *diffYImg, FPTYPE *dirHist, int nx, int ny) {
	int ixy, nxy = nx * ny, nxy2 = 2 * nxy, stCos, stSin, order;
	FPTYPE diffX, diffY, diffPower, diffSqrt, ansCos, ansSin, diffCos, diffSin;

	for (ixy = 0 ; ixy < nxy ; ++ixy) {
		diffX     = diffXImg[ixy];
		diffY     = diffYImg[ixy];
		diffPower = diffX * diffX + diffY * diffY;
#ifdef HISTREG
		dirHist[ixy] = diffSqrt = sqrt(diffPower + EDGEPOWER);
		stCos = nxy; stSin = nxy2;
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
#else
		if (diffPower > EDGEPOWER) {
			dirHist[ixy] = diffSqrt = sqrt(diffPower);
			stCos = nxy; stSin = nxy2;
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
		} else {
			dirHist[ixy] = 0.0;
			stCos = nxy; stSin = nxy2;
			for (order = 1 ; order <= PANGLE ; ++order) {
				dirHist[ixy + stCos] = 0.0;
				dirHist[ixy + stSin] = 0.0;
				stCos = stCos + nxy2, stSin = stSin + nxy2;
			}
		}
#endif
	}
	return 0;
}

/* Search direction with the max value in histogram */
int maxDirection(FPTYPE *histFourS1, WorkMT *workMT, int *maxITheta) {
	int nMaxTheta = 0, iTheta, iThetaP, order, posCoef;
	int nData = workMT->nx * workMT->ny;
	FPTYPE vals[NMAXTHETA + 2], val, maxVal1, maxVal2;
	int pos, maxPos1, maxPos2;
	pos = 0;
	for (iTheta = 0 ; iTheta < NMAXTHETA ; ++iTheta) {
		val     = 0.0;
		posCoef = 0;
		for (order = 0 ; order < 2 * PANGLE + 1 ; ++order) {
			val += workMT->invFourTbl[pos++] * histFourS1[posCoef];
			//printf("%d  %f\n", order, blurLargeDirHist[posCoef]);
			posCoef += nData;
		}
		vals[iTheta + 1] = val;
		// printf("iTheta = %d val = %f \n", iTheta, vals[iTheta + 1]);
	}
	vals[0]             = vals[NMAXTHETA];
	vals[NMAXTHETA + 1] = vals[1];
	maxPos1 = maxPos2 = 0;
	for (iThetaP = 1 ; iThetaP <= NMAXTHETA ; ++iThetaP) {
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
	// printf("nMaxTheta = %d, maxITheta1 = %d maxTheta2 = %d  maxVal1 = %.1f maxVal2 = %.1f \n", nMaxTheta, maxITheta[0], maxITheta[1], maxVal1,  maxVal2 );
	return nMaxTheta;
}

/* Rotation of directional histogram */
int rotateDirHist(FPTYPE *vectorIn, FPTYPE *vectorOut, int stRotTbl, WorkMT *workMT) {
	int order, posRotTbl, posVect;
	FPTYPE vCos, vSin, rCos, rSin;

	posRotTbl  = stRotTbl + 1;
	posVect    = 0;
	vectorOut[posVect] = vectorIn[posVect]; ++posVect;
	for (order = 1 ; order <= PANGLE ; ++order) {
		vCos = vectorIn[posVect];
		vSin = vectorIn[posVect + 1];
		rCos = workMT->rotTbl[posRotTbl++];
		rSin = workMT->rotTbl[posRotTbl++];
		vectorOut[posVect++] =  vCos * rCos - vSin * rSin;
		vectorOut[posVect++] =  vCos * rSin + vSin * rCos;
	}
	return 0;
}

/* Large scale smoothing is approximated by small scale smoothing */
int approxLargeScale(FPTYPE *histFourS4, FPTYPE *histFourAS1, int ix, int iy, WorkMT *workMT) {
	int pInd, pos, tPos, tPos2, order, orderP;
	int nxy = workMT->nx * workMT->ny;
	FPTYPE a;

	/* Initialize the output */
	for (order = 0 ; order < 2 * PANGLE + 1 ; ++order) 	histFourAS1[order] = 0.0;

	tPos = tPos2 = 0;
	for (pInd = 0 ; pInd < NAPPROPOINT ; ++pInd) {
		//printf("(%d, %d) (%d, %d) a = %f \n", ix, iy, ix + workMT->largeScaleTbl[tPos2], iy + workMT->largeScaleTbl[tPos2+1], workMT->largeScaleATbl[tPos]);
		pos  = (ix + workMT->largeScaleTbl[tPos2]);
		pos += workMT->nx * (iy + workMT->largeScaleTbl[tPos2++]);
		a   = workMT->largeScaleATbl[tPos++];
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
int gaussDiff(FPTYPE *inImg, FPTYPE *diffXImg,  FPTYPE *diffYImg, WorkImg *workImg, WorkIIR *workIIR) {
	/* X direction */;
	gaussSmooth1(2, inImg, workImg->xBlurImg, workImg->xDiffXImg, workImg->nx, workImg->ny,  workIIR);
	/* Y direction */
	transposition(workImg->xBlurImg, workImg->trans0, workImg->nx, workImg->ny);
	gaussSmooth1(1, workImg->trans0, workImg->trans1, workImg->trans2, workImg->ny, workImg->nx,  workIIR);
	transposition(workImg->trans2, diffYImg, workImg->ny, workImg->nx);

	transposition(workImg->xDiffXImg, workImg->trans0, workImg->nx, workImg->ny);
	gaussSmooth1(0, workImg->trans0, workImg->trans1, workImg->trans2, workImg->ny, workImg->nx,  workIIR);
	transposition(workImg->trans1, diffXImg, workImg->ny, workImg->nx);

	return 0;
}

/* Gauss smoothing */
int gaussSmooth(FPTYPE *inImg, FPTYPE *blurImg, WorkImg *workImg, WorkIIR *workIIR) {

	/* X direction */;
	gaussSmooth1(0, inImg, workImg->xBlurImg, workImg->xDiffXImg, workImg->nx, workImg->ny,  workIIR);

	/* Y direction */
	transposition(workImg->xBlurImg, workImg->trans0, workImg->nx, workImg->ny);
	gaussSmooth1(0, workImg->trans0, workImg->trans1, workImg->trans2, workImg->ny, workImg->nx,  workIIR);
	transposition(workImg->trans1, blurImg,  workImg->ny, workImg->nx);
	return 0;
}

/* One dimensional transform type = 0 : only blur, type = 1 : only diff, type = 2 : blur and diff */
int gaussSmooth1(int type, FPTYPE *inImg, FPTYPE *bImg, FPTYPE *dImg, int nx, int ny,  WorkIIR *workIIR) {
	int ix, iy, pos, posS, loopOrd;
	int K = workIIR->K, K2 = 2 * K;
	int nOrd = workIIR->nOrd;
	FPTYPE valB, valD;

	workIIR->nInt = nx + K2;
	pos = 0;
	for (iy = 0 ; iy < ny ; ++iy) {
		for (ix = 0 ; ix < nx ; ++ix) {
			workIIR->lineExt[ix + K] = inImg[ix + pos];
		}
		sincosIIRF(nx, workIIR);
		for (ix = 0 ; ix < nx ; ++ix) {
			posS = ix;

			if (type == 0) { /* Only blurring */
				valB = 0.0;
				for (loopOrd = 0 ; loopOrd <= nOrd ; ++loopOrd) {
					valB += workIIR->blurCoef[loopOrd] * (workIIR->intCos[posS + K2] - workIIR->intCos[posS]);
					posS += workIIR->nInt;
				}
				bImg[ix + pos] = valB + workIIR->blurCoef[nOrd + 1] * workIIR->lineExt[ix];
			} else	if (type == 1) { /* Only differential */
				 valD = 0.0;
				for (loopOrd = 1 ; loopOrd <= nOrd ; ++loopOrd) {
					posS += workIIR->nInt;
					valD += workIIR->diffCoef[loopOrd] * (workIIR->intSin[posS + K2] - workIIR->intSin[posS]);
				}
				dImg[ix + pos] = valD;
			} else {
				loopOrd = 0;
				valB = workIIR->blurCoef[loopOrd] * (workIIR->intCos[posS + K2] - workIIR->intCos[posS]);
				valD = 0.0;
				for (loopOrd = 1 ; loopOrd <= nOrd ; ++loopOrd) {
					posS += workIIR->nInt;
					valB += workIIR->blurCoef[loopOrd] * (workIIR->intCos[posS + K2] - workIIR->intCos[posS]);
					valD += workIIR->diffCoef[loopOrd] * (workIIR->intSin[posS + K2] - workIIR->intSin[posS]);
				}
				bImg[ix + pos] = valB + workIIR->blurCoef[nOrd + 1] * workIIR->lineExt[ix];
				dImg[ix + pos] = valD;
			}
		}
		pos += nx;
	}
	return 0;
}

/* Make sinusoidal integral image */
int sincosIIRF(int nData, WorkIIR *workIIR) {
	int pos, posInt, loopOrd;
	int K        = workIIR->K;
	FPTYPE *cosL = workIIR->cosL;
	FPTYPE *sinL = workIIR->sinL;

	switch (workIIR->extType) {
	case 0:
		for (pos = 0         ; pos < K             ; ++pos) workIIR->lineExt[pos] = 0.0;
		for (pos = K + nData ; pos < workIIR->nInt ; ++pos) workIIR->lineExt[pos] = 0.0;
		break;
	case 1:
		for (pos = 0         ; pos < K             ; ++pos) workIIR->lineExt[pos] = workIIR->lineExt[K];
		for (pos = K + nData ; pos < workIIR->nInt ; ++pos) workIIR->lineExt[pos] = workIIR->lineExt[K + nData - 1];
		break;
	}

	posInt = 0;
	workIIR->intCos[posInt] = workIIR->lineExt[0];
	for (loopOrd = 0 ; loopOrd < workIIR->nOrd ; ++loopOrd) {
		posInt += workIIR->nInt;
		workIIR->intCos[posInt] = workIIR->lineExt[0];
		workIIR->intSin[posInt] = 0.0;
	}

	for (pos = 0 ; pos < workIIR->nInt - 1 ; ++pos) {
		posInt = pos;
		workIIR->intCos[posInt + 1] = workIIR->intCos[posInt] + workIIR->lineExt[pos + 1]; 	/* For 0 th order */
		for (loopOrd = 0 ; loopOrd < workIIR->nOrd ; ++loopOrd) {
			posInt += workIIR->nInt;
			workIIR->intCos[posInt + 1] = cosL[loopOrd] * workIIR->intCos[posInt] - sinL[loopOrd] * workIIR->intSin[posInt] + workIIR->lineExt[pos + 1];
			workIIR->intSin[posInt + 1] = sinL[loopOrd] * workIIR->intCos[posInt] + cosL[loopOrd] * workIIR->intSin[posInt];
		}
	}
	return 0;
}

/* Transposition of an image */
int transposition(FPTYPE *inImg, FPTYPE *outImg, int nx, int ny) {
	int ix, iy, posT, pos = 0;
	for (iy = 0 ; iy < ny ; ++iy) {
		posT = iy;
		for (ix = 0 ; ix < nx ; ++ix) {
			outImg[posT] = inImg[pos];
			posT += ny; ++pos;
		}
	}
	return 0;
}
