/*
 * fastGMatch.c
 *
 *  Created on: 2018/03/29
 *      Author: yamasita
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "fastGMatch.h"
#include "calFeature.h"
#include "utilities.h"
#include <time.h>

int init(int nx, int ny, int nAngleCoef, int mkTemplate, Features *features,  Features *featuresTmpl,
		WorkImg *workImg, WorkIIR *workIIR2, WorkIIR *workIIR4, WorkMT *workMT, GptTbl *gptTbl);
int initByK(int K, WorkIIR *workIIR);

int main(int argc, char *argv[]) {
	int nx = NX, ny = NY;
	int nAngleCoef = 2 * PANGLE + 1;
	int K2, K4;
	int step, order, pos, stPos, stX, stY, endX, endY;
	int ix, iy, iFeature, transType = 0;
	int nxy = nx * ny;
	int nMaxTheta;
	int maxITheta[2];
	int mkTemplate = 0, nTrans;
	double scale, scaleRatio, scaleMax;
	char fileName[MAX_FILENAME];
	char dFileName[MAX_FILENAME];     /* Data file name*/
	char tmplDFileName[MAX_FILENAME];  /* template data file name */
	int printTime = 1;

	clock_t start, now;

	int    *inImgI         = (int *) malloc(sizeof(int) * nxy);
	FPTYPE *inImgOrg       = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	FPTYPE *inImg          = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	FPTYPE *diffXImg       = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	FPTYPE *diffYImg       = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	FPTYPE *dirHist        = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy * nAngleCoef);
	FPTYPE *histFourS1     = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy * nAngleCoef);
	FPTYPE *histFourS4     = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy * nAngleCoef);
	Features *features     = (Features *) malloc(sizeof(Features));
	Features *featuresTmpl = (Features *) malloc(sizeof(Features));
	Features *featuresP;
	WorkIIR *workIIR2      = (WorkIIR *) malloc(sizeof(WorkIIR));
	WorkIIR *workIIR4      = (WorkIIR *) malloc(sizeof(WorkIIR));
	WorkImg *workImg       = (WorkImg *) malloc(sizeof(WorkImg));
	WorkMT  *workMT        = (WorkMT *) malloc(sizeof(WorkMT));
	GptTbl  *gptTbl        = (GptTbl *) malloc(sizeof(GptTbl));

	/* Set file names */
	strcpy(fileName, WORKBASE); strcat(fileName, IMGIN);
	strcpy(dFileName, WORKBASE); strcat(dFileName, IMGOUT);
	strcpy(tmplDFileName, WORKBASE); strcat(tmplDFileName, IMGOUT);
	if (argc == 4) {
		mkTemplate = atoi(argv[1]);
		strcat(fileName, argv[2]);
		strcat(dFileName, argv[2]);
		strcat(tmplDFileName, argv[3]);
	} else {
		mkTemplate = MKTEMPLATE;
		strcat(fileName, IMAGENAME);
		strcat(dFileName, IMAGENAME);
		strcat(tmplDFileName, TEMPLATEDNAME);
	}
	strcat(fileName, ".pgm"); 	strcat(tmplDFileName, ".dat"); strcat(dFileName, ".dat");

	init(nx, ny, nAngleCoef, mkTemplate, features, featuresTmpl, workImg, workIIR2, workIIR4, workMT, gptTbl);

	/* Input image */
	loadImageFile(fileName, inImgI, nx, ny);
	for (ix = 0 ; ix < nxy ; ++ix) inImg[ix] = inImgOrg[ix] = inImgI[ix];

	if (mkTemplate) {
		featuresP           = featuresTmpl;
		featuresP->nFeature = 0;
		nTrans              = gptTbl->nTran;
		scaleMax            = TEMPLSCALEMAX;
		scaleRatio          = TEMPLSCALERATIO;
	} else {
		featuresP           = features;
		featuresP->nFeature = 0;
		nTrans              = 1;
		scaleMax            = SCALEMAX;
		scaleRatio          = SCALERATIO;
	}
	printf("nTrans = %d  %d \n", nTrans, mkTemplate);
	start = clock();
	for (transType = 0 ; transType < nTrans ; ++transType)	{
		printf("nTrans = %d  %d \n", nTrans, transType);
		if (mkTemplate) {
			gptTransformImage(& (gptTbl->gptTbl[transType * 9]), inImgOrg, inImg, nx, ny, TEMPLATEX, TEMPLATEY);
			sprintf(fileName, "%s%stranImg%02d.pgm", WORKBASE, IMGOUT, transType); saveImageFileFp(fileName, inImg, nx, ny);
			scale = TEMPLSCALEINIT;
		} else {
			scale = SCALEINIT;
		}
		while(scale <= scaleMax) {
			/* Calculate gradient images */
			K4 = (int) (PI * scale * GRADRATIO / (GRIDSIZE * SIGMA4) + ROUNDFRAC);
			workIIR4->extType = 1;
			initByK(K4, workIIR4);
			gaussDiff(inImg, diffXImg, diffYImg, workImg, workIIR4);
			if (printTime == 1) printf("End of diff image %lf (sec)\n", (double)((now = clock())- start)/CLOCKS_PER_SEC);

			/* Output differential image */
			strcpy(fileName, WORKBASE); strcat(fileName, IMGOUT); strcat(fileName, "diffX.pgm"); saveImageFileFp(fileName, diffXImg, nx, ny);
			strcpy(fileName, WORKBASE); strcat(fileName, IMGOUT); strcat(fileName, "diffY.pgm"); saveImageFileFp(fileName, diffYImg, nx, ny);

			/* Fourier series expression of directions */
			calDirHistPoint(diffXImg, diffYImg, dirHist, nx, ny);

#ifdef USELARGESCALE
			/* Blurred directions for large scale */
			K2 = (int) (PI * scale / SIGMA2 + ROUNDFRAC);
			workIIR2->extType = 1;
			initByK(K2, workIIR2);
			stPos = 0;
			for (order = 0 ; order < 2 * PANGLE + 1 ; ++order) {
				gaussSmooth(&(dirHist[stPos]), &(histFourS1[stPos]), workImg, workIIR2);
				//printf("order %d  %f   %f\n", order, dirHist[stPos + nx / 2 + nx * ny /2], blurLargeDirHist[stPos + nx / 2 + nx * ny /2] );
				stPos += nxy;
			}
			if (printTime == 1) printf("End of large blur %lf (sec)\n", (double)((now = clock()) - start)/CLOCKS_PER_SEC);
#else
			for (pos = 0 ; pos < 2 * NAPPROPOINT ; ++pos) {
				workMT->largeScaleTbl[pos] = (int) (scale * workMT->largeScaleRelTbl[pos] + ROUNDFRAC);
			}
#endif

			/* Blurred directions for small scale */
			K2 = (int) (PI * scale / (SIGMA2 * GRIDSIZE) + ROUNDFRAC);
			workIIR2->extType = 1;
			initByK(K2, workIIR2);
			stPos = 0;
			for (order = 0 ; order < 2 * PANGLE + 1; ++order) {
				gaussSmooth(&(dirHist[stPos]), &(histFourS4[stPos]), workImg, workIIR2);
				stPos += nxy;
			}
			if (printTime == 1) printf("End of small blur %lf (sec)\n", (double)((now = clock()) - start)/CLOCKS_PER_SEC);

			for (pos = 0 ; pos < NMAXTHETA * GRIDSIZE * GRIDSIZE * 2 ; ++pos) {
				featuresP->relativePosScL[pos] = (int) (scale * featuresP->relativePosL[pos] + ROUNDFRAC);
				// printf("relativePosScL[%d] = %d\n", pos, featuresP->relativePosScL[pos]);
			}
			/* Process for each sample point */
			step  = (int) (scale * STEPRATIO + ROUNDFRAC);

			if (mkTemplate) {
				stX = endX = TEMPLATEX; stY = endY = TEMPLATEY;
			} else {
				stX = stY = (int) (1.5 * scale + ROUNDFRAC);
				endX = nx - stX; endY = ny - stY;
			}

			for (iy = stY ; iy <= endY ; iy += step) {
				for (ix = stX ; ix <= endX ; ix += step) {
					// printf("(ix, iy) = (%d, %d)\n", ix, iy);
					pos = ix + nx * iy;
					/* Local direction using first mode */
#ifndef USELARGESCALE
					approxLargeScale(histFourS4, & (histFourS1[pos]), ix, iy, workMT);
#endif
					nMaxTheta = maxDirection(& (histFourS1[pos]), workMT, maxITheta);

					if (mkTemplate == 1 && nMaxTheta >= 2) nMaxTheta = 1; /* At most one feature for template */
					/* Rotate and construct feature vectors  */
					//printf("theta = %f, maxTheta = %f  %d  %f\n", 2.0 * PI * transType / NTEMPLATEROT, 2.0 * PI * maxITheta[0] / NMAXTHETA, nMaxTheta, 2.0 * PI * maxITheta[1] / NMAXTHETA);
					//printf("%f  %f \n", 2.0 * PI * transType / NTEMPLATEROT, 2.0 * PI * maxITheta[0] / NMAXTHETA);
					for (iFeature = 0 ; iFeature < nMaxTheta ; ++iFeature) {
						Feature *feature = & ((featuresP->featureL)[featuresP->nFeature]);
						feature->ix        = ix;
						feature->iy        = iy;
						feature->iTheta    = maxITheta[iFeature];
						feature->ordHist   = iFeature;
						feature->scale     = scale;
						feature->transType = transType;
						mkFeature(histFourS4, workMT, featuresP, feature);
						++(featuresP->nFeature);
					}
				}
			}
			//printf("scale = %6.2f, step = %d, nFeature = %d \n", scale, step, featuresP->nFeature);
			scale = scale * scaleRatio;
			if (printTime == 1) printf("End of cal feature %lf (sec)\n", (double)((now = clock()) - start)/CLOCKS_PER_SEC);
		} /* End of scale loop */
	} /* End of transformation loop */

	if (mkTemplate) {
		saveFeatures(tmplDFileName, featuresTmpl);
		printf("Output %s and end of make tmplate nFeature = %d \n", tmplDFileName, featuresTmpl->nFeature);
	} else {
		saveFeatures(dFileName, features);
		loadFeatures(tmplDFileName, featuresTmpl); /* Read features of template*/
		/* Calc matching */
		int iFeatureTmpl, iVec;
		FPTYPE *record = (FPTYPE *) malloc(sizeof(FPTYPE) * features->nFeature * featuresTmpl->nFeature);
		int    *index = (int *)     malloc(sizeof(int)     * features->nFeature * featuresTmpl->nFeature);
		pos = 0;
		printf("nFeature %d  %d \n",featuresTmpl->nFeature, features->nFeature);
		for (iFeature = 0 ; iFeature < features->nFeature ; ++iFeature) {
			FPTYPE *vector = (features->featureL)[iFeature].vector;
			for (iFeatureTmpl = 0 ; iFeatureTmpl < featuresTmpl->nFeature ; ++iFeatureTmpl) {
				FPTYPE *vectorTmpl = (featuresTmpl->featureL)[iFeatureTmpl].vector;
				FPTYPE val = 0.0;
				for (iVec = 0 ; iVec < (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE ; ++iVec) {
					val += vector[iVec] * vectorTmpl[iVec];
				}
				index[pos]  = pos;
				record[pos] = val;
				++pos;
			}
		}
		heapSort(record, index,  features->nFeature * featuresTmpl->nFeature);
		iFeature     = index[0] / featuresTmpl->nFeature;
		iFeatureTmpl = index[0] % featuresTmpl->nFeature;
		printf("Corelation = %f,  femplate No. = %2d,  ", record[0], iFeatureTmpl);
		prFeatureInf(& (features->featureL[iFeature]));
		printf("Output %s and end of matching \n", dFileName);
	}
	now = clock(); printf("End of calculation %lf (sec)\n", (double)(now - start)/CLOCKS_PER_SEC);
	return 0;
}

int init(int nx, int ny, int nAngleCoef, int mkTemplate, Features *features, Features *featuresTmpl,
		WorkImg *workImg, WorkIIR *workIIR2, WorkIIR *workIIR4, WorkMT *workMT, GptTbl *gptTbl) {
	int nxy = nx * ny;
	int maxNxy, pos, posTheta, order, iTheta, ix, iy;
	double dTheta, cosTheta, sinTheta, x, y, scale;
	FPTYPE sign;

	/* Coefficients P = 2, sigma = 1.1 */
	FPTYPE coefG2[]   = {1.5846315202e-01, 1.7508079293e-01, 2.7323762984e-02};
	FPTYPE coefDG2[]  = {0.0000000000e+00, -1.7506270969e-01, -5.4683692457e-02};

	/* Coefficients P = 4, sigma = 0.8 */
	FPTYPE coefG4[]   = {1.5914081930e-01,  2.3116780826e-01,  8.8476942331e-02,  1.7890179087e-02,  1.8836364529e-03};
	FPTYPE coefDG4[]  = {0.0000000000e+00, -2.3116693143e-01, -1.7695563833e-01, -5.3667906760e-02, -7.5380531476e-03};

	features->relativePosL   = (FPTYPE *) malloc(sizeof(FPTYPE) * NMAXTHETA * GRIDSIZE * GRIDSIZE * 2);
	features->relativePosScL = (int *)    malloc(sizeof(int) * NMAXTHETA * GRIDSIZE * GRIDSIZE * 2);
	features->tmpVector      = (FPTYPE *) malloc(sizeof(FPTYPE) * (2 * PANGLE + 1) * 4 * 4);
	features->nx             = nx;
	features->ny             = ny;
	features->nxy            = nxy;
	features->nAngleCoef     = nAngleCoef;
	features->nEstFeature    = 0;

	featuresTmpl->relativePosL   = (FPTYPE *) malloc(sizeof(FPTYPE) * NMAXTHETA * GRIDSIZE * GRIDSIZE * 2);
	featuresTmpl->relativePosScL = (int *) malloc(sizeof(int) * NMAXTHETA * GRIDSIZE * GRIDSIZE * 2);
	featuresTmpl->tmpVector      = (FPTYPE *) malloc(sizeof(FPTYPE) * (2 * PANGLE + 1) * 4 * 4);
	featuresTmpl->nx             = nx;
	featuresTmpl->ny             = ny;
	featuresTmpl->nxy            = nxy;
	featuresTmpl->nAngleCoef     = nAngleCoef;
	featuresTmpl->nEstFeature    = 0;

	scale = TEMPLSCALEINIT;
	while(scale <= TEMPLSCALEMAX) {
		featuresTmpl->nEstFeature += NTEMPLATEROT * NTEMPLATEENLONG + 1;
		scale = (int) scale * TEMPLSCALERATIO;
	}
	featuresTmpl->featureL = (Feature *) malloc(sizeof(Feature) * featuresTmpl->nEstFeature);
	for (pos = 0 ; pos < featuresTmpl->nEstFeature ; ++pos) {
		featuresTmpl->featureL[pos].vector = (FPTYPE *) malloc(sizeof(FPTYPE) * (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE);
	}

	scale = SCALEINIT;
	while(scale <= SCALEMAX) {
		features->nEstFeature += (int) (nxy / (scale * STEPRATIO * scale * STEPRATIO) + ROUNDFRAC);
		scale = (int) scale * SCALERATIO;
	}
	features->nEstFeature *= 2;
	// printf("nEstFeature %d \n", features->nEstFeature);
	features->featureL = (Feature *) malloc(sizeof(Feature) * features->nEstFeature);
	for (pos = 0 ; pos < features->nEstFeature ; ++pos) {
		features->featureL[pos].vector = (FPTYPE *) malloc(sizeof(FPTYPE) * (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE);
	}

	/* Make feature of template with various transformation */
if (mkTemplate) {
#ifdef ONLYROT /* Only rotatation */
	FPTYPE cosTran, sinTran;
	gptTbl->nTran   = NTEMPLATEROT;
	gptTbl->gptTbl  = (FPTYPE *) malloc(sizeof(FPTYPE) * 9 * gptTbl->nTran);
	pos = 0;
	for (iTheta = 0 ; iTheta < NTEMPLATEROT ; ++iTheta) {
		cosTran = cos((2.0 * PI / NTEMPLATEROT) * iTheta); 		sinTran = sin((2.0 * PI / NTEMPLATEROT) * iTheta);
		gptTbl->gptTbl[C11 + pos] = cosTran;
		gptTbl->gptTbl[C22 + pos] = cosTran;
		gptTbl->gptTbl[C12 + pos] = - sinTran; /* Note that upside down */
		gptTbl->gptTbl[C21 + pos] = sinTran;
		gptTbl->gptTbl[C13 + pos] = gptTbl->gptTbl[C23 + pos] = gptTbl->gptTbl[C31 + pos]  = gptTbl->gptTbl[C32 + pos] = 0.0;
		gptTbl->gptTbl[C33 + pos] = 1.0;
		//gptPr(&(gptTbl->gptTbl[pos]), "GPTInit ");
		pos += 9;
	}
#else
	int iEnlong;
	FPTYPE cosTran, sinTran, alpha, beta;
	gptTbl->nTran   = NTEMPLATEROT * NTEMPLATEENLONG + 1;
	gptTbl->gptTbl  = (FPTYPE *) malloc(sizeof(FPTYPE) * 9 * gptTbl->nTran);
	gptTbl->gptTbl[C11] = gptTbl->gptTbl[C22] = 1.0;
	gptTbl->gptTbl[C12] = gptTbl->gptTbl[C21] = 0.0;
	gptTbl->gptTbl[C11] = gptTbl->gptTbl[C22] = 1.0;
	gptTbl->gptTbl[C13] = gptTbl->gptTbl[C23] = gptTbl->gptTbl[C31]  = gptTbl->gptTbl[C32] = 0.0;
	gptTbl->gptTbl[C33] = 1.0;
	pos = 0;
	for (iTheta = 0 ; iTheta < NTEMPLATEROT ; ++iTheta) {
		cosTran = cos((PI / NTEMPLATEROT) * iTheta); 		sinTran = sin((PI / NTEMPLATEROT) * iTheta);
		for (iEnlong = 1 ; iEnlong <= NTEMPLATEENLONG ; ++iEnlong) {
			pos += 9;
			alpha = 1.0 + TEMPLATEENCOEF * iEnlong; beta = 1.0 - TEMPLATEENCOEF * iEnlong;
			//printf("alpha %f  beta %f\n", alpha, beta);
			gptTbl->gptTbl[C11 + pos] = alpha * cosTran * cosTran + beta  * sinTran * sinTran;
			gptTbl->gptTbl[C22 + pos] = beta  * cosTran * cosTran + alpha * sinTran * sinTran;
			gptTbl->gptTbl[C12 + pos] = gptTbl->gptTbl[C21 + pos] = (alpha - beta) * cosTran * sinTran;
			gptTbl->gptTbl[C13 + pos] = gptTbl->gptTbl[C23 + pos] = gptTbl->gptTbl[C31 + pos]  = gptTbl->gptTbl[C32 + pos] = 0.0;
			gptTbl->gptTbl[C33 + pos] = 1.0;
			//gptPr(&(gptTbl[pos]), "GPTInit ");
		}
	}
#endif
}

	/* Table for grid relative position */
	pos = 0;
	for (iTheta = 0 ; iTheta < NMAXTHETA ; ++iTheta) {
		cosTheta = cos(2 * PI * iTheta / NMAXTHETA);
		sinTheta = sin(2 * PI * iTheta / NMAXTHETA);
		y = - (GRIDSIZE - 1.0) / GRIDSIZE;
		for (iy = 0 ; iy < GRIDSIZE ; ++iy) {
			x = - (GRIDSIZE - 1.0) / GRIDSIZE;
			for (ix = 0 ; ix < GRIDSIZE ; ++ix) {
				featuresTmpl->relativePosL[pos] = features->relativePosL[pos] = x * cosTheta - y * sinTheta; ++pos;// x
				featuresTmpl->relativePosL[pos] = features->relativePosL[pos] = x * sinTheta + y * cosTheta; ++pos;// y
				x = x + 2.0 / GRIDSIZE;
			}
			y = y + 2.0 / GRIDSIZE;
		}
	}

	if (nx > ny) maxNxy = nx;
	else maxNxy = ny;

	/* Initialize image regions */
	workImg->nx        = nx;
	workImg->ny        = ny;
	workImg->xBlurImg  = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	workImg->xDiffXImg = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	workImg->trans0    = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	workImg->trans1    = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	workImg->trans2    = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);

	/* Initialize image regions P = 2 */
	workIIR2->nOrd = 2;
	workIIR2->blurCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR2->nOrd + 1));
	workIIR2->diffCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR2->nOrd + 1));
	workIIR2->blurCoef[workIIR2->nOrd + 1] = 0.0;
	sign = 1.0;
	for (pos = 0; pos <= workIIR2->nOrd ; ++pos) {
		workIIR2->blurCoef[pos]                 = sign * coefG2[pos];
		workIIR2->diffCoef[pos]                 = sign * coefDG2[pos];
		printf("diff Coef %d = %lf\n", pos, workIIR2->diffCoef[pos]);
		workIIR2->blurCoef[workIIR2->nOrd + 1] += workIIR2->blurCoef[pos];
		sign *= -1;
	}
	workIIR2->maxNInt = maxNxy + 2 * MAXK;
	workIIR2->cosL    = (FPTYPE *) malloc(sizeof(FPTYPE) * workIIR2->nOrd);
	workIIR2->sinL    = (FPTYPE *) malloc(sizeof(FPTYPE) * workIIR2->nOrd);
	workIIR2->lineExt = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR2->maxNInt));
	workIIR2->intCos  = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR2->maxNInt) * (workIIR2->nOrd + 2));
	workIIR2->intSin  = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR2->maxNInt) * (workIIR2->nOrd + 1));
	workIIR2->intCos[0] = 0.0;
	workIIR2->intSin[0] = 0.0;

	/* Initialize image regions P = 4 */
	workIIR4->nOrd = 4;
	workIIR4->blurCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR4->nOrd + 2));
	workIIR4->diffCoef = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR4->nOrd + 1));
	workIIR4->blurCoef[workIIR4->nOrd + 1] = 0.0;
	sign = 1.0;
	for (pos = 0; pos <= workIIR4->nOrd ; ++pos) {
		workIIR4->blurCoef[pos]                 = sign * coefG4[pos];
		workIIR4->diffCoef[pos]                 = sign * coefDG4[pos];
		workIIR4->blurCoef[workIIR4->nOrd + 1] += workIIR4->blurCoef[pos];
		sign *= -1;
	}
	workIIR4->maxNInt = maxNxy + 2 * MAXK;
	workIIR4->cosL    = (FPTYPE *) malloc(sizeof(FPTYPE) * workIIR4->nOrd);
	workIIR4->sinL    = (FPTYPE *) malloc(sizeof(FPTYPE) * workIIR4->nOrd);
	workIIR4->lineExt = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR4->maxNInt));
	workIIR4->intCos  = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR4->maxNInt) * (workIIR4->nOrd + 2));
	workIIR4->intSin  = (FPTYPE *) malloc(sizeof(FPTYPE) * (workIIR4->maxNInt) * (workIIR4->nOrd + 1));
	workIIR4->intCos[0] = 0.0;
	workIIR2->intSin[0] = 0.0;

	/* Initialize workMT */
	workMT->nx               = nx;
	workMT->ny               = ny;
	workMT->invFourTbl       = (FPTYPE *) malloc(sizeof(FPTYPE) * 2 * NMAXTHETA * (2 * PANGLE + 1));
	workMT->rotTbl           = (FPTYPE *) malloc(sizeof(FPTYPE) * 2 * NMAXTHETA * (2 * PANGLE));
	workMT->largeScaleRelTbl = (FPTYPE *) malloc(sizeof(FPTYPE) * 2 * NAPPROPOINT);
	workMT->largeScaleTbl    = (int *)    malloc(sizeof(int)    * 2 * NAPPROPOINT);
	workMT->largeScaleATbl   = (FPTYPE *) malloc(sizeof(FPTYPE) * NAPPROPOINT);
	dTheta = 2 * PI / NMAXTHETA;

	/* For inverse Fourier transform of max angle*/
	pos = 0;
	for (posTheta = 0 ; posTheta < NMAXTHETA ; ++posTheta) {
		workMT->invFourTbl[pos++] = 0.5 / PI;
		for (order = 1 ; order <= PANGLE ; ++order) {
			workMT->invFourTbl[pos++] =   cos(dTheta * order * posTheta) / PI;
			workMT->invFourTbl[pos++] = - sin(dTheta * order * posTheta) / PI;
		}
	}

	/* For rotation of direction histogram */
	pos = 0;
	for (posTheta = 0 ; posTheta <= NMAXTHETA ; ++posTheta) {
		workMT->rotTbl[pos++] = 1.0;
		for (order = 1 ; order <= PANGLE ; ++order) {
			workMT->rotTbl[pos++] = cos(dTheta * order * posTheta);
			workMT->rotTbl[pos++] = sin(dTheta * order * posTheta);
		}
	}

	/* For approximation of Gaussian function of large variance by small variance */
#ifndef USELARGESCALE
	pos = order = 0;
	workMT->largeScaleRelTbl[pos++] = 0;
	workMT->largeScaleRelTbl[pos++] = 0;
	workMT->largeScaleATbl[order++] = A0AP;
	dTheta = 2 * PI / N1AP;
	for (posTheta = 0 ; posTheta < N1AP ; ++posTheta) {
		workMT->largeScaleRelTbl[pos++] = R1AP * cos(dTheta * posTheta);
		workMT->largeScaleRelTbl[pos++] = R1AP * sin(dTheta * posTheta);
		workMT->largeScaleATbl[order++] = A1AP;
	}
	dTheta = 2 * PI / N2AP;
	for (posTheta = 0 ; posTheta < N2AP ; ++posTheta) {
		workMT->largeScaleRelTbl[pos++] = R2AP * cos(dTheta * posTheta);
		workMT->largeScaleRelTbl[pos++] = R2AP * sin(dTheta * posTheta);
		workMT->largeScaleATbl[order++] = A2AP;
	}
	dTheta = 2 * PI / N3AP;
	for (posTheta = 0 ; posTheta < N3AP ; ++posTheta) {
		workMT->largeScaleRelTbl[pos++] = R3AP * cos(dTheta * posTheta);
		workMT->largeScaleRelTbl[pos++] = R3AP * sin(dTheta * posTheta);
		workMT->largeScaleATbl[order++] = A3AP;
	}
#endif
	return(0);
}

int initByK(int K, WorkIIR *workIIR) {
	// printf("Init by K %d\n", K);
	int order;
	if (K < 2) K = 2;
	workIIR->K = K;
	for (order = 1 ; order <= workIIR->nOrd ; ++order) {
		workIIR->cosL[order - 1] = cos(PI * order / K);
		workIIR->sinL[order - 1] = sin(PI * order / K);
	}
	return(0);
}
