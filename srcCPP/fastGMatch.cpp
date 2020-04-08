/*
 * fastGMatch.cpp
 *
 *  Created on: 2019/03/13
 *      Author: yamasita
 */


#include <cmath>
#include <string>
#include <iostream>

#include "parameter.h"
#include "utilities.h"
#include "extFeature.h"
#include "fastGMatch.h"

int main(int argc, char* argv[]) {
	int nx = NX, ny = NY, mkTemplate;

	if (argc == 4) {
		mkTemplate = atoi(argv[1]);
	} else {
		mkTemplate = MKTEMPLATE;
	}
	FastGMatch *fm = new FastGMatch(nx, ny, mkTemplate);

	fm->fileName      = fm->fileName      + WORKBASE + IMGPRE + BLURPRE + INPRE;
	fm->dFileName     = fm->dFileName     + WORKBASE + IMGPRE + BLURPRE + OUTPRE;
	fm->tmplDFileName = fm->tmplDFileName + WORKBASE + IMGPRE + TMPLPRE;
	/* Set file names */

	if (argc == 4) {
		fm->fileName      = fm->fileName      + argv[2] + ".pgm";
		fm->dFileName     = fm->dFileName     + argv[2] + ".dat";
		fm->tmplDFileName = fm->tmplDFileName + argv[3] + ".dat";
	} else {
		fm->fileName      = fm->fileName      + IMAGENAME + ".pgm";
		fm->dFileName     = fm->dFileName     + IMAGENAME + ".dat";
		fm->tmplDFileName = fm->tmplDFileName + TEMPLATEDNAME + ".dat";
	}


	/* Excute matching */
	fm->fastGMatch();
	return 0;
}

int FastGMatch::fastGMatch() {
	char fileNameC[MAX_FILENAME];

	if (mkTemplate) {
		featuresP           = featuresTmpl;
		nFeature            = 0;
		nTrans              = gptTbl->nTran;
		scaleMax            = TEMPLSCALEMAX;
		scaleRatio          = TEMPLSCALERATIO;
	} else {
		featuresP           = features;
		nFeature            = 0;
		nTrans              = 1;
		scaleMax            = SCALEMAX;
		scaleRatio          = SCALERATIO;
	}

	/* Input image */
	printf("file %s \n", fileName.c_str());
	imageIO->loadImageFile(fileName.c_str(), inImgI);
#ifdef INIMAGEDEBUG
	for (ix = 0 ; ix < nx ; ++ix) {
		for (iy = 0 ; iy < ny ; ++iy) {
			if (ix < 300) inImgI[ix + nx * iy] = 0;
			else inImgI[ix + nx * iy] = 128;
		}
	}
#endif

	for (ix = 0 ; ix < nxy ; ++ix) inImg[ix] = inImgOrg[ix] = inImgI[ix];
	printf("nTrans = %d mkTemplate = %d \n", nTrans, mkTemplate);

	start = clock();
	nFeature = 0;
	for (transType = 0 ; transType < nTrans ; ++transType)	{
		/* Calculate transformation of template */
		if (mkTemplate) {
			gptTbl->gptPara[transType]->gptTransformImage(inImgOrg, inImg, nx, ny, TEMPLATEX, TEMPLATEY);
			sprintf(fileNameC, "%s%s%stranImg%02d.pgm", WORKBASE, IMGPRE, TMPLPRE, transType);
			imageIO->saveImageFileFp(fileNameC, inImg);
			scale = TEMPLSCALEINIT;
		} else {
			scale = SCALEINIT;
		}

		while(scale <= scaleMax) {
			printf("nTrans = %d transType = %d scale = %f \n", nTrans, transType, scale);
			/* Process for each sample point */
			step  = (int) (scale * STEPRATIO + ROUNDFRAC);
			if (mkTemplate) {
				stX = endX = TEMPLATEX; stY = endY = TEMPLATEY;
			} else {
//				stX = stY = (int) (1. * scale + ROUNDFRAC);
				stX = stY = (int) (1.1 * scale + ROUNDFRAC);
				endX = nx - stX; endY = ny - stY;
			}
			//printf("st (%d, %d) end (%d, %d)", stX, stY, endX, endY);
			int cont = 0;
			for (iy = stY ; iy <= endY ; iy += step) {
				for (ix = stX ; ix <= endX ; ix += step) {
	      //printf("(ix, iy) = (%d, %d)\n", ix, iy);
					Feature *feature1, *feature2;
					nMaxTheta = extFeature->extract(cont, inImg, ix, iy, scale, & feature1, & feature2);
					if (nMaxTheta >= 1) {
						feature1->transType = transType;
						features[nFeature] = feature1;
						++nFeature;
					}
					if (mkTemplate == 0 && nMaxTheta >= 2) {
						feature2->transType = transType;
						features[nFeature] = feature2;
						++nFeature;
					}
					cont = 1;
				}
			}
			// printf("scale = %6.2f, step = %d, nFeature = %d \n", scale, step, featuresP->nFeature);
			scale = scale * scaleRatio;
			if (printTime == 1) printf("End of cal feature %lf (sec)\n", (double)((now = clock()) - start)/CLOCKS_PER_SEC);
		} /* End of scale loop */
	} /* End of transformation loop for make template */

	if (mkTemplate) {
		extFeature->saveFeatures(tmplDFileName.c_str(), nFeature, nAngleCoef, features);
		printf("Output %s and end of make template nFeature = %d \n", tmplDFileName.c_str(), nFeature);

	} else {
		int nFeatureTmpl, nAngleCoefTmpl;
		extFeature->saveFeatures(dFileName.c_str(), nFeature, nAngleCoef, features);
		extFeature->loadFeatures(tmplDFileName.c_str(), & nFeatureTmpl, & nAngleCoefTmpl, featuresTmpl); /* Read features of template*/

		/* Calc matching */
		int iFeatureTmpl, iVec;
		FPTYPE *record = (FPTYPE *) malloc(sizeof(FPTYPE) * nFeatureTmpl * nFeature);
		int    *index = (int *)     malloc(sizeof(int)    * nFeatureTmpl * nFeature);
		// printf("nFeature %d  %d \n",featuresTmpl->nFeature, features->nFeature);

		pos = 0;
		for (iFeature = 0 ; iFeature < nFeature ; ++iFeature) {
			FPTYPE *vector = features[iFeature]->vector;
			for (iFeatureTmpl = 0 ; iFeatureTmpl < nFeatureTmpl ; ++iFeatureTmpl) {
				FPTYPE *vectorTmpl = featuresTmpl[iFeatureTmpl]->vector;
				FPTYPE val = 0.0;
				for (iVec = 0 ; iVec < (2 * PANGLE + 1) * GRIDSIZE * GRIDSIZE ; ++iVec) {
					//printf("iF %d iFT %d iV %d v %f vT %f\n", iFeature, iFeatureTmpl, iVec, vector[iVec], vectorTmpl[iVec]);
					val += vector[iVec] * vectorTmpl[iVec];
				}
				index[pos]  = pos;
				record[pos] = val;
				++pos;
			}
		}
		utilities->heapSort(record, index,  nFeature * nFeatureTmpl);
		iFeature     = index[0] / nFeatureTmpl;
		iFeatureTmpl = index[0] % nFeatureTmpl;
		printf("Corelation = %f,  template No. = %2d,  %d  ", record[0], iFeatureTmpl, nFeatureTmpl);
		features[iFeature]->prFeatureInf();
		printf("Output %s and end of matching \n", dFileName.c_str());
	}
	now = clock(); printf("End of calculation %lf (sec)\n", (double)(now - start)/CLOCKS_PER_SEC);
	return 0;
}

FastGMatch::FastGMatch(int nx, int ny, int mkTemplate) {
	this->nx = nx; this->ny = ny;
	this->mkTemplate = mkTemplate;
	nxy = nx * ny;

	nAngleCoef = 2 * PANGLE + 1;

	/* Image data region */
	inImgI         = (int *) malloc(sizeof(int) * nxy);
	inImgOrg       = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);
	inImg          = (FPTYPE *) malloc(sizeof(FPTYPE) * nxy);


	int nEstFeatureTmpl = 0, nEstFeature = 0;
	scale = TEMPLSCALEINIT;
	while(scale <= TEMPLSCALEMAX) {
		nEstFeatureTmpl += NTEMPLATEROT * NTEMPLATEENLONG + 1;
		scale = (int) scale * TEMPLSCALERATIO;
	}
	featuresTmpl = (Feature **) malloc(sizeof(Feature *) * nEstFeatureTmpl);

	scale = SCALEINIT;
	while(scale <= SCALEMAX) {
		nEstFeature += (int) (nxy / (scale * STEPRATIO * scale * STEPRATIO) + ROUNDFRAC);
		scale = (int) scale * SCALERATIO;
	}

	nEstFeature *= 2;
	features     = (Feature **) malloc(sizeof(Feature *) * nEstFeature);
	extFeature  = new ExtFeature(nx, ny);
	imageIO     = new ImageIO(nx, ny);
	utilities   = new Utilities();
	extFeature->imageIO = imageIO; /* For debug */

	/* Make feature of template with various transformation */
	if (mkTemplate) gptTbl = new GptTbl();

}

/* Save edge histogram as images*/
int FastGMatch::histSave(FPTYPE *dirHist) {
	int nxy = nx * ny;
	std::string fileOutPre(WORKBASE), fileName;
	fileOutPre += IMGPRE; fileOutPre += BLURPRE; fileOutPre += OUTPRE;
	fileName = fileOutPre + "histc0.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), dirHist);
	fileName = fileOutPre + "histc1.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy]));
	fileName = fileOutPre + "histc2.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy * 3]));
	fileName = fileOutPre + "histc3.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy * 5]));
	fileName = fileOutPre + "histc4.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy * 7]));
	fileName = fileOutPre + "hists1.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy * 2]));
	fileName = fileOutPre + "hists2.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy * 4]));
	fileName = fileOutPre + "hists3.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy * 6]));
	fileName = fileOutPre + "hists4.pgm";
	imageIO->saveImageFileFp(fileName.c_str(), &(dirHist[nxy * 8]));
	return 0;
}

