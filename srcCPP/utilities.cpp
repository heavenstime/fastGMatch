/*
 * utilities.cpp
 *
 *  Created on: 2019/03/12
 *      Author: yamasita
 */

#include <cmath>
#include <string>

#include "parameter.h"
#include "utilities.h"
#include "extFeature.h"
#include "fastGMatch.h"

Utilities::Utilities() {

}
/* Heap sort */
void Utilities::heapSort(FPTYPE *record, int *index, int nData) {
	int    root, leaf, tmpIndex;
	FPTYPE tmp;
	leaf = nData - 1;      /* Leaf (position is the index of array) */
	root = nData / 2 - 1;  /* Root*/

	/* Initial semi-ordered tree construction */
	while (root >= 0) {
		mkHeap(record, index, leaf, root);
		--root;
	}
	/* Extract the smallest element from semi-ordered tree */
	while (leaf > 0) {
		/* Exchange root and leaf */
		tmp           = record [0];
		record [0]    = record [leaf];
		record [leaf] = tmp ;
		tmpIndex      = index[0];
		index[0]      = index[leaf];
		index[leaf]   = tmpIndex;
		/* Reconstruct semi-ordered tree */
		--leaf;
		mkHeap(record , index, leaf, 0);
	}
}

/* Make semi-ordered tree when root is not correct position for heap sort*/
void Utilities::mkHeap(FPTYPE *record, int *index, int leaf, int root) {
	int min, tmpIndex;
	int child = root * 2 + 1;
	FPTYPE tmp;

	while(child <= leaf) {
		if (child < leaf && record[child + 1] < record[child]) {
			min = child + 1;
		} else {
			min = child;
		}
		if (record[root] <= record[min]) break;
		tmp          = record[root];
		record[root] = record[min];
		record[min]  = tmp;
		tmpIndex     = index[root];
		index[root]  = index[min];
		index[min]   = tmpIndex;
		root          = min;
		child         = root * 2 + 1;
	}
}

/* save FP data as text */
void  Utilities::saveData(const char *fileName, FPTYPE *data, int nData){
	FILE *fp; /* File pointer */
	int  pos; /* Loop variable */

	fp = fopen(fileName, "w");
	if (fp == NULL) {
		printf("     Cannot open! %s \n\n", fileName);
		exit(1);
	}
	/* Output of image data */
	for (pos = 0 ; pos < nData ; ++pos) {
		fprintf(fp, "%d  %e\n", pos, data[pos]);
	}
	fclose(fp);
}

ImageIO::ImageIO(int nx, int ny) {
	this->nx = nx;
	this->ny = ny;
}


/* Rotate image (This can be done by gptTransformImage() )*/
int ImageIO::rotImg(FPTYPE *inImg, FPTYPE *outImg, FPTYPE theta, int nxOut, int nyOut) {
	int cx = nx / 2, cy = ny / 2;
	int cxOut = nxOut / 2, cyOut = nyOut / 2;
	int ixOut, iyOut;
	int itx, ity, pos;
	FPTYPE tx, ty, fx, fy;
	FPTYPE cost = cos(theta), sint = sin(theta);


	for (iyOut = 0 ; iyOut < nyOut ; ++iyOut) {
		for (ixOut = 0 ; ixOut < nxOut ; ++ixOut) {
			tx =   cost * (ixOut - cxOut) + sint * (iyOut - cyOut) + cx;
			ty = - sint * (ixOut - cxOut) + cost * (iyOut - cyOut) + cy;
			itx = (int) tx;
			ity = (int) ty;
			fx  = tx - itx; fy = ty - ity;
			pos = itx + nx * ity;
			outImg[ixOut + nxOut * iyOut] = (1.0 - fx) * (1.0 - fy) * inImg[pos] + fx * (1.0 - fy) * inImg[pos + 1]
			                              + (1.0 - fx) * fy * inImg[pos + nx] + fx * fy * inImg[pos + nx + 1];
		}
	}
	return 0;
}

/* Load image (Input of header & body information of pgm file) */
void ImageIO::loadImageFile(const char *fileName, int *img){
	/* unsigned char buffer[MAX_BUFFERSIZE];*/
	char buffer[MAX_BUFFERSIZE];
	FILE *fp;         /* File pointer */
	int ntx, nty;     /* Original image size */
	int max_gray;    /* Maximum gray level */
	int pos, ix, iy; /* Loop variable */

	/* Input file open */
	//printf("Input file = %s\n", fileName);
	fp = fopen(fileName, "rb");
	if (fp == NULL) {
		printf("     The file doesn't exist! %s \n\n", fileName);
		exit(1);
	}

	/* Check of file-type ---P5 */
	if (fgets(buffer, MAX_BUFFERSIZE, fp) == 0) {
		printf("     Read error of image file %s \n\n", fileName);
		exit(1);
	}
	if (buffer[0] != 'P' || buffer[1] != '5') {
		printf("     Mistaken file format, not P5! %s \n\n", fileName);
		exit(1);
	}

	/* input of x_size1, y_size1 */
	ntx = 0; nty = 0;
	while (ntx == 0 || nty == 0) {
		if (fgets(buffer, MAX_BUFFERSIZE, fp) == 0) {
			printf("     Read error of image file %s \n\n", fileName);
			exit(1);
		}
		if (buffer[0] != '#')  sscanf(buffer, "%d %d", & ntx, & nty);
	}

	if (ntx != nx || nty != ny) {
		printf("     Image size is not correct ! %s  (%d, %d) != (%d, %d) \n\n", fileName, ntx, nty, nx, ny);
		exit(1);
	}

	/* input of max_gray */
	max_gray = 0;
	while (max_gray == 0) {
		if (fgets(buffer, MAX_BUFFERSIZE, fp) == 0) {
			printf("     Read error of image file %s \n\n", fileName);
			exit(1);
		}
		if (buffer[0] != '#')  sscanf(buffer, "%d", &max_gray);
	}

	if (max_gray != WHITE) {
		printf("     Invalid value of maximum gray level! %s \n\n", fileName);
		exit(1);
	}

	/* Input of image data */
	pos = 0;
	for (iy = 0 ; iy < nty ; ++iy) {
		for (ix = 0 ; ix < ntx ; ++ix) {
			img[pos++] = (int) ((unsigned char) fgetc(fp));
		}
	}
	fclose(fp);
//#define SIMPLEIMAGE
#ifdef SIMPLEIMAGE
	for (iy = 0 ; iy < ny ; ++iy) {
		for (ix = 0 ; ix < nx ; ++ix) {
//		if (ix % 64 <32)  	img[ix + 2 * nx * iy] = 1; else img[ix + 2 * nx * iy] = 0;
//			if ((iy - 255) > 0 && (ix - 255) > 0) 	img[ix + nx * iy] = 1.0; else img[ix + 2 * nx * iy] = 0;/* Small square */
			 if (ix * 2 >= ny)  	img[ix + nx * iy] = 1; else img[ix + nx * iy] = 0;
		}
	}
#endif
}

/* Save image (FP). Values are normalized */
void ImageIO::saveImageFileFp(const char *fileName, FPTYPE *img){
	FILE *fp; /* File pointer */
	int  pos; /* Loop variable */
	FPTYPE maxVal, minVal, coef;

	fp = fopen(fileName, "wb");
	if (fp == NULL) {
		printf("     Cannot open! %s \n\n", fileName);
		exit(1);
	}
	/* output of pgm file header information */
	fputs("P5\n", fp);
	fputs("# Created by Image Processing\n", fp);
	fprintf(fp, "%d %d\n", nx, ny);
	fprintf(fp, "%d\n", WHITE);

	/* Output of image data */
	maxVal = minVal = img[0] ;
	for (pos = 1 ; pos < nx * ny ; ++pos) {
		if(maxVal < img[pos]) maxVal = img[pos];
		else if(minVal > img[pos]) minVal = img[pos];
	}
	coef = (WHITE - BLACK) / (maxVal - minVal);
	// printf("%d  %d   %f  %f \n", nx, ny, maxVal, minVal);
	for (pos = 0 ; pos < nx * ny ; ++pos) {
		fputc((int) (coef * (img[pos] - minVal)), fp);
	}
	fclose(fp);
}

/* Save image (int) */
void ImageIO::saveImageFile(const char *fileName, int *img){
	FILE *fp; /* File pointer */
	int  pos; /* Loop variable */

	fp = fopen(fileName, "wb");
	if (fp == NULL) {
		printf("     Cannot open! %s \n\n", fileName);
		exit(1);
	}
	/* output of pgm file header information */
	fputs("P5\n", fp);
	fputs("# Created by Image Processing\n", fp);
	fprintf(fp, "%d %d\n", nx, ny);
	fprintf(fp, "%d\n", WHITE);

	/* Output of image data */
	for (pos = 0 ; pos < nx * ny ; ++pos) {
		if(img[pos] > WHITE) fputc(WHITE, fp);
		else if(img[pos] < BLACK) fputc(BLACK, fp);
		else fputc(img[pos], fp);
	}
	fclose(fp);
}

GptTbl::GptTbl() {
//#define ONLYROT
#ifdef ONLYROT /* Only rotatation */
		FPTYPE cosTran, sinTran;
		nTran    = NTEMPLATEROT;
		gptPara  = (GptPara **) malloc(sizeof(GptPara *) * nTran);
		int pos = 0;
		for (int iTheta = 0 ; iTheta < NTEMPLATEROT ; ++iTheta) {
			gptPara[iTheta] = new GptPara();
			gptPara[iTheta]->gpt = (FPTYPE *) malloc(sizeof(FPTYPE *) * 9);
			cosTran = cos((2.0 * PI / NTEMPLATEROT) * iTheta); 		sinTran = sin((2.0 * PI / NTEMPLATEROT) * iTheta);
			gptPara[pos]->gpt[C11] = cosTran;
			gptPara[pos]->gpt[C22] = cosTran;
			gptPara[pos]->gpt[C12] = - sinTran; /* Note that upside down */
			gptPara[pos]->gpt[C21] = sinTran;
			gptPara[pos]->gpt[C13] = gptPara[pos]->gpt[C23] = gptPara[pos]->gpt[C31]  = gptPara[pos]->gpt[C32] = 0.0;
			gptPara[pos]->gpt[C33] = 1.0;
			//gptPr(&(gptTbl->gptTbl[pos]), "GPTInit ");
			++pos;
		}
#else
	FPTYPE cosTran, sinTran, alpha, beta;
	nTran   = NTEMPLATEROT * NTEMPLATEENLONG + 1;
	gptPara  = (GptPara **) malloc(sizeof(GptPara *) * nTran);
	for (int iTran = 0 ; iTran < nTran ; ++iTran) 	gptPara[iTran] = new GptPara();

	gptPara[0]->gpt[C11] = gptPara[0]->gpt[C22] = 1.0;
	gptPara[0]->gpt[C12] = gptPara[0]->gpt[C21] = 0.0;
	gptPara[0]->gpt[C13] = gptPara[0]->gpt[C23] = gptPara[0]->gpt[C31]  = gptPara[0]->gpt[C32] = 0.0;
	gptPara[0]->gpt[C33] = 1.0;
	int pos = 0;
	for (int iTheta = 0 ; iTheta < NTEMPLATEROT ; ++iTheta) {
		cosTran = cos((PI / NTEMPLATEROT) * iTheta); 		sinTran = sin((PI / NTEMPLATEROT) * iTheta);
		for (int iEnlong = 1 ; iEnlong <= NTEMPLATEENLONG ; ++iEnlong) {
			++pos;
			alpha = 1.0 + TEMPLATEENCOEF * iEnlong; beta = 1.0 - TEMPLATEENCOEF * iEnlong;
			//printf("alpha %f  beta %f\n", alpha, beta);
			gptPara[pos]->gpt[C11] = alpha * cosTran * cosTran + beta  * sinTran * sinTran;
			gptPara[pos]->gpt[C22] = beta  * cosTran * cosTran + alpha * sinTran * sinTran;
			gptPara[pos]->gpt[C12] = gptPara[pos]->gpt[C21] = (alpha - beta) * cosTran * sinTran;
			gptPara[pos]->gpt[C13] = gptPara[pos]->gpt[C23] = gptPara[pos]->gpt[C31]  = gptPara[pos]->gpt[C32] = 0.0;
			gptPara[pos]->gpt[C33] = 1.0;
			//gptPr(&(gptTbl[pos]), "GPTInit ");
		}
	}
#endif
}

GptPara::GptPara() {
	gpt = (FPTYPE *) malloc(sizeof(FPTYPE *) * 9);
}
/* Inverse projection transformation */
void GptPara::gptInverse(GptPara *iGpt) {
	FPTYPE det;
	det =   gpt[C11] * gpt[C22] * gpt[C33]
	      + gpt[C21] * gpt[C32] * gpt[C13]
	      + gpt[C31] * gpt[C12] * gpt[C23]
	      - gpt[C11] * gpt[C32] * gpt[C23]
	      - gpt[C21] * gpt[C12] * gpt[C33]
	      - gpt[C31] * gpt[C22] * gpt[C13];
	if (fabs(det) < EPS) {
		printf("Singular GPT is generated!!! in gptInverse() \n");
	}
	iGpt->gpt[C11] =  (gpt[C22] * gpt[C33] - gpt[C32] * gpt[C23]) / det;
	iGpt->gpt[C21] = -(gpt[C21] * gpt[C33] - gpt[C31] * gpt[C23]) / det;
	iGpt->gpt[C31] =  (gpt[C21] * gpt[C32] - gpt[C31] * gpt[C22]) / det;

	iGpt->gpt[C12] = -(gpt[C12] * gpt[C33] - gpt[C32] * gpt[C13]) / det;
	iGpt->gpt[C22] =  (gpt[C11] * gpt[C33] - gpt[C31] * gpt[C13]) / det;
	iGpt->gpt[C32] = -(gpt[C11] * gpt[C32] - gpt[C31] * gpt[C12]) / det;

	iGpt->gpt[C13] =  (gpt[C12] * gpt[C23] - gpt[C22] * gpt[C13]) / det;
	iGpt->gpt[C23] = -(gpt[C11] * gpt[C23] - gpt[C21] * gpt[C13]) / det;
	iGpt->gpt[C33] =  (gpt[C11] * gpt[C22] - gpt[C21] * gpt[C12]) / det;
}

void GptPara::gptTransformPoint(FPTYPE *inP, FPTYPE *outP) {
	int i, j;
	FPTYPE sum;
	for(i = 0 ; i < 3 ; ++i) {
		sum = 0.0;
		for(j = 0 ; j < 3 ; ++j) 	sum += gpt[i + j * 3] * inP[j];
		outP[i] = sum;
	}
}

/* Projection transformation of the image by bilinear interpolation */
void GptPara::gptTransformImage(FPTYPE *inImg, FPTYPE *outImg, int nx, int ny, int cx, int cy) {
	int x, y, xOrg, yOrg, posOrg;
	FPTYPE xOrgFl, xOrgFrac, yOrgFl, yOrgFrac;
	FPTYPE inVect[3], outVect[3];
	GptPara invGpt;

	/* Inverse transform of GPT */
	gptInverse(&invGpt);

	/* Output image generation by bilinear interpolation */
	inVect[2] = 1.0;
	for (y = 0 ; y < ny ; y++) {
		inVect[1] = y - cy;
		for (x = 0 ; x < nx ; x++) {
			inVect[0] = x - cx;
			invGpt.gptTransformPoint(inVect, outVect);
			xOrgFl = outVect[0] / outVect[2] + cx;
			yOrgFl = outVect[1] / outVect[2] + cy;
			xOrg   = (int) floor(xOrgFl);
			yOrg   = (int) floor(yOrgFl);
			xOrgFrac = xOrgFl - xOrg;
			yOrgFrac = yOrgFl - yOrg;
			if (xOrg >= 0 && xOrg + 1 < nx && yOrg >= 0 && yOrg + 1 < ny) {
				posOrg = yOrg * nx + xOrg;
				outImg[y * nx + x] = (
						(1.0 - yOrgFrac) * ((1.0 - xOrgFrac) * inImg[posOrg] + xOrgFrac * inImg[posOrg + 1])
						+ yOrgFrac * ((1.0 - xOrgFrac) * inImg[posOrg + nx] + xOrgFrac * inImg[posOrg + nx + 1]));
			} else {
				outImg[y * nx + x] = BLANK;
			}
		}
	}
}

/* Show GPT parameter */
void GptPara::gptPr(char *st) {
	int i;
	printf("%s \n", st);
	for(i = 0 ; i < 3 ; ++i) {
		printf("%10.6f  %10.6f  %10.6f\n", gpt[i], gpt[i + 3], gpt[i + 6]);
	}
}
