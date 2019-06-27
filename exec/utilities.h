/*
 * utilities.h
 *
 *  Created on: 2018/05/01
 *      Author: yamasita
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_


#define B11 0
#define B21 1
#define B12 2
#define B22 3

#define C11 0
#define C21 1
#define C31 2
#define C12 3
#define C22 4
#define C32 5
#define C13 6
#define C23 7
#define C33 8
#define C14 9
#define C24 10
#define C34 11

void gptTransformImage(FPTYPE *inGpt, FPTYPE *inImg, FPTYPE *outImg, int nx, int ny, int cx, int cy);
void heapSort(FPTYPE *record, int *index, int nData);
void saveFeatures(char *fileName, Features *features);
void loadFeatures(char *fileName, Features *features);
void prFeatureInf(Feature *feature);
void gptPr(FPTYPE *gpt, char *st);
int rotImg(FPTYPE *inImg, FPTYPE *outImg, FPTYPE theta, int nx, int ny, int nxOut, int nyOut);
void saveImageFile(char *fileName, int *img, int nx, int ny);
void saveImageFileFp(char *fileName, FPTYPE *img, int nx, int ny);
void loadImageFile(char *fileName, int *img, int nx, int ny);
void saveData(char *fileName, FPTYPE *data, int nData);

int histSave(FPTYPE *dirHist, int nx, int ny);


#endif /* UTILITIES_H_ */
