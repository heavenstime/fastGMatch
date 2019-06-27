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

class GptPara {
public:
	GptPara();
	FPTYPE *gpt;
	void gptTransformImage(FPTYPE *inImg, FPTYPE *outImg, int nx, int ny, int cx, int cy);
	void gptPr(char *st);
	private:
	void gptInverse(GptPara *iGpt);
	void gptTransformPoint(FPTYPE *inP, FPTYPE *outP);
}; /* Work max for theta*/

class GptTbl {
public:
	GptTbl();
	int       nTran;
	GptPara **gptPara;
}; /* Work max for theta*/

class ImageIO {
	int nx, ny;
	char *iffX, *diffY;
public:
	ImageIO(int nx, int ny);
	int  rotImg(FPTYPE *inImg, FPTYPE *outImg, FPTYPE theta, int nxOut, int nyOut);
	void loadImageFile(const char *fileName, int *img);
	void saveImageFile(const char *fileName, int *img);
	void saveImageFileFp(const char *fileName, FPTYPE *img);
};

class Utilities {
public:
	Utilities();
	void heapSort(FPTYPE *record, int *index, int nData);
	void saveData(const char *fileName, FPTYPE *data, int nData);
	private:
	void mkHeap(FPTYPE *record, int *index, int leaf, int root);
};

#endif /* UTILITIES_H_ */
