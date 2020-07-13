#include <stdint.h>
#include <stdio.h>

#include "flow_test.h"

int main()
{
	int IMAGE_WIDTH = 40, IMAGE_HEIGHT = 40;
	int imageSize = IMAGE_WIDTH * IMAGE_HEIGHT;
	uint8_t image1[50][50], image2[50][50];

	for (int i = 0; i < IMAGE_HEIGHT; i++) {
		for (int j = 0; j < IMAGE_WIDTH; j++) {
			if (i > 15 && i < 25 && j > 15 && j < 25)
				image1[i][j] = 100;
			else
				image1[i][j] = 10;
		}
	}
	for (int i = 0; i < IMAGE_HEIGHT; i++) {
		for (int j = 0; j < IMAGE_WIDTH; j++) {
			if (i > 22 && i < 32 && j > 22 && j < 32)
				image2[i][j] = 100;
			else
				image2[i][j] = 10;
		}
	}
	//for (int i = 0; i < IMAGE_HEIGHT; i++) {
	//	for (int j = 0; j < IMAGE_WIDTH; j++) {
	//		printf("%d ", image1[i][j]);
	//	}
	//	printf("\n");
	//}
	//printf("\n");
	
	uint8_t *img1 = (uint8_t *)malloc(sizeof(uint8_t) * 40*40);
	uint8_t *img2 = (uint8_t *)malloc(sizeof(uint8_t) * 40*40);
	img1 = &image1;
	img2 = &image2;
	for (int i = 0; i < IMAGE_HEIGHT; i++) {
		for (int j = 0; j < IMAGE_WIDTH; j++) {
			*(img1 + i * 40 + j) = image1[i][j];
			*(img2 + i * 40 + j) = image2[i][j];
		}
	}


	//for (int i = 0; i < IMAGE_HEIGHT; i++) {
	//	for (int j = 0; j < IMAGE_WIDTH; j++) {
	//		printf("%d ",*(img1 + i*40 +j));
	//	}
	//	printf("\n");
	//}
	//printf("\n\n");
	//for (int i = 0; i < IMAGE_HEIGHT; i++) {
	//	for (int j = 0; j < IMAGE_WIDTH; j++) {
	//		printf("%d ", *(img2 + i * 40 + j));
	//	}
	//	printf("\n");
	//}

	float dx, dy;
	uint8_t qual = computeFlow(img1, img2, &dx, &dy);

	printf("qual: \ndx: %d, dy: %d \n", qual, dx, dy);

	system("pause");
	return 0;
}