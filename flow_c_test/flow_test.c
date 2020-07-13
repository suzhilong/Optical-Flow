#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>


#define IMAGE_WIDTH	40//853
#define IMAGE_HEIGHT 40//640
#define SEARCH_SIZE	8
#define TILE_SIZE	8
#define NUM_BLOCKS	5
#define DIFF_THRESHOLD 40
#define SAD_THRESHOLD 5000
#define HIST_FILTER true

#define HISTSIZE 2 * (SEARCH_SIZE+SEARCH_SIZE+1) + 1


uint32_t __USAD8(uint32_t val1, uint32_t val2)
{
	uint32_t res = 0;
	uint8_t *val1_bytes = (uint8_t *)(&val1);
	uint8_t *val2_bytes = (uint8_t *)(&val2);

	for (int i = 0; i < 4; i++) {
		int16_t v1 = val1_bytes[i];
		int16_t v2 = val2_bytes[i];
		res += (uint32_t)(abs(v1 - v2));
	}

	return res;
}

uint32_t __USADA8(uint32_t val1, uint32_t val2, uint32_t val3)
{
	uint32_t res = val3;
	uint8_t *val1_bytes = (uint8_t *)(&val1);
	uint8_t *val2_bytes = (uint8_t *)(&val2);

	for (int i = 0; i < 4; i++) {
		int16_t v1 = val1_bytes[i];
		int16_t v2 = val2_bytes[i];
		res += (uint32_t)(abs(v1 - v2));
	}

	return res;
}

uint32_t __UHADD8(uint32_t val1, uint32_t val2)
{
	uint32_t res = 0;
	uint8_t *res_bytes = (uint8_t *)(&res);
	uint8_t *val1_bytes = (uint8_t *)(&val1);
	uint8_t *val2_bytes = (uint8_t *)(&val2);

	for (int i = 0; i < 4; i++) {
		res_bytes[i] = (val1_bytes[i] + val2_bytes[i]) >> 1;
	}

	return res;
}

/**
* @brief Compute the average pixel gradient of all horizontal and vertical steps
*
* TODO compute_diff is not appropriate for low-light mode images
*
* @param image  the array holding pixel data
* @param offX   x coordinate of upper left corner of 8x8 pattern in image
* @param offY   y coordinate of upper left corner of 8x8 pattern in image
*
* @return       value indicating how suitable a pattern is for optical flow calculation
*/
static inline uint32_t computeDiff(uint8_t *image, uint16_t offX, uint16_t offY, uint16_t row_size)
{//计算是否为特征点
	/* calculate position in image buffer */
	uint16_t off = (offY + 2) * row_size + (offX + 2);
	/* we calculate only the center 4x4 pattern */
	uint32_t acc;

	/* calculate vertical gradient */
	acc = __USAD8(*((uint32_t*)&image[off + 0 + 0 * row_size]), *((uint32_t*)&image[off + 0 + 1 * row_size]));
	acc = __USADA8(*((uint32_t*)&image[off + 0 + 1 * row_size]), *((uint32_t*)&image[off + 0 + 2 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image[off + 0 + 2 * row_size]), *((uint32_t*)&image[off + 0 + 3 * row_size]), acc);

	/* we need to get columns */
	uint32_t col1 = (image[off + 0 + 0 * row_size] << 24) | image[off + 0 + 1 * row_size] << 16 | image[off + 0 + 2 * row_size] << 8 | image[off + 0 + 3 * row_size];
	uint32_t col2 = (image[off + 1 + 0 * row_size] << 24) | image[off + 1 + 1 * row_size] << 16 | image[off + 1 + 2 * row_size] << 8 | image[off + 1 + 3 * row_size];
	uint32_t col3 = (image[off + 2 + 0 * row_size] << 24) | image[off + 2 + 1 * row_size] << 16 | image[off + 2 + 2 * row_size] << 8 | image[off + 2 + 3 * row_size];
	uint32_t col4 = (image[off + 3 + 0 * row_size] << 24) | image[off + 3 + 1 * row_size] << 16 | image[off + 3 + 2 * row_size] << 8 | image[off + 3 + 3 * row_size];

	/* calculate horizontal gradient */
	acc = __USADA8(col1, col2, acc);
	acc = __USADA8(col2, col3, acc);
	acc = __USADA8(col3, col4, acc);

	return acc;

}


static inline uint32_t computeSubpixel(uint8_t *image1, uint8_t *image2, uint16_t off1X, uint16_t off1Y, uint16_t off2X, uint16_t off2Y, uint32_t *acc, uint16_t row_size)
{
	/* calculate position in image buffer */
	uint16_t off1 = off1Y * row_size + off1X; // image1
	uint16_t off2 = off2Y * row_size + off2X; // image2

	uint32_t s0, s1, s2, s3, s4, s5, s6, s7, t1, t3, t5, t7;

	for (uint16_t i = 0; i < 8; i++) {
		acc[i] = 0;
	}


	/*
	* calculate for each pixel in the 8x8 field with upper left corner (off1X / off1Y)
	* every iteration is one line of the 8x8 field.
	*
	*  + - + - + - + - + - + - + - + - +
	*  |   |   |   |   |   |   |   |   |
	*  + - + - + - + - + - + - + - + - +
	*/

	for (uint16_t i = 0; i < 8; i++) {
		/*
		* first column of 4 pixels:
		*
		*  + - + - + - + - + - + - + - + - +
		*  | x | x | x | x |   |   |   |   |
		*  + - + - + - + - + - + - + - + - +
		*
		* the 8 s values are from following positions for each pixel (X):
		*  + - + - + - +
		*  +   5   7   +
		*  + - + 6 + - +
		*  +   4 X 0   +
		*  + - + 2 + - +
		*  +   3   1   +
		*  + - + - + - +
		*
		*  variables (s1, ...) contains all 4 results (32bit -> 4 * 8bit values)
		*/

		/* compute average of two pixel values */
		s0 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 + 1 + (i + 0) * row_size])));
		s1 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i + 1) * row_size]), *((uint32_t*)&image2[off2 + 1 + (i + 1) * row_size])));
		s2 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 + 0 + (i + 1) * row_size])));
		s3 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i + 1) * row_size]), *((uint32_t*)&image2[off2 - 1 + (i + 1) * row_size])));
		s4 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 - 1 + (i + 0) * row_size])));
		s5 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i - 1) * row_size]), *((uint32_t*)&image2[off2 - 1 + (i - 1) * row_size])));
		s6 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 + 0 + (i - 1) * row_size])));
		s7 = (__UHADD8(*((uint32_t*)&image2[off2 + 0 + (i - 1) * row_size]), *((uint32_t*)&image2[off2 + 1 + (i - 1) * row_size])));

		/* these 4 t values are from the corners around the center pixel */
		t1 = (__UHADD8(s0, s1));
		t3 = (__UHADD8(s3, s4));
		t5 = (__UHADD8(s4, s5));
		t7 = (__UHADD8(s7, s0));

		/*
		* finally we got all 8 subpixels (s0, t1, s2, t3, s4, t5, s6, t7):
		*  + - + - + - +
		*  |   |   |   |
		*  + - 5 6 7 - +
		*  |   4 X 0   |
		*  + - 3 2 1 - +
		*  |   |   |   |
		*  + - + - + - +
		*/


		/* fill accumulation vector */
		acc[0] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), s0, acc[0]);
		acc[1] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), t1, acc[1]);
		acc[2] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), s2, acc[2]);
		acc[3] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), t3, acc[3]);
		acc[4] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), s4, acc[4]);
		acc[5] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), t5, acc[5]);
		acc[6] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), s6, acc[6]);
		acc[7] = __USADA8((*((uint32_t*)&image1[off1 + 0 + i * row_size])), t7, acc[7]);

		/*
		* same for second column of 4 pixels:
		*
		*  + - + - + - + - + - + - + - + - +
		*  |   |   |   |   | x | x | x | x |
		*  + - + - + - + - + - + - + - + - +
		*/

		s0 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 + 5 + (i + 0) * row_size])));
		s1 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i + 1) * row_size]), *((uint32_t*)&image2[off2 + 5 + (i + 1) * row_size])));
		s2 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 + 4 + (i + 1) * row_size])));
		s3 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i + 1) * row_size]), *((uint32_t*)&image2[off2 + 3 + (i + 1) * row_size])));
		s4 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 + 3 + (i + 0) * row_size])));
		s5 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i - 1) * row_size]), *((uint32_t*)&image2[off2 + 3 + (i - 1) * row_size])));
		s6 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i + 0) * row_size]), *((uint32_t*)&image2[off2 + 4 + (i - 1) * row_size])));
		s7 = (__UHADD8(*((uint32_t*)&image2[off2 + 4 + (i - 1) * row_size]), *((uint32_t*)&image2[off2 + 5 + (i - 1) * row_size])));

		t1 = (__UHADD8(s0, s1));
		t3 = (__UHADD8(s3, s4));
		t5 = (__UHADD8(s4, s5));
		t7 = (__UHADD8(s7, s0));

		acc[0] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), s0, acc[0]);
		acc[1] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), t1, acc[1]);
		acc[2] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), s2, acc[2]);
		acc[3] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), t3, acc[3]);
		acc[4] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), s4, acc[4]);
		acc[5] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), t5, acc[5]);
		acc[6] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), s6, acc[6]);
		acc[7] = __USADA8((*((uint32_t*)&image1[off1 + 4 + i * row_size])), t7, acc[7]);
	}

	return 0;
}

static inline uint32_t computeSAD(uint8_t *image1, uint8_t *image2, uint16_t off1X, uint16_t off1Y, uint16_t off2X, uint16_t off2Y, uint16_t row_size)
{
	/* calculate position in image buffer */
	uint16_t off1 = off1Y * row_size + off1X; // image1
	uint16_t off2 = off2Y * row_size + off2X; // image2

	uint32_t acc;
	acc = __USAD8(*((uint32_t*)&image1[off1 + 0 + 0 * row_size]), *((uint32_t*)&image2[off2 + 0 + 0 * row_size]));
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 0 * row_size]), *((uint32_t*)&image2[off2 + 4 + 0 * row_size]), acc);

	acc = __USADA8(*((uint32_t*)&image1[off1 + 0 + 1 * row_size]), *((uint32_t*)&image2[off2 + 0 + 1 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 1 * row_size]), *((uint32_t*)&image2[off2 + 4 + 1 * row_size]), acc);

	acc = __USADA8(*((uint32_t*)&image1[off1 + 0 + 2 * row_size]), *((uint32_t*)&image2[off2 + 0 + 2 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 2 * row_size]), *((uint32_t*)&image2[off2 + 4 + 2 * row_size]), acc);

	acc = __USADA8(*((uint32_t*)&image1[off1 + 0 + 3 * row_size]), *((uint32_t*)&image2[off2 + 0 + 3 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 3 * row_size]), *((uint32_t*)&image2[off2 + 4 + 3 * row_size]), acc);

	acc = __USADA8(*((uint32_t*)&image1[off1 + 0 + 4 * row_size]), *((uint32_t*)&image2[off2 + 0 + 4 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 4 * row_size]), *((uint32_t*)&image2[off2 + 4 + 4 * row_size]), acc);

	acc = __USADA8(*((uint32_t*)&image1[off1 + 0 + 5 * row_size]), *((uint32_t*)&image2[off2 + 0 + 5 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 5 * row_size]), *((uint32_t*)&image2[off2 + 4 + 5 * row_size]), acc);

	acc = __USADA8(*((uint32_t*)&image1[off1 + 0 + 6 * row_size]), *((uint32_t*)&image2[off2 + 0 + 6 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 6 * row_size]), *((uint32_t*)&image2[off2 + 4 + 6 * row_size]), acc);

	acc = __USADA8(*((uint32_t*)&image1[off1 + 0 + 7 * row_size]), *((uint32_t*)&image2[off2 + 0 + 7 * row_size]), acc);
	acc = __USADA8(*((uint32_t*)&image1[off1 + 4 + 7 * row_size]), *((uint32_t*)&image2[off2 + 4 + 7 * row_size]), acc);

	return acc;
}


uint8_t computeFlow(uint8_t *image1, uint8_t *image2, float *pixel_flow_x, float *pixel_flow_y)
{
	/* constants */
	const int16_t winmin = -SEARCH_SIZE;
	const int16_t winmax = SEARCH_SIZE;
	const uint16_t hist_size = 2 * (winmax - winmin + 1) + 1;//2 * (8+8) + 1

	/* variables */
	uint16_t pixLo = SEARCH_SIZE + 1;//图片最左边一块中间的像素
	uint16_t pixHi = IMAGE_WIDTH - (SEARCH_SIZE + 1) - TILE_SIZE;//图片最右边一块中间的像素 //IMAGE_WIDTH是图片宽度
	uint16_t pixStep = (pixHi - pixLo) / NUM_BLOCKS + 1;
	uint16_t i, j;
	uint32_t acc[8]; // subpixels
	//uint16_t histX[hist_size]; // counter for x shift
	//uint16_t histY[hist_size]; // counter for y shift
	uint16_t histX[HISTSIZE];
	uint16_t histY[35];
	int8_t dirsX[64]; // shift directions in x
	int8_t dirsy[64]; // shift directions in y
	uint8_t subdirs[64]; // shift directions of best subpixels
	float meanflowX = 0.0f;
	float meanflowY = 0.0f;
	uint16_t meancount = 0;
	float histflowX = 0.0f;
	float histflowY = 0.0f;

	/* initialize with 0 */
	for (j = 0; j < hist_size; j++) { histX[j] = 0; histY[j] = 0; }

	/* iterate over all patterns*/
	for (j = pixLo; j < pixHi; j += pixStep) {
		for (i = pixLo; i < pixHi; i += pixStep) {
			/* test pixel if it is suitable for flow tracking */
			uint32_t diff = computeDiff(image1, i, j, (uint16_t)IMAGE_WIDTH);//计算是否为特征点
			if (diff < DIFF_THRESHOLD)//diff太小，不满足特征点的要求
			{
				continue;
			}
			uint32_t dist = 0xFFFFFFFF; // set initial distance to "infinity"
			int8_t sumX = 0;
			int8_t sumY = 0;
			int8_t ii, jj;
			uint8_t *base1 = image1 + j * (uint16_t)IMAGE_WIDTH + i;
			for (jj = winmin; jj <= winmax; jj++) {//得到winmin~winmax内最相似的偏移量
				uint8_t *base2 = image2 + (j + jj) * (uint16_t)IMAGE_WIDTH + i;
				for (ii = winmin; ii <= winmax; ii++) {
					 uint32_t temp_dist = computeSAD(image1, image2, i, j, i + ii, j + jj, (uint16_t) IMAGE_WIDTH);

					if (temp_dist < dist) {
						sumX = ii;//x偏移量
						sumY = jj;//y偏移量
						dist = temp_dist;
					}
				}
			}
			/* acceptance SAD distance threshold */
			if (dist < SAD_THRESHOLD) {
				meanflowX += (float)sumX;
				meanflowY += (float)sumY;
				//8个方向的半像素寻找最佳匹配，把八个方向的半像素SAD保存到acc
				computeSubpixel(image1, image2, i, j, i + sumX, j + sumY, acc, (uint16_t)IMAGE_WIDTH);
				uint32_t mindist = dist; // best SAD until now
				uint8_t mindir = 8; // direction 8 for no direction
				for (uint8_t k = 0; k < 8; k++) {
					if (acc[k] < mindist) {
						// SAD becomes better in direction k
						mindist = acc[k];
						mindir = k;
					}
				}
				dirsX[meancount] = sumX;
				dirsy[meancount] = sumY;
				subdirs[meancount] = mindir;
				meancount++;//记录达到SAD要求的特征点个数
							/* feed histogram filter*/
				uint8_t hist_index_x = 2 * sumX + (winmax - winmin + 1);
				if (mindir == 0 || mindir == 1 || mindir == 7) hist_index_x += 1;
				if (mindir == 3 || mindir == 4 || mindir == 5) hist_index_x += -1;
				uint8_t hist_index_y = 2 * sumY + (winmax - winmin + 1);
				if (mindir == 5 || mindir == 6 || mindir == 7) hist_index_y += -1;
				if (mindir == 1 || mindir == 2 || mindir == 3) hist_index_y += 1;
				// 4个方向的直方图
				histX[hist_index_x]++;
				histY[hist_index_y]++;
			}
		}
	}

/* create flow image if needed (image1 is not needed anymore)
	* -> can be used for debugging purpose
	*
	if (FLOAT_AS_BOOL(global_data.param[PARAM_USB_SEND_VIDEO]))//&& global_data.param[PARAM_VIDEO_USB_MODE] == FLOW_VIDEO)
	{
		for (j = pixLo; j < pixHi; j += pixStep) {
			for (i = pixLo; i < pixHi; i += pixStep) {
				uint32_t diff = compute_diff(image1, i, j, (uint16_t)IMAGE_WIDTH);
				if (diff > DIFF_THRESHOLD) {
					image1[j * ((uint16_t)IMAGE_WIDTH) + i] = 255;
				}
			}
		}
	}
*/

	/* evaluate flow calculation */
	if (meancount > 10) {
		meanflowX /= meancount;
		meanflowY /= meancount;
		int16_t maxpositionx = 0;
		int16_t maxpositiony = 0;
		uint16_t maxvaluex = 0;
		uint16_t maxvaluey = 0;
		/* position of maximal histogram peek */
		for (j = 0; j < hist_size; j++) {
			if (histX[j] > maxvaluex) {
				maxvaluex = histX[j];
				maxpositionx = j;
			}
			if (histY[j] > maxvaluey) {
				maxvaluey = histY[j];
				maxpositiony = j;
			}
		}
		/* check if there is a peak value in histogram */
		if (1) //(histX[maxpositionx] > meancount / 6 && histY[maxpositiony] > meancount / 6)
		{
			if (HIST_FILTER) {//滤波法
																				  /* use histogram filter peek value */
				uint16_t hist_x_min = maxpositionx;
				uint16_t hist_x_max = maxpositionx;
				uint16_t hist_y_min = maxpositiony;
				uint16_t hist_y_max = maxpositiony;
				/* x direction */
				if (maxpositionx > 1 && maxpositionx < hist_size - 2) {
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx + 2;
				} else if (maxpositionx == 0) {
					hist_x_min = maxpositionx;
					hist_x_max = maxpositionx + 2;
				} else if (maxpositionx == hist_size - 1) {
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx;
				} else if (maxpositionx == 1) {
					hist_x_min = maxpositionx - 1;
					hist_x_max = maxpositionx + 2;
				} else if (maxpositionx == hist_size - 2) {
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx + 1;
				}
				/* y direction */
				if (maxpositiony > 1 && maxpositiony < hist_size - 2) {
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony + 2;
				} else if (maxpositiony == 0) {
					hist_y_min = maxpositiony;
					hist_y_max = maxpositiony + 2;
				} else if (maxpositiony == hist_size - 1) {
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony;
				} else if (maxpositiony == 1) {
					hist_y_min = maxpositiony - 1;
					hist_y_max = maxpositiony + 2;
				} else if (maxpositiony == hist_size - 2) {
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony + 1;
				}
				float hist_x_value = 0.0f;
				float hist_x_weight = 0.0f;
				float hist_y_value = 0.0f;
				float hist_y_weight = 0.0f;
				for (uint8_t h = hist_x_min; h < hist_x_max + 1; h++) {
					hist_x_value += (float)(h*histX[h]);
					hist_x_weight += (float)histX[h];
				}
				for (uint8_t h = hist_y_min; h<hist_y_max + 1; h++) {
					hist_y_value += (float)(h*histY[h]);
					hist_y_weight += (float)histY[h];
				}
				histflowX = (hist_x_value / hist_x_weight - (winmax - winmin + 1)) / 2.0f;
				histflowY = (hist_y_value / hist_y_weight - (winmax - winmin + 1)) / 2.0f;
			} else {//平均法
					/* use average of accepted flow values */
				uint32_t meancount_x = 0;
				uint32_t meancount_y = 0;
				for (uint8_t h = 0; h < meancount; h++) {
					float subdirx = 0.0f;
					if (subdirs[h] == 0 || subdirs[h] == 1 || subdirs[h] == 7) subdirx = 0.5f;
					if (subdirs[h] == 3 || subdirs[h] == 4 || subdirs[h] == 5) subdirx = -0.5f;
					histflowX += (float)dirsX[h] + subdirx;
					meancount_x++;
					float subdiry = 0.0f;
					if (subdirs[h] == 5 || subdirs[h] == 6 || subdirs[h] == 7) subdiry = -0.5f;
					if (subdirs[h] == 1 || subdirs[h] == 2 || subdirs[h] == 3) subdiry = 0.5f;
					histflowY += (float)dirsy[h] + subdiry;
					meancount_y++;
				}
				histflowX /= meancount_x;
				histflowY /= meancount_y;
			}
			

			///*-----------从这里是角速度补偿，没有旋转数据就不用看------------------*/
			///* compensate rotation */
			///* calculate focal length in pixels */
			//const float focal_length_px = (global_data.param[PARAM_FOCAL_LENGTH_MM]) / (4.0f * 6.0f) * 1000.0f; //original focal lenght: 12mm pixelsize: 6um, binning 4 enabled
			//																									/*
			//																									* gyro compensation
			//																									* the compensated value is clamped to
			//																									* the maximum measurable flow value (PARAM_MAX_FLOW_PIXEL) +0.5
			//																									* (sub pixel flow can add half pixel to the value)
			//																									*
			//																									* -y_rate gives x flow
			//																									* x_rates gives y_flow
			//																									*/
			//if (FLOAT_AS_BOOL(global_data.param[PARAM_BOTTOM_FLOW_GYRO_COMPENSATION])) {//角速度补偿，需要用到陀螺仪的数据
			//	if (fabsf(y_rate) > global_data.param[PARAM_GYRO_COMPENSATION_THRESHOLD]) {
			//		/* calc pixel of gyro */
			//		float y_rate_pixel = y_rate * (get_time_between_images() / 1000000.0f) * focal_length_px;
			//		float comp_x = histflowX + y_rate_pixel;
			//		/* clamp value to maximum search window size plus half pixel from subpixel search */
			//		if (comp_x < (-SEARCH_SIZE - 0.5f))
			//			*pixel_flow_x = (-SEARCH_SIZE - 0.5f);
			//		else if (comp_x >(SEARCH_SIZE + 0.5f))
			//			*pixel_flow_x = (SEARCH_SIZE + 0.5f);
			//		else
			//			*pixel_flow_x = comp_x;
			//	} else {
			//		*pixel_flow_x = histflowX;
			//	}
			//	if (fabsf(x_rate) > global_data.param[PARAM_GYRO_COMPENSATION_THRESHOLD]) {
			//		/* calc pixel of gyro */
			//		float x_rate_pixel = x_rate * (get_time_between_images() / 1000000.0f) * focal_length_px;
			//		float comp_y = histflowY - x_rate_pixel;
			//		/* clamp value to maximum search window size plus/minus half pixel from subpixel search */
			//		if (comp_y < (-SEARCH_SIZE - 0.5f))
			//			*pixel_flow_y = (-SEARCH_SIZE - 0.5f);
			//		else if (comp_y >(SEARCH_SIZE + 0.5f))
			//			*pixel_flow_y = (SEARCH_SIZE + 0.5f);
			//		else
			//			*pixel_flow_y = comp_y;
			//	} else {
			//		*pixel_flow_y = histflowY;
			//	}
			//	/* alternative compensation */
			//	/* compensate y rotation */
			//	/* *pixel_flow_x = histflowX + y_rate_pixel;*/
			//	/* compensate x rotation */
			//	/* *pixel_flow_y = histflowY - x_rate_pixel;*/
			//} else {//角速度补偿
			//		/* without gyro compensation */
			//	*pixel_flow_x = histflowX;
			//	*pixel_flow_y = histflowY;
			//}
			///*-----------角速度补偿一直到这里--------------*/


		} else {// if(1)
			*pixel_flow_x = 0.0f;
			*pixel_flow_y = 0.0f;
			return 0;
		}
	} else {//if (meancount > 10)
		*pixel_flow_x = 0.0f;
		*pixel_flow_y = 0.0f;
		return 0;
	}
	/* calculate quality */
	uint8_t qual = (uint8_t)(meancount * 255 / (NUM_BLOCKS*NUM_BLOCKS));
	return qual;
}