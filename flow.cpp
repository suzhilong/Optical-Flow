#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;

int SEARCH_SIZE = 5;
int TILE_SIZE = 5;//5*5
int THRESHOLD_DIFF = 40;
int THRESHOLD_SAD = 5000;
bool HIST_FILTER = true;

int compute_diff(const cv::Mat &img1, int row, int col);
int comput_SAD(const cv::Mat &img1, const cv::Mat &img2, int row1, int col1, int row2, int col2);
void compute_subpixel(const cv::Mat &img1, const cv::Mat &img2, int row1, int col1, int row2, int col2, vector<int> &acc);
float compute_flow(const cv::Mat &img1, const cv::Mat &img2, float &dx, float &dy);



int main()
{
	cv::Mat img1, img2;
	std::string path1 = "D:\\szl\\pic\\image\\100.png";
	std::string path2 = "D:\\szl\\pic\\image\\101.png";

	img1 = cv::imread(path1, CV_LOAD_IMAGE_GRAYSCALE);
	img2 = cv::imread(path2, CV_LOAD_IMAGE_GRAYSCALE);

	cout << "size of img1: " << img1.rows << " * " << img1.cols << endl;
	cout << "size of img2: " << img2.rows << " * " << img2.cols << endl;

	cout << "img1:" << endl;
	for(int i=0;i<20;i++)
		cout << int(img1.at<uchar>(400,100+i)) << " ";
	cout << endl;
	cout << "img2:" << endl;
	for (int i = 0; i<20; i++)
		cout << int(img2.at<uchar>(400, 100 + i)) << " ";
	cout << endl;
	
	/*cv::imshow("im1", img1);
	cv::imshow("im2", img2);
	cv::waitKey(0);*/

	float dx = 0, dy = 0;
	float qual = compute_flow(img1, img2, dx, dy);

	cout << dx << dy << endl;

	system("pause");
	return 0;
}

int compute_diff(const cv::Mat &img1, int row, int col)
{
	int diff = 0;
	int winmin = row - (TILE_SIZE / 2);
	int winmax = row + (TILE_SIZE / 2);

	if (winmin < 0 || winmax > img1.cols || winmax>img1.rows)
		return THRESHOLD_DIFF;

	for (int r = winmin; r <= winmax; r++) {
		for (int c = winmin + 1; c <= winmax; c++) {
			diff += abs(img1.at<uchar>(r, c) - img1.at<uchar>(r , c-1));
		}
	}
	for (int r = winmin + 1; r <= winmax; r++) {
		for (int c = winmin; c <= winmax; c++) {
			diff += abs(img1.at<uchar>(r, c) - img1.at<uchar>(r-1, c));
		}
	}
	return diff;
}

int comput_SAD(const cv::Mat &img1, const cv::Mat &img2, int row1, int col1, int row2, int col2)
{
	int winmin1 = row1 - SEARCH_SIZE, winmax1 = row1 + SEARCH_SIZE;
	int winmin2 = row2 - SEARCH_SIZE, winmax2 = row2 + SEARCH_SIZE;

	if (winmin1<0 || winmin2<0 || winmax1>img1.cols || winmax1>img1.rows || winmax2>img2.cols || winmax2>img2.rows)
		return INT_MAX;

	int sad = 0;
	for (int i = winmin1; i <= winmax1; i++) {
		for (int j = winmin2; j <= winmax2; j++) {
			sad = abs(img1.at<uchar>(j, i) - img2.at<uchar>(j, i));
		}
	}
	return sad;
}

void compute_subpixel(const cv::Mat &img1, const cv::Mat &img2, int row1, int col1, int row2, int col2, vector<int> &acc)
{
	/*	+	  -		  +  -  +     -	    +
	*	|			 s5    s7			|
	*	+  leftip	 t5 s6 t7  rightup  +
	*	|			 s4  X s0			|
	*	+ leftdown   t3 s2 t1 rightdown +
	*	|		     s3	   s1		    |
	*	+	  -		 +   -  +	 -	    + */
	int s0_1, s1_1, t1_1, s2_1, s3_1, t3_1, s4_1, s5_1, t5_1, s6_1, s7_1, t7_1;
	s0_1 = (img1.at<uchar>(row1, col1) + img1.at<uchar>(row1, col1 + 1)) / 2;
	s1_1 = (img1.at<uchar>(row1 + 1, col1) + img1.at<uchar>(row1 + 1, col1 + 1)) / 2;
	s2_1 = (img1.at<uchar>(row1, col1) + img1.at<uchar>(row1 + 1, col1)) / 2;
	s3_1 = (img1.at<uchar>(row1 + 1, col1) + img1.at<uchar>(row1 + 1, col1 - 1)) / 2;
	s4_1 = (img1.at<uchar>(row1, col1) + img1.at<uchar>(row1, col1 - 1)) / 2;
	s5_1 = (img1.at<uchar>(row1 - 1, col1) + img1.at<uchar>(row1 - 1, col1 - 1)) / 2;
	s6_1 = (img1.at<uchar>(row1, col1) + img1.at<uchar>(row1 - 1, col1)) / 2;
	s7_1 = (img1.at<uchar>(row1 - 1, col1) + img1.at<uchar>(row1 - 1, col1 + 1)) / 2;
	t1_1 = (s0_1 + s1_1) / 2;
	t3_1 = (s4_1 + s3_1) / 2;
	t5_1 = (s4_1 + s5_1) / 2;
	t7_1 = (s0_1 + s7_1) / 2;

	int s0_2, s1_2, t1_2, s2_2, s3_2, t3_2, s4_2, s5_2, t5_2, s6_2, s7_2, t7_2;
	s0_2 = (img2.at<uchar>(row2, col2) + img2.at<uchar>(row2, col2 + 1)) / 2;
	s1_2 = (img2.at<uchar>(row2 + 1, col2) + img2.at<uchar>(row2 + 1, col2 + 1)) / 2;
	s2_2 = (img2.at<uchar>(row2, col2) + img2.at<uchar>(row2 + 1, col2)) / 2;
	s3_2 = (img2.at<uchar>(row2 + 1, col2) + img2.at<uchar>(row2 + 1, col2 - 1)) / 2;
	s4_2 = (img2.at<uchar>(row2, col2) + img2.at<uchar>(row2, col2 - 1)) / 2;
	s5_2 = (img2.at<uchar>(row2 - 1, col2) + img2.at<uchar>(row2 - 1, col2 - 1)) / 2;
	s6_2 = (img2.at<uchar>(row2, col2) + img2.at<uchar>(row2 - 1, col2)) / 2;
	s7_2 = (img2.at<uchar>(row2 - 1, col2) + img2.at<uchar>(row1 - 1, col2 + 1)) / 2;
	t1_2 = (s0_2 + s1_2) / 2;
	t3_2 = (s4_2 + s3_2) / 2;
	t5_2 = (s4_2 + s5_2) / 2;
	t7_2 = (s0_2 + s7_2) / 2;

	int rightup_1 = (img1.at<uchar>(row1 - 1, col1 + 1) + img1.at<uchar>(row1, col1 + 1)) / 2;
	int rightdown_1 = (img1.at<uchar>(row1 + 1, col1 + 1) + img1.at<uchar>(row1, col1 + 1)) / 2;
	int leftup_1 = (img1.at<uchar>(row1 - 1, col1 - 1) + img1.at<uchar>(row1, col1 - 1)) / 2;
	int leftdown_1 = (img1.at<uchar>(row1 + 1, col1 - 1) + img1.at<uchar>(row1, col1 - 1)) / 2;
	int rightup_2 = (img2.at<uchar>(row2 - 1, col2 + 1) + img2.at<uchar>(row2, col2 + 1)) / 2;
	int rightdown_2 = (img2.at<uchar>(row2 + 1, col1 + 1) + img2.at<uchar>(row2, col1 + 1)) / 2;
	int leftup_2 = (img2.at<uchar>(row2 - 1, col1 - 1) + img2.at<uchar>(row2, col1 - 1)) / 2;
	int leftdown_2 = (img2.at<uchar>(row2 + 1, col1 - 1) + img2.at<uchar>(row2, col1 - 1)) / 2;

	acc[0] = abs(s6_1 - s6_2) + abs(t7_1 - t7_2) + abs(rightup_1 - rightdown_2) +
		abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2)) + abs(s0_1 - s0_2) + abs(t1_1 - t1_2) + abs(img1.at<uchar>(row1, col1 + 1) - img2.at<uchar>(row2, col2 + 1)) +
		abs(s2_1 - s2_2) + abs(t1_1 - t1_2) + abs(rightdown_1 - rightdown_2);
	acc[1] = abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2)) + abs(s0_1 - s0_2) + abs(img1.at<uchar>(row1, col1 + 1) - img2.at<uchar>(row2, col2 + 1)) +
		abs(s2_1 - s2_2) + abs(t1_1 - t1_2) + abs(rightdown_1 - rightdown_2) +
		abs(img1.at<uchar>(row1 + 1, col1) - img2.at<uchar>(row2 + 1, col2)) + abs(s1_1 - s1_2) + abs(img1.at<uchar>(row1 + 1, col1 + 1) - img2.at<uchar>(row2 + 1, col2 + 1));
	acc[2] = abs(s4_1 - s4_2) + abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2)) + abs(s0_1 - s0_2) +
		abs(t3_1 - t3_2) + abs(s2_1 - s2_2) + abs(t1_1 - t1_2) +
		abs(s3_1 - s3_2) + abs(img1.at<uchar>(row1 + 1, col1) - img2.at<uchar>(row2 + 1, col2)) + abs(s1_1 - s1_2);
	acc[3] = abs(img1.at<uchar>(row1, col1 - 1) - img2.at<uchar>(row2, col2 - 1)) + abs(s4_1 - s4_2) + abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2)) +
		abs(leftdown_1 - leftdown_2) + abs(t3_1 - t3_2) + abs(s2_1 - s2_2) +
		abs(img1.at<uchar>(row1 + 1, col1 - 1) - img2.at<uchar>(row2 + 1, col2 - 1)) + abs(s3_1 - s3_2) + abs(img1.at<uchar>(row1 + 1, col1) - img2.at<uchar>(row2 + 1, col2));
	acc[4] = abs(leftup_1 - leftdown_2) + abs(t5_1 - t5_2) + abs(s6_1 - s6_2) +
		abs(img1.at<uchar>(row1, col1 - 1) - img2.at<uchar>(row2, col2 - 1)) + abs(s4_1 - s4_2) + abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2)) +
		abs(leftdown_1 - leftdown_2) + abs(t3_1 - t3_2) + abs(s2_1 - s2_2);
	acc[5] = abs(img1.at<uchar>(row1 - 1, col1 - 1) - img2.at<uchar>(row2 - 1, col2 - 1)) + abs(s5_1 - s5_2) + abs(img1.at<uchar>(row1 - 1, col1) - img2.at<uchar>(row2 - 1, col2)) +
		abs(leftup_1 - leftdown_2) + abs(t5_1 - t5_2) + abs(s6_1 - s6_2) +
		abs(img1.at<uchar>(row1, col1 - 1) - img2.at<uchar>(row2, col2 - 1)) + abs(s4_1 - s4_2) + abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2));
	acc[6] = abs(s5_1 - s5_2) + abs(img1.at<uchar>(row1 - 1, col1) - img2.at<uchar>(row2 - 1, col2)) + abs(s7_1 - s7_2) +
		abs(t5_1 - t5_2) + abs(s6_1 - s6_2) + abs(t7_1 - t7_2) +
		abs(s4_1 - s4_2) + abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2)) + abs(s0_1 - s0_2);
	acc[7] = abs(img1.at<uchar>(row1 - 1, col1) - img2.at<uchar>(row2 - 1, col2)) + abs(s7_1 - s7_2) + abs(img1.at<uchar>(row1 - 1, col1 + 1) - img2.at<uchar>(row2 - 1, col2 + 1)) +
		abs(s6_1 - s6_2) + abs(t7_1 - t7_2) + abs(rightup_1 - rightup_2) +
		abs(img1.at<uchar>(row1, col1) - img2.at<uchar>(row2, col2)) + abs(s0_1 - s0_2) + abs(img1.at<uchar>(row1, col1 + 1) - img2.at<uchar>(row2, col2 + 1));

}


float compute_flow(const cv::Mat &img1, const cv::Mat &img2, float &dx, float &dy)
{
	int NUM_BLOCKS_WIDTH = img1.cols / TILE_SIZE;
	int NUM_BLOCKS_HEIGHT = img1.rows / TILE_SIZE;
	int hist_size = 2 * (SEARCH_SIZE * 2 + 1) + 1;
	int meancount = 0;
	float meanflowX = 0.0f;
	float meanflowY = 0.0f;
	float histflowX = 0.0f;
	float histflowY = 0.0f;

	vector<int> acc(8);//8个半像素方向的SAD值
	vector<int> dirsX(SEARCH_SIZE*SEARCH_SIZE,0);//每个特征点x方向的位移
	vector<int> dirsY(SEARCH_SIZE*SEARCH_SIZE,0);//每个特征点y方向上位移
	vector<int> subdirs(SEARCH_SIZE*SEARCH_SIZE,0);//每个特征点半像素的方向
	vector<int> histX(hist_size,0);//x方向位移（半像素）的统计
	vector<int> histY(hist_size,0);//y方向位移（半像素）的统计


	for (int row1 = TILE_SIZE/2; row1 < img1.rows; row1 += TILE_SIZE) {
		for (int col1 = TILE_SIZE/2; col1 < img1.cols; col1 += TILE_SIZE) {
			int diff = compute_diff(img1, row1, col1);
			if (diff < THRESHOLD_DIFF) continue;

			int dist = INT_MAX;
			int sumX = 0;
			int sumY = 0;
			for (int row2 = max(0,row1 - SEARCH_SIZE); row2 < min(row1 + SEARCH_SIZE,img2.rows); row2++) {
				for (int col2 = max(0, row1 - SEARCH_SIZE); col2 < min(col1 + SEARCH_SIZE,img2.cols); col2++) {
					int temp_dist = comput_SAD(img1, img2, row1, col1, row2, col2);
					if (temp_dist < dist) {
						sumX = col2 - col1;
						sumY = row2 - row1;
						dist = temp_dist;
					}
				}
			}

			float meanflowX = 0, meanflowY = 0;
			if (dist < THRESHOLD_SAD) {
				meanflowX += (float)sumX;
				meanflowY += (float)sumY;

				compute_subpixel(img1, img2, row1, col1, row1 + sumX, col1 + sumY, acc);
				int mindist = dist;
				int mindir = 8;
				for (int k = 0; k < 8; k++) {
					if (acc[k] < mindist) {
						// SAD becomes better in direction k
						mindist = acc[k];
						mindir = k;
					}
				}
				dirsX[meancount] = sumX;
				dirsY[meancount] = sumY;
				subdirs[meancount] = mindir;
				meancount++;

				/* 
				* histX的索引：
				* 因为半像素，所以偏移量sumX要*2；histX索引 从winmin到winmax，所以要有一个(SEARCH_SIZE*2)的偏移。
				*
				* 比如 sumX = -5 个像素，偏移到了最左边边缘的像素，
				* 对应的 histX 的 hist_index_x 为 1 = 2 * -5 + 5*2+1。
				* 当mindir即半像素位置偏左半个像素时，对应的 histX 的 hist_index_x 就为 0
				*
				*/
				int hist_index_x = 2 * sumX + (SEARCH_SIZE*2 + 1);
				if (mindir == 0 || mindir == 1 || mindir == 7) hist_index_x += 1;
				if (mindir == 3 || mindir == 4 || mindir == 5) hist_index_x += -1;
				int hist_index_y = 2 * sumY + (SEARCH_SIZE*2 + 1);
				if (mindir == 5 || mindir == 6 || mindir == 7) hist_index_y += -1;
				if (mindir == 1 || mindir == 2 || mindir == 3) hist_index_y += 1;
				histX[hist_index_x]++;
				histY[hist_index_y]++;
			}
		}
	}

	if (meancount > 10) {
		meanflowX /= meancount;
		meanflowY /= meancount;

		int maxpositionX = 0;
		int maxpositionY = 0;
		int maxvalueX = 0;
		int maxvalueY = 0;

		for (int j = 0; j < hist_size; j++) {
			if (histX[j] > maxvalueX) {
				maxvalueX = histX[j];
				maxpositionX = j;
			}
			if (histY[j] > maxvalueY) {
				maxvalueY = histY[j];
				maxpositionY = j;
			}
		}

		if (1) 
		{
			if (HIST_FILTER) {/*直方图均衡*/
				int hist_x_min = maxpositionX;
				int hist_x_max = maxpositionX;
				int hist_y_min = maxpositionY;
				int hist_y_max = maxpositionY;

				/* x direction */
				if (maxpositionX > 1 && maxpositionX < hist_size - 2) {
					hist_x_min = maxpositionX - 2;
					hist_x_max = maxpositionX + 2;
				} else if (maxpositionX == 0) {
					hist_x_min = maxpositionX;
					hist_x_max = maxpositionX + 2;
				} else  if (maxpositionX == hist_size - 1) {
					hist_x_min = maxpositionX - 2;
					hist_x_max = maxpositionX;
				} else if (maxpositionX == 1) {
					hist_x_min = maxpositionX - 1;
					hist_x_max = maxpositionX + 2;
				} else  if (maxpositionX == hist_size - 2) {
					hist_x_min = maxpositionX - 2;
					hist_x_max = maxpositionX + 1;
				}

				/* y direction */
				if (maxpositionY > 1 && maxpositionY < hist_size - 2) {
					hist_y_min = maxpositionY - 2;
					hist_y_max = maxpositionY + 2;
				} else if (maxpositionY == 0) {
					hist_y_min = maxpositionY;
					hist_y_max = maxpositionY + 2;
				} else if (maxpositionY == hist_size - 1) {
					hist_y_min = maxpositionY - 2;
					hist_y_max = maxpositionY;
				} else if (maxpositionY == 1) {
					hist_y_min = maxpositionY - 1;
					hist_y_max = maxpositionY + 2;
				} else if (maxpositionY == hist_size - 2) {
					hist_y_min = maxpositionY - 2;
					hist_y_max = maxpositionY + 1;
				}

				float hist_x_value = 0.0f;
				float hist_x_weight = 0.0f;

				float hist_y_value = 0.0f;
				float hist_y_weight = 0.0f;

				for (int h = hist_x_min; h < hist_x_max + 1; h++) {
					hist_x_value += (float)(h*histX[h]);
					hist_x_weight += (float)histX[h];
				}

				for (int h = hist_y_min; h<hist_y_max + 1; h++) {
					hist_y_value += (float)(h*histY[h]);
					hist_y_weight += (float)histY[h];
				}

				histflowX = (hist_x_value / hist_x_weight - (SEARCH_SIZE*2 + 1)) / 2.0f;
				histflowY = (hist_y_value / hist_y_weight - (SEARCH_SIZE*2 + 1)) / 2.0f;

			} else {/*平均法*/
				int meancount_x = 0;
				int meancount_y = 0;

				for (int h = 0; h < meancount; h++) {
					float subdirX = 0.0f;
					if (subdirs[h] == 0 || subdirs[h] == 1 || subdirs[h] == 7) subdirX = 0.5f;
					if (subdirs[h] == 3 || subdirs[h] == 4 || subdirs[h] == 5) subdirX = -0.5f;
					histflowX += (float)dirsX[h] + subdirX;
					meancount_x++;

					float subdirY = 0.0f;
					if (subdirs[h] == 5 || subdirs[h] == 6 || subdirs[h] == 7) subdirY = -0.5f;
					if (subdirs[h] == 1 || subdirs[h] == 2 || subdirs[h] == 3) subdirY = 0.5f;
					histflowY += (float)dirsY[h] + subdirY;
					meancount_y++;
				}

				histflowX /= meancount_x;
				histflowY /= meancount_y;

			}

		} else {// if(1)
			dx = 0.0f;
			dy = 0.0f;
			return 0;
		}
	} else {//if (meancount > 10)
		dx = 0.0f;
		dy = 0.0f;
		return 0;
	}

	
	float qual = (float)(meancount * 255 / (NUM_BLOCKS_WIDTH*NUM_BLOCKS_HEIGHT));
	
	return qual;
}