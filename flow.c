uint8_t compute_flow(uint8_t *image1, uint8_t *image2, float x_rate, float y_rate, float z_rate, float *pixel_flow_x, float *pixel_flow_y) {
	/* constants */
	const int16_t winmin = -SEARCH_SIZE;
	const int16_t winmax = SEARCH_SIZE;
	const uint16_t hist_size = 2*(winmax-winmin+1)+1;
	/* variables */
	uint16_t pixLo = SEARCH_SIZE + 1;//图片最左边一块中间的像素
	uint16_t pixHi = FRAME_SIZE - (SEARCH_SIZE + 1) - TILE_SIZE;//图片最右边一块中间的像素 //FRAME_SIZE是图片宽度
	uint16_t pixStep = (pixHi - pixLo) / NUM_BLOCKS + 1;
	uint16_t i, j;
	uint32_t acc[8];           // subpixels
	uint16_t histx[hist_size]; // counter for x shift
	uint16_t histy[hist_size]; // counter for y shift
	int8_t  dirsx[64];         // shift directions in x
	int8_t  dirsy[64];         // shift directions in y
	uint8_t  subdirs[64];      // shift directions of best subpixels
	float meanflowx = 0.0f;
	float meanflowy = 0.0f;
	uint16_t meancount = 0;
	float histflowx = 0.0f;
	float histflowy = 0.0f;
	/* initialize with 0 */
	for (j = 0; j < hist_size; j++) { histx[j] = 0; histy[j] = 0; }
	/* iterate over all patterns*/
	for (j = pixLo; j < pixHi; j += pixStep)
	{
		for (i = pixLo; i < pixHi; i += pixStep)
		{
			/* test pixel if it is suitable for flow tracking */
			uint32_t diff = compute_diff(image1, i, j, (uint16_t) global_data.param[PARAM_IMAGE_WIDTH]);//计算SAD
			if (diff < global_data.param[PARAM_BOTTOM_FLOW_FEATURE_THRESHOLD])//SAD太小，不满足特征点的要求
			{
				continue;
			}
			uint32_t dist = 0xFFFFFFFF; // set initial distance to "infinity"
			int8_t sumx = 0;
			int8_t sumy = 0;
			int8_t ii, jj;
			uint8_t *base1 = image1 + j * (uint16_t) global_data.param[PARAM_IMAGE_WIDTH] + i;
			for (jj = winmin; jj <= winmax; jj++)
			{//得到winmin~winmax内最相似的偏移量
				uint8_t *base2 = image2 + (j+jj) * (uint16_t) global_data.param[PARAM_IMAGE_WIDTH] + i;
				for (ii = winmin; ii <= winmax; ii++)
				{
//					uint32_t temp_dist = compute_sad_8x8(image1, image2, i, j, i + ii, j + jj, (uint16_t) global_data.param[PARAM_IMAGE_WIDTH]);
					uint32_t temp_dist = ABSDIFF(base1, base2 + ii);
					/*
					int absdiff(unsigned char *base1,unsigned char *base2)
					{
						int result = 0;
						for(int i = 0;i <8;i++)
						{
							for(j = 0;j < 8;j++)
							{
								result += abs(base1[i*64 + j] - base2[i*64 + j]);
							}
						}
						return result;
					}
					*/
					if (temp_dist < dist)
					{
						sumx = ii;//x偏移量
						sumy = jj;//y偏移量
						dist = temp_dist;
					}
				}
			}
			/* acceptance SAD distance threshhold */
			if (dist < global_data.param[PARAM_BOTTOM_FLOW_VALUE_THRESHOLD])
			{
				meanflowx += (float) sumx;
				meanflowy += (float) sumy;
				//8个方向的半像素寻找最佳匹配，把八个方向的半像素SAD保存到acc
				compute_subpixel(image1, image2, i, j, i + sumx, j + sumy, acc, (uint16_t) global_data.param[PARAM_IMAGE_WIDTH]);
				uint32_t mindist = dist; // best SAD until now
				uint8_t mindir = 8; // direction 8 for no direction
				for(uint8_t k = 0; k < 8; k++)
				{
					if (acc[k] < mindist)
					{
						// SAD becomes better in direction k
						mindist = acc[k];
						mindir = k;
					}
				}
				dirsx[meancount] = sumx;
				dirsy[meancount] = sumy;
				subdirs[meancount] = mindir;
				meancount++;//记录达到SAD要求的特征点个数
				/* feed histogram filter*/
				uint8_t hist_index_x = 2*sumx + (winmax-winmin+1);
				if (mindir == 0 || mindir == 1 || mindir == 7) hist_index_x += 1;
				if (mindir == 3 || mindir == 4 || mindir == 5) hist_index_x += -1;
				uint8_t hist_index_y = 2*sumy + (winmax-winmin+1);
				if (mindir == 5 || mindir == 6 || mindir == 7) hist_index_y += -1;
				if (mindir == 1 || mindir == 2 || mindir == 3) hist_index_y += 1;
				// 4个方向的直方图
				histx[hist_index_x]++;
				histy[hist_index_y]++;
			}
		}
	}
	/* create flow image if needed (image1 is not needed anymore)
	 * -> can be used for debugging purpose
	 */
	if (FLOAT_AS_BOOL(global_data.param[PARAM_USB_SEND_VIDEO]))//&& global_data.param[PARAM_VIDEO_USB_MODE] == FLOW_VIDEO)
	{
		for (j = pixLo; j < pixHi; j += pixStep)
		{
			for (i = pixLo; i < pixHi; i += pixStep)
			{
				uint32_t diff = compute_diff(image1, i, j, (uint16_t) global_data.param[PARAM_IMAGE_WIDTH]);
				if (diff > global_data.param[PARAM_BOTTOM_FLOW_FEATURE_THRESHOLD])
				{
					image1[j * ((uint16_t) global_data.param[PARAM_IMAGE_WIDTH]) + i] = 255;
				}
			}
		}
	}
	/* evaluate flow calculation */
	if (meancount > 10)
	{
		meanflowx /= meancount;
		meanflowy /= meancount;
		int16_t maxpositionx = 0;
		int16_t maxpositiony = 0;
		uint16_t maxvaluex = 0;
		uint16_t maxvaluey = 0;
		/* position of maximal histogram peek */
		for (j = 0; j < hist_size; j++)
		{
			if (histx[j] > maxvaluex)
			{
				maxvaluex = histx[j];
				maxpositionx = j;
			}
			if (histy[j] > maxvaluey)
			{
				maxvaluey = histy[j];
				maxpositiony = j;
			}
		}
		/* check if there is a peak value in histogram */
		if (1) //(histx[maxpositionx] > meancount / 6 && histy[maxpositiony] > meancount / 6)
		{
			if (FLOAT_AS_BOOL(global_data.param[PARAM_BOTTOM_FLOW_HIST_FILTER]))
			{//滤波法
				/* use histogram filter peek value */
				uint16_t hist_x_min = maxpositionx;
				uint16_t hist_x_max = maxpositionx;
				uint16_t hist_y_min = maxpositiony;
				uint16_t hist_y_max = maxpositiony;
				/* x direction */
				if (maxpositionx > 1 && maxpositionx < hist_size-2)
				{
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx + 2;
				}
				else if (maxpositionx == 0)
				{
					hist_x_min = maxpositionx;
					hist_x_max = maxpositionx + 2;
				}
				else  if (maxpositionx == hist_size-1)
				{
					hist_x_min = maxpositionx - 2;  
					hist_x_max = maxpositionx;
				}
				else if (maxpositionx == 1)
				{
					hist_x_min = maxpositionx - 1;
					hist_x_max = maxpositionx + 2;
				}
				else  if (maxpositionx == hist_size-2)
				{
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx + 1;
				}
				/* y direction */
				if (maxpositiony > 1 && maxpositiony < hist_size-2)
				{
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony + 2;
				}
				else if (maxpositiony == 0)
				{
					hist_y_min = maxpositiony;
					hist_y_max = maxpositiony + 2;
				}
				else if (maxpositiony == hist_size-1)
				{
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony;
				}
				else if (maxpositiony == 1)
				{
					hist_y_min = maxpositiony - 1;
					hist_y_max = maxpositiony + 2;
				}
				else if (maxpositiony == hist_size-2)
				{
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony + 1;
				}
				float hist_x_value = 0.0f;
				float hist_x_weight = 0.0f;
				float hist_y_value = 0.0f;
				float hist_y_weight = 0.0f;
				for (uint8_t h = hist_x_min; h < hist_x_max+1; h++)
				{
					hist_x_value += (float) (h*histx[h]);
					hist_x_weight += (float) histx[h];
				}
				for (uint8_t h = hist_y_min; h<hist_y_max+1; h++)
				{
					hist_y_value += (float) (h*histy[h]);
					hist_y_weight += (float) histy[h];
				}
				histflowx = (hist_x_value/hist_x_weight - (winmax-winmin+1)) / 2.0f ;
				histflowy = (hist_y_value/hist_y_weight - (winmax-winmin+1)) / 2.0f;
			}
			else
			{//平均法
				/* use average of accepted flow values */
				uint32_t meancount_x = 0;
				uint32_t meancount_y = 0;
				for (uint8_t h = 0; h < meancount; h++)
				{
					float subdirx = 0.0f;
					if (subdirs[h] == 0 || subdirs[h] == 1 || subdirs[h] == 7) subdirx = 0.5f;
					if (subdirs[h] == 3 || subdirs[h] == 4 || subdirs[h] == 5) subdirx = -0.5f;
					histflowx += (float)dirsx[h] + subdirx;
					meancount_x++;
					float subdiry = 0.0f;
					if (subdirs[h] == 5 || subdirs[h] == 6 || subdirs[h] == 7) subdiry = -0.5f;
					if (subdirs[h] == 1 || subdirs[h] == 2 || subdirs[h] == 3) subdiry = 0.5f;
					histflowy += (float)dirsy[h] + subdiry;
					meancount_y++;
				}
				histflowx /= meancount_x;
				histflowy /= meancount_y;
			}
/*-----------从这里是角速度补偿，没有旋转数据就不用看------------------*/
			/* compensate rotation */
			/* calculate focal length in pixels */
			const float focal_length_px = (global_data.param[PARAM_FOCAL_LENGTH_MM]) / (4.0f * 6.0f) * 1000.0f; //original focal lenght: 12mm pixelsize: 6um, binning 4 enabled
			/*
			 * gyro compensation
			 * the compensated value is clamped to
			 * the maximum measurable flow value (PARAM_MAX_FLOW_PIXEL) +0.5
			 * (sub pixel flow can add half pixel to the value)
			 *
			 * -y_rate gives x flow
			 * x_rates gives y_flow
			 */
			if (FLOAT_AS_BOOL(global_data.param[PARAM_BOTTOM_FLOW_GYRO_COMPENSATION]))
			{//角速度补偿，需要用到陀螺仪的数据
				if(fabsf(y_rate) > global_data.param[PARAM_GYRO_COMPENSATION_THRESHOLD])
				{
					/* calc pixel of gyro */
					float y_rate_pixel = y_rate * (get_time_between_images() / 1000000.0f) * focal_length_px;
					float comp_x = histflowx + y_rate_pixel;
					/* clamp value to maximum search window size plus half pixel from subpixel search */
					if (comp_x < (-SEARCH_SIZE - 0.5f))
						*pixel_flow_x = (-SEARCH_SIZE - 0.5f);
					else if (comp_x > (SEARCH_SIZE + 0.5f))
						*pixel_flow_x = (SEARCH_SIZE + 0.5f);
					else
						*pixel_flow_x = comp_x;
				}
				else
				{
					*pixel_flow_x = histflowx;
				}
				if(fabsf(x_rate) > global_data.param[PARAM_GYRO_COMPENSATION_THRESHOLD])
				{
					/* calc pixel of gyro */
					float x_rate_pixel = x_rate * (get_time_between_images() / 1000000.0f) * focal_length_px;
					float comp_y = histflowy - x_rate_pixel;
					/* clamp value to maximum search window size plus/minus half pixel from subpixel search */
					if (comp_y < (-SEARCH_SIZE - 0.5f))
						*pixel_flow_y = (-SEARCH_SIZE - 0.5f);
					else if (comp_y > (SEARCH_SIZE + 0.5f))
						*pixel_flow_y = (SEARCH_SIZE + 0.5f);
					else
						*pixel_flow_y = comp_y;
				}
				else
				{
					*pixel_flow_y = histflowy;
				}
					/* alternative compensation */
					/* compensate y rotation */
					/* *pixel_flow_x = histflowx + y_rate_pixel;*/
	
					/* compensate x rotation */
					/* *pixel_flow_y = histflowy - x_rate_pixel;*/
			}
			else
			{//角速度补偿
				/* without gyro compensation */
				*pixel_flow_x = histflowx;
				*pixel_flow_y = histflowy;
			}
/*-----------角速度补偿一直到这里--------------*/
		}
		else
		{// if(1)
			*pixel_flow_x = 0.0f;
			*pixel_flow_y = 0.0f;
			return 0;
		}
	}
	else
	{//if (meancount > 10)
		*pixel_flow_x = 0.0f;
		*pixel_flow_y = 0.0f;
		return 0;
	}
	/* calculate quality */
	uint8_t qual = (uint8_t)(meancount * 255 / (NUM_BLOCKS*NUM_BLOCKS));
	return qual;
}