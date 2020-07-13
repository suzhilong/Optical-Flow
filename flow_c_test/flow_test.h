#ifndef FLOW_H_
#define FLOW_H_

#include <stdint.h>

/**
* @brief Computes pixel flow from image1 to image2
*/
uint8_t computeFlow(uint8_t *image1, uint8_t *image2, float *pixel_flow_x, float *pixel_flow_y);

#endif
