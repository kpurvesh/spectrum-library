/*
 * dlpspec_scan_chemo.c
 *
 *  Created on: Sep 6, 2016
 *      Author: Purvesh
 */
#include <string.h>
#include "dlpspec_util.h"
#include "dlpspec_scan_chemo.h"

DLPSPEC_ERR_CODE dlpspec_scan_chemo_interpret(const uScanData *puScanData, scanResults *pResults)
{
	patDefChemo patDef;
	int i;
	chemoScanConfig cfg;
	const chemoScanData *pScanData = &puScanData->chemo_data;

	dlpspec_chemo_subtract_remove_dc_level(pScanData, pResults);
	dlpspec_copy_chemoScanData_hdr_to_scanResults(pScanData, pResults);

//	for(i=0; i < pScanData->chemoCfg.head.num_sections; i++)
//	{
//		cfg.head.scan_type = pScanData->chemoCfg.head.scan_type;
//		cfg.head.scanConfigIndex = pScanData->chemoCfg.head.scanConfigIndex;
//		cfg.chemoSection[i].num_patterns = pScanData->chemoCfg.chemoSection[i].num_patterns;
//		cfg.head.num_repeats = pScanData->chemoCfg.head.num_repeats;
//		for(i = 0; i < pScanData->chemoCfg.chemoSection[i].num_patterns; i++) {
//			cfg.chemoSection[i].heights_px[i] = pScanData->chemoCfg.chemoSection[i].heights_px[i];
//			cfg.chemoSection[i].wavelengths[i] = pScanData->chemoCfg.chemoSection[i].wavelengths[i];
//		}
//	}

	for(i=0; i < pResults->length; i++)
	{
		pResults->wavelength[i] = i+1;
	}

	pResults->pga = pScanData->pga;
	return DLPSPEC_PASS;
}

uint8_t dlpspec_scan_chemo_genPatDef(const chemoScanConfig *pScanConfig, const calibCoeffs *pCoeffs,
			patDefChemo *patDefCh, const FrameBufferDescriptor *pFB)
{
	int i,j;
	double mid_x;
	int32_t ret;

	DLPSPEC_ERR_CODE ret_val = (DLPSPEC_PASS);
	if ((pScanConfig == NULL) || (pCoeffs == NULL) || (patDefCh == NULL))
		return (ERR_DLPSPEC_NULL_POINTER);

	for(j = 0; j < pScanConfig->head.num_sections; j++)
	{
		patDefCh->numPatterns = pScanConfig->chemoSection[j].num_patterns;
		patDefCh->colWidth = pScanConfig->chemoSection[j].width_px;
		for(i = 0; i < pScanConfig->chemoSection[j].num_patterns; i++) {
			ret_val = dlpspec_util_nmToColumn(pScanConfig->chemoSection[j].wavelengths[i],pCoeffs->PixelToWavelengthCoeffs, &mid_x);
			if (ret_val < 0)
			{
				return (ERR_DLPSPEC_INVALID_INPUT);
			}
			patDefCh->colMidPix[i] = mid_x;
			patDefCh->colHeight[i] = pScanConfig->chemoSection[j].heights_px[i];
		}
		ret = dlpspec_scan_chemo_genSinglePattern(patDefCh,pFB,0,j);
		if(ret < 0)
			return ret;
	}

	return pScanConfig->head.num_sections;
}

int32_t dlpspec_scan_chemo_genSinglePattern(const patDefChemo *patDefCh,
		const FrameBufferDescriptor *pFB, uint32_t startPattern, uint32_t currSect)
{
	int i;
	RectangleDescriptor rect;
	int curPattern;
	int patterns_per_image;
	uint32_t curBuffer=0;
	int frameBufferSz = (pFB->width * pFB->height * (pFB->bpp/8));
	FrameBufferDescriptor frameBuffer;

	if ((patDefCh == NULL) || (pFB == NULL))
		return (ERR_DLPSPEC_NULL_POINTER);

	memcpy(&frameBuffer, pFB, sizeof(FrameBufferDescriptor));

	if(frameBuffer.bpp == 16)
		patterns_per_image=16;
	else
		patterns_per_image=24;

	/* Depending on startPattern, skip N buffers */
	curBuffer = startPattern/patterns_per_image;
	frameBuffer.frameBuffer += ((frameBufferSz/4)*curBuffer);
	curPattern = startPattern - curBuffer*patterns_per_image;

	for(i=0; i < patDefCh->numPatterns; i++)
	{
		if(currSect % patterns_per_image == 0)
		{
			//First clear the area of interest
			rect.startX = 0;
			rect.startY = 0;
			rect.height = frameBuffer.height;
			rect.width = frameBuffer.width;
			rect.pixelVal = 0;
			DrawRectangle(&rect, &frameBuffer, true);
		}

		//Guard against rectangles drawn out of the left bound of the frame
		if((patDefCh->colMidPix[i] - patDefCh->colWidth/2) < 0)
			rect.startX = 0;
		else
			rect.startX = patDefCh->colMidPix[i] - patDefCh->colWidth/2;

		rect.startY = 0;
		rect.height = patDefCh->colHeight[i];

		//Guard against rectangles drawn out of the right bound of the frame
		if((rect.startX + patDefCh->colWidth) > pFB->width)
			rect.width = pFB->width - rect.startX;
		else
			rect.width = patDefCh->colWidth;

		rect.pixelVal = 1 << (currSect%patterns_per_image);

		DrawRectangle(&rect, &frameBuffer, false);
		curPattern++;
		if(currSect % patterns_per_image == 0)
		{
			//Advance frame buffer pointer
			frameBuffer.frameBuffer += frameBufferSz/4;
			curBuffer++;
			if(curBuffer == frameBuffer.numFBs)
				break;
		}
	}

	return 1;//(patDefCh->numPatterns);
}

