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
	double mid_px_f;
	chemoScanConfig cfg;
	const chemoScanData *pScanData;

	DLPSPEC_ERR_CODE ret_val = (DLPSPEC_PASS);

	if ((pResults == NULL) || (puScanData == NULL))
		return (ERR_DLPSPEC_NULL_POINTER);

	pScanData = &puScanData->chemo_data;

	/*
	    Procedure:
	    1. Compute DC detector level during black patterns (every 25th, as defined
															by sequence, #def.)
	    2. Subtract found DC detector level from remaining measurements.
	    3. Compute wavelength centers for the scanData configuration, using genPatDef
	    4. Copy computed wavelengths, intensities, and header info from scanData
																to scanResults.
	*/

	dlpspec_copy_chemoScanData_hdr_to_scanResults(pScanData, pResults);

	dlpspec_chemo_subtract_remove_dc_level(pScanData, pResults);

	cfg.scan_type = pScanData->chemoCfg.scan_type;
	cfg.scanConfigIndex = pScanData->chemoCfg.scanConfigIndex;
	cfg.num_patterns = pScanData->chemoCfg.num_patterns;
	cfg.num_repeats = pScanData->chemoCfg.num_repeats;
	for(i = 0; i < pScanData->chemoCfg.num_patterns; i++) {
		cfg.widths_px[i] = pScanData->chemoCfg.widths_px[i];
		cfg.heights_px[i] = pScanData->chemoCfg.heights_px[i];
	}
	/* Compute wavelength centers for the scanData configuration, using genPatDef */
	ret_val = dlpspec_scan_chemo_genPatDef(&cfg, &(pScanData->calibration_coeffs),
			&patDef);
	if (ret_val < 0)
	{
		return ret_val;
	}

	for(i=0; i < pResults->length; i++)
	{
		if (pScanData->chemoCfg.widths_px[i] % 2 != 0)
			mid_px_f = patDef.colMidPix[i];
		else
			mid_px_f = patDef.colMidPix[i] - 0.5;

		ret_val = dlpspec_util_columnToNm(mid_px_f,
				&(pScanData->calibration_coeffs.PixelToWavelengthCoeffs[0]),
				&pResults->wavelength[i]);
		if (ret_val < 0)
		{
			return ret_val;
		}
	}

	pResults->pga = pScanData->pga;
	return ret_val;
}

DLPSPEC_ERR_CODE dlpspec_scan_chemo_genPatDef(const chemoScanConfig *pScanConfig, const calibCoeffs *pCoeffs, patDefChemo *patDefCh)
{
	int i;
	double mid_x;

	DLPSPEC_ERR_CODE ret_val = (DLPSPEC_PASS);
	if ((pScanConfig == NULL) || (pCoeffs == NULL) || (patDefCh == NULL))
		return (ERR_DLPSPEC_NULL_POINTER);

	patDefCh->numPatterns = pScanConfig->num_patterns;
	for(i = 0; i < pScanConfig->num_patterns; i++) {
		ret_val = dlpspec_util_nmToColumn(pScanConfig->wavelengths[i],pCoeffs->PixelToWavelengthCoeffs, &mid_x);
		if (ret_val < 0)
		{
			return (ERR_DLPSPEC_INVALID_INPUT);
		}
		patDefCh->colMidPix[i] = mid_x;
		patDefCh->colWidth[i] = pScanConfig->widths_px[i];
		patDefCh->colHeight[i] = pScanConfig->heights_px[i];
	}

	return (DLPSPEC_PASS);
}

int32_t dlpspec_scan_chemo_genPatterns(const patDefChemo *patDefCh,
		const FrameBufferDescriptor *pFB, uint32_t startPattern)
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
		if(curPattern % patterns_per_image == 0)
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
		if((patDefCh->colMidPix[i] - patDefCh->colWidth[i]/2) < 0)
			rect.startX = 0;
		else
			rect.startX = patDefCh->colMidPix[i] - patDefCh->colWidth[i]/2;

		rect.startY = 0;
		rect.height = patDefCh->colHeight[i];

		//Guard against rectangles drawn out of the right bound of the frame
		if((rect.startX + patDefCh->colWidth[i]) > pFB->width)
			rect.width = pFB->width - rect.startX;
		else
			rect.width = patDefCh->colWidth[i];

		rect.pixelVal = 1 << (curPattern%patterns_per_image);

		DrawRectangle(&rect, &frameBuffer, false);
		curPattern++;
		if(curPattern % patterns_per_image == 0)
		{
			//Advance frame buffer pointer
			frameBuffer.frameBuffer += frameBufferSz/4;
			curBuffer++;
			if(curBuffer == frameBuffer.numFBs)
				break;
		}
	}

	return (patDefCh->numPatterns);
}

