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
	int i,j,ind;
	chemoScanConfig cfg;
	const chemoScanData *pScanData;
	int average_intensities[CHEMO_SCAN_MAX_SEQS/2];

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
	cfg.num_seqs = pScanData->chemoCfg.num_seqs;
	cfg.num_repeats = pScanData->chemoCfg.num_repeats;
	cfg.repeat_count = pScanData->chemoCfg.repeat_count;
	cfg.width_px = pScanData->chemoCfg.width_px;
	for(i = 0; i < cfg.num_seqs; i++)
	{
		cfg.seq_id[i] = pScanData->chemoCfg.seq_id[i];
	}

	i=0;
	for(i = 0; i < CHEMO_SCAN_MAX_SEQS; i++)
	{
		average_intensities[i] = 0;
	}

	i = 0;
	ind = 0;
	while(i < pResults->length)
	{
		average_intensities[ind%(cfg.num_seqs/2)] += (pResults->intensity[i] -
												 pResults->intensity[i+1]);
		i+=2;
		ind++;
	}

	pResults->length = (cfg.num_seqs / 2);
	for(i = 0; i < pResults->length; i++)
	{
		pResults->wavelength[i] = 900 + i;
		pResults->intensity[i] = average_intensities[i] / cfg.repeat_count;
	}

	pResults->pga = pScanData->pga;
	return ret_val;
}

chemoSelection dlpspec_get_wavelengths_heights(uint8_t seqId)
{
	int i;
	chemoSelection results;
	uint16_t heights_0[MAX_CHEMO_PATTERNS_PER_SCAN] = {0,1,1,2,2,2,2,2,4,4,5,6,6,7,8,9,10,11,11,12,12,12,12,14,14,14,14,15,15,15,16,17,18,21,22,23,23,24,24,25,25,26,26,27,28,28,29,29,29,30,30,32,32,34,35,36,36,37,37,37,37,39,40,40,40,41,41,42,42,42,42,44,46,48,48,52,52,52,53,55,59,61,63,66,71,75,75,75,76,77,79,80,84,87,93,103};
	uint16_t wavelengths_0[MAX_CHEMO_PATTERNS_PER_SCAN] = {1312,1034,1588,1185,1549,1301,1010,1673,1218,1515,1087,1484,986,1682,1529,1162,1602,1121,1454,907,1365,1279,1701,1413,1640,1317,1668,1649,1469,1686,1202,1081,943,1433,1635,950,1449,1212,1235,1328,1597,1443,1499,1224,1334,1418,1306,1386,1116,1029,1495,919,956,1144,1474,1569,1052,1607,1040,1208,1093,1360,1257,901,1190,1381,937,1654,1268,1127,1535,1544,998,1338,968,1631,1005,1583,1241,1539,1439,1290,1022,1075,1678,1322,1252,1098,1285,1376,1274,1464,1110,961,1295,1104};
	uint16_t heights_1[MAX_CHEMO_PATTERNS_PER_SCAN] = {93,92,68,61,60,59,51,50,47,47,46,41,40,40,38,37,37,36,33,33,32,32,29,28,27,26,25,25,24,23,22,21,20,17,16,15,13,12,12,12,12,11,10,6,5,5,4,3,3,2,2,1,1,0};
	uint16_t wavelengths_1[MAX_CHEMO_PATTERNS_PER_SCAN] = {1659,1593,925,912,1691,1574,1246,1196,1578,1344,992,1070,1428,1696,1046,1167,1489,1504,1371,1626,1173,1402,1524,1016,1519,1156,974,1349,1509,1664,1616,1392,1230,1355,1558,1479,1263,1423,1554,1179,1133,1459,1612,1645,1063,1397,1150,1407,932,1057,1139,981,1564,1621};
	uint16_t heights_2[MAX_CHEMO_PATTERNS_PER_SCAN] = {1,2,2,3,4,4,4,6,6,7,8,8,8,8,8,8,9,9,9,11,11,13,13,14,14,14,15,16,17,17,18,18,20,20,20,21,22,22,23,23,24,27,27,27,27,28,28,28,30,31,31,31,32,35,35,38,39,40,40,41,41,45,45,46,47,49,49,51,52,53,53,55,56,58,59,62,63,64,66,66,66,67,76,76,84};
	uint16_t wavelengths_2[MAX_CHEMO_PATTERNS_PER_SCAN] = {1070,1224,1144,1162,1554,1230,1005,1402,912,1701,1057,956,1208,1549,1150,1338,1344,1212,1597,1696,1673,1016,1098,1252,1578,1649,1040,1574,1691,1371,1544,1010,1196,1022,1659,1443,1469,1484,1127,1268,1274,1075,1418,1218,1167,1602,1558,1621,1524,1179,1246,950,1295,1202,1397,1029,981,992,1034,1241,986,1190,1474,1682,1612,925,1640,1285,907,1386,1046,1052,1423,1509,1489,1156,1110,1081,1263,919,1407,932,1635,1349,1504};
	uint16_t heights_3[MAX_CHEMO_PATTERNS_PER_SCAN] = {100,85,79,77,74,71,70,69,69,67,66,64,64,62,61,57,53,52,51,47,46,45,44,42,40,38,35,34,30,29,27,27,26,25,23,21,21,21,18,18,17,16,16,16,15,15,14,14,13,13,12,11,9,9,6,6,6,5,3,3,3,2,2,1,1};
	uint16_t wavelengths_3[MAX_CHEMO_PATTERNS_PER_SCAN] = {1301,1322,1664,1439,901,1116,1312,1607,1616,1376,998,1433,1454,1104,1654,1499,1392,1290,1087,1529,1334,1668,1257,1686,1569,1365,1306,1449,1495,1564,1279,1133,1464,1515,1459,1328,1593,1678,1185,1360,1093,1626,937,1645,1588,1173,968,1583,1381,1317,943,1355,1535,974,1428,1063,1235,1139,1413,1631,1539,1519,961,1479,1121};

	switch(seqId)
	{
		case 0:
			results.num_patterns = 96;
			for(i = 0; i < results.num_patterns; i++)
				{
					results.wavelengths[i] = wavelengths_0[i];
					results.heights[i] = heights_0[i];
				}
			break;
		case 1:
			results.num_patterns = 54;
			for(i = 0; i < results.num_patterns; i++)
				{
					results.wavelengths[i] = wavelengths_1[i];
					results.heights[i] = heights_1[i];
				}
			break;
		case 2:
			results.num_patterns = 85;
			for(i = 0; i < results.num_patterns; i++)
				{
					results.wavelengths[i] = wavelengths_2[i];
					results.heights[i] = heights_2[i];
				}
			break;
		case 3:
			results.num_patterns = 65;
			for(i = 0; i < results.num_patterns; i++)
				{
					results.wavelengths[i] = wavelengths_3[i];
					results.heights[i] = heights_3[i];
				}
			break;
	}

	return results;
}

DLPSPEC_ERR_CODE dlpspec_scan_chemo_genPatDef(const chemoScanConfig *pScanConfig, const calibCoeffs *pCoeffs, patDefChemo *patDefCh,uint32_t seqId)
{
	int i;
	double mid_x;

	DLPSPEC_ERR_CODE ret_val = (DLPSPEC_PASS);
	if ((pScanConfig == NULL) || (pCoeffs == NULL) || (patDefCh == NULL))
		return (ERR_DLPSPEC_NULL_POINTER);

	chemoSelection selection = dlpspec_get_wavelengths_heights(seqId);
	patDefCh->numPatterns = selection.num_patterns;
	patDefCh->colWidth = pScanConfig->width_px;

	for(i = 0; i < patDefCh->numPatterns; i++) {
		ret_val = dlpspec_util_nmToColumn(selection.wavelengths[i],pCoeffs->PixelToWavelengthCoeffs, &mid_x);
		if (ret_val < 0)
		{
			return (ERR_DLPSPEC_INVALID_INPUT);
		}
		patDefCh->colMidPix[i] = mid_x;
		patDefCh->colHeight[i] = selection.heights[i];
	}

	return (DLPSPEC_PASS);
}

int32_t dlpspec_scan_chemo_genSinglePatterns(const patDefChemo *patDefCh,
		const FrameBufferDescriptor *pFB, uint32_t seqId)
{
	int i;
	RectangleDescriptor rect;
	int curSeqId;
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
	curBuffer = seqId/patterns_per_image;
	frameBuffer.frameBuffer += ((frameBufferSz/4)*curBuffer);
	curSeqId = seqId - curBuffer*patterns_per_image;

	if(curBuffer == frameBuffer.numFBs)
		return 0;

	if(curSeqId % patterns_per_image == 0)
	{
		//First clear the area of interest
		rect.startX = 0;
		rect.startY = 0;
		rect.height = frameBuffer.height;
		rect.width = frameBuffer.width;
		rect.pixelVal = 0;
		DrawRectangle(&rect, &frameBuffer, true);
	}

	for(i=0; i < patDefCh->numPatterns; i++)
	{
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

		rect.pixelVal = 1 << (curSeqId%patterns_per_image);

		DrawRectangle(&rect, &frameBuffer, false);
	}

	return 1;
}

