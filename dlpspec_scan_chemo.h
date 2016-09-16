/*
 * dlpspec_scan_chemo.h
 *
 *  Created on: Sep 6, 2016
 *      Author: Purvesh
 */
#include <stdint.h>
#include "dlpspec_scan.h"
#include "dlpspec_types.h"
#include "dlpspec_helper.h"

#ifndef DLPSPEC_SCAN_CHEMO_H_
#define DLPSPEC_SCAN_CHEMO_H_

typedef struct
{
	uint16_t numPatterns; /**< Number of binary DMD patterns required for the scan */
	uint16_t colHeight[MAX_CHEMO_PATTERNS_PER_SCAN]; /** Height in pixels of each group of on-state pixels in each pattern */
	uint16_t colWidth[MAX_CHEMO_PATTERNS_PER_SCAN]; /**< Width in pixels of each group of on-state pixels in each pattern */
	uint16_t colMidPix[MAX_CHEMO_PATTERNS_PER_SCAN]; /**< The DMD column number corresponding to the center of the on-state group of pixels for each pattern */
} patDefChemo;

// Function prototypes
DLPSPEC_ERR_CODE dlpspec_scan_chemo_interpret(const uScanData *pScanData,scanResults *pResults);
DLPSPEC_ERR_CODE dlpspec_scan_chemo_genPatDef(const chemoScanConfig *pScanConfig,
		const calibCoeffs *pCoeffs, patDefChemo *patDefCh);
int32_t dlpspec_scan_chemo_genPatterns(const patDefChemo *patDefCh,
		const FrameBufferDescriptor *pFB, uint32_t startPattern);

#endif /* DLPSPEC_SCAN_CHEMO_H_ */
