/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "globals.h"
#include "tightString.h"
#include "readSet.h"
#include "utility.h"
#include "binarySequences.h"

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#include "../third-party/zlib-1.2.3/Win32/include/zlib.h"
#else
#include "../third-party/zlib-1.2.3/zlib.h"
#endif

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

typedef struct {
	Category    m_numCategories;
	uint32_t	m_magic;
	boolean		m_bColor;
	uint64_t	m_sequenceCnt;
	uint64_t	m_timeStamp;
	uint64_t	m_seqNuclStoreSize;
	uint64_t	m_minSeqLen;
	uint64_t	m_maxSeqLen;
	uint64_t	m_totalSeqLen;
	boolean		m_bFileWriteCompleted;
} CnyUnifiedSeqFileHeader;

typedef struct {
	Category    m_numCategories;
	Category    m_currCategory;
	uint64_t	m_minSeqLen;
	uint64_t	m_maxSeqLen;
	uint8_t *			m_pReadBuffer;
	uint8_t *			m_pReadBufEnd;
	uint64_t			m_readBufPos;
	uint8_t *			m_pCurrentReadPtr;
	uint8_t *			m_pNextReadPtr;
	uint64_t			m_currentNuclReadIdx;
	uint64_t			m_currentReadLength;
	int32_t         m_referenceID;
	int32_t         m_pos;
	boolean         m_bIsRef;
} CnyUnifiedSeqReadInfo;

#define USF_READ_BUF_SIZE				(64*1024)

// write defines, typedefs, and protos
#define WRITE_BUF_SHFT					16							// byte shift and mask
#define WRITE_BUF_SIZE					(1<<WRITE_BUF_SHFT)
#define WRITE_BUF_MASK					(WRITE_BUF_SIZE-1)
#define SHORT_NUCL_LENGTH				128							// Nucleotide length (2-bits each)

struct binarySequencesWriter_st {
        FILE *		m_pFile;
        CnyUnifiedSeqFileHeader m_unifiedSeqFileHeader;
	uint64_t	m_insertStartIndex;
	uint64_t	m_insertLength;
	uint64_t	m_insertLengthIndex;
	uint64_t	m_insertCurrentIndex;
	uint32_t	m_hostBuffersInUse;
	uint32_t	m_fileSegmentWriteIdx;
	uint8_t	*	m_pWriteBuffer[3];
	uint8_t *	m_pHostBufPtr;
	uint8_t *	m_pHostLengthBufPtr;
	uint8_t *	m_pHostLengthBufPtrMax;
	uint8_t *	m_pHostBufPtrMax;
	int64_t		m_hostBufferFilePos[3];
	int32_t         m_referenceID;
	int32_t         m_pos;
	boolean         m_bIsRef;
};

FILE *openCnySeqForRead(const char *fileName, CnyUnifiedSeqFileHeader *seqFileHeader)
{
	FILE *pFile;
	if ((pFile = fopen(fileName, "rb")) == 0) {
		velvetLog("Unable to open %s for reading\n", fileName);
		return NULL;
	}

	if (fread(seqFileHeader, sizeof(*seqFileHeader), 1, pFile) != 1) {
		velvetLog("Unable to read file %s\n", fileName);
		fclose(pFile);
		return NULL;
	}

	if (strncmp((char *)&seqFileHeader->m_magic, "CSQ0", 4) != 0) {
		velvetLog("Unknown format for file %s\n", fileName);
		fclose(pFile);
		return NULL;
	}

	if (seqFileHeader->m_bFileWriteCompleted == false) {
		velvetLog("Corrupted file, %s\n", fileName);
		fclose(pFile);
		return NULL;
	}

	if (seqFileHeader->m_numCategories > CATEGORIES) {
		velvetLog("File %s has %d categories, please rebuild velvet to match\n", fileName, seqFileHeader->m_numCategories);
		fclose(pFile);
		return NULL;
	}

#ifdef COLOR
	if (!seqFileHeader->m_bColor) {
		velvetLog("File %s does not specify color, please rebuild velvet to match\n", fileName);
		fclose(pFile);
		return NULL;
	}
#else
	if (seqFileHeader->m_bColor) {
		velvetLog("File %s specifies color, please rebuild velvet to match\n", fileName);
		fclose(pFile);
		return NULL;
	}
#endif
	return pFile;
}

static boolean refillCnySeqReadBuffer(FILE *pFile, CnyUnifiedSeqFileHeader *seqFileHeader, CnyUnifiedSeqReadInfo *pReadInfo) 
{
	uint64_t readLen = (USF_READ_BUF_SIZE < seqFileHeader->m_seqNuclStoreSize - pReadInfo->m_readBufPos) ?
		USF_READ_BUF_SIZE : seqFileHeader->m_seqNuclStoreSize - pReadInfo->m_readBufPos;

	if (readLen == 0)
		return false;

	if (fread(pReadInfo->m_pReadBuffer, (uint32_t)readLen, 1, pFile) != 1) {
		velvetLog("Unable to read file\n");
		exit (1);
	}

	pReadInfo->m_pCurrentReadPtr = pReadInfo->m_pReadBuffer;
	pReadInfo->m_pReadBufEnd = pReadInfo->m_pReadBuffer + readLen;
	pReadInfo->m_readBufPos += readLen;
	if (pReadInfo->m_pNextReadPtr >= pReadInfo->m_pReadBufEnd) {
		pReadInfo->m_pNextReadPtr -= USF_READ_BUF_SIZE;
	}

	return true;
}

static int32_t readCnySeqUint8(FILE *pFile, CnyUnifiedSeqFileHeader *pHeader,
		CnyUnifiedSeqReadInfo *pReadInfo)
{
	if (pReadInfo->m_pCurrentReadPtr == pReadInfo->m_pReadBufEnd && !refillCnySeqReadBuffer(pFile, pHeader, pReadInfo))
	{
		return -1;
	}

	return *pReadInfo->m_pCurrentReadPtr++;
}

static uint32_t readCnySeqUint32(FILE *pFile, CnyUnifiedSeqFileHeader *pHeader, CnyUnifiedSeqReadInfo *pReadInfo) 
{
	uint32_t data;
	data = 0;
	int i;
	for (i = 0; i < 4; i += 1)
		data |= readCnySeqUint8(pFile, pHeader, pReadInfo) << (i*8);
	return data;
}

static boolean advanceCnySeqCurrentRead(FILE *pFile, CnyUnifiedSeqFileHeader *pHeader, CnyUnifiedSeqReadInfo *pReadInfo)
{
	// Perform consistency check, unused bits of previous sequence should have a fixed pattern
	uint8_t finalNuclOffset = 1;
	if (pReadInfo->m_bIsRef) {
		finalNuclOffset += (sizeof(pReadInfo->m_referenceID) + sizeof(pReadInfo->m_pos));
	}
	if ((pReadInfo->m_currentReadLength & 3) != 0 && ((pReadInfo->m_pNextReadPtr - finalNuclOffset) >= pReadInfo->m_pReadBuffer)) {
		uint8_t mask = 0xFF << (pReadInfo->m_currentReadLength & 3) * 2;
		if ((*(pReadInfo->m_pNextReadPtr - finalNuclOffset) & mask) != (0xAA & mask)) {
			velvetLog("Cny seq consistency check failed in advance\n");
			exit(1);
		}
	}

	pReadInfo->m_pCurrentReadPtr = pReadInfo->m_pNextReadPtr;
	pReadInfo->m_currentNuclReadIdx = 0;

	// clear ref flag before each code check
	pReadInfo->m_bIsRef = false;

	for(;;) {
		int32_t	code = readCnySeqUint8(pFile, pHeader, pReadInfo);
		// printf("checking code %d\n", code);
		switch (code & 0xc0) {
			case 0x00:	// short sequence
			case 0x40:
				pReadInfo->m_currentReadLength = code & 0x7f;
				// printf("short len %d\n", (int32_t) pReadInfo->m_currentReadLength);
				pReadInfo->m_pNextReadPtr = pReadInfo->m_pCurrentReadPtr + ((pReadInfo->m_currentReadLength + 3) >> 2);
				break;
			case 0x80:	// long sequence
				pReadInfo->m_currentReadLength = readCnySeqUint32(pFile, pHeader, pReadInfo);
				// printf("long len %d\n", (int32_t) pReadInfo->m_currentReadLength);
				pReadInfo->m_pNextReadPtr = pReadInfo->m_pCurrentReadPtr + ((pReadInfo->m_currentReadLength + 3) >> 2);
				if (code & 0x20) {
				    // ref info present
				    pReadInfo->m_bIsRef = true;
				    pReadInfo->m_pNextReadPtr += (sizeof(pReadInfo->m_referenceID) + sizeof(pReadInfo->m_pos));
				}
				break;
			case 0xc0:	// new file / category
				if (code == EOF) {
					return false;
				}
				pReadInfo->m_currCategory = (Category) readCnySeqUint32(pFile, pHeader, pReadInfo);
				if (pReadInfo->m_currCategory < 0 || pReadInfo->m_currCategory > REFERENCE) {
					velvetLog("Illegal category %d\n", (int32_t) pReadInfo->m_currCategory);
					exit(1);
				}
				continue;
		}

		if (pReadInfo->m_currentReadLength > pReadInfo->m_maxSeqLen ||
				pReadInfo->m_currentReadLength < pReadInfo->m_minSeqLen) {
			velvetLog("Cny seq consistency check failed, len mismatch\n");
			exit(1);
		}

		return true;

	}
}

static void resetCnySeqCurrentRead(FILE *pFile, CnyUnifiedSeqFileHeader *pHeader, CnyUnifiedSeqReadInfo *pReadInfo) 
{
	pReadInfo->m_pReadBufEnd = pReadInfo->m_pReadBuffer;
	pReadInfo->m_pNextReadPtr = pReadInfo->m_pReadBuffer;
	pReadInfo->m_pCurrentReadPtr = pReadInfo->m_pReadBuffer;
	pReadInfo->m_currentReadLength = 0;
	pReadInfo->m_readBufPos = 0;

	if (fseek(pFile, sizeof(CnyUnifiedSeqFileHeader), SEEK_SET) < 0) {
		perror("Unable to seek\n");
		exit(1);
	}

	advanceCnySeqCurrentRead(pFile, pHeader, pReadInfo);
}

static void getCnySeqNucl(FILE *pFile, CnyUnifiedSeqFileHeader *pHeader,
		CnyUnifiedSeqReadInfo *pReadInfo, uint8_t *sequence) {
	uint32_t nuclIdx;
	for (nuclIdx = 0; nuclIdx < pReadInfo->m_currentReadLength; nuclIdx += 1) {
		uint8_t nucleotide;
		uint8_t byte;
		if ((nuclIdx % 4) == 0) {
			byte = (uint8_t)readCnySeqUint8(pFile, pHeader, pReadInfo);
		}
		nucleotide = (byte >> ((nuclIdx & 3) * 2)) & 0x03;
		sequence[nuclIdx] = "ACGT"[nucleotide];
	}
}

ReadSet *importCnyReadSet(char *filename)
{
	FILE *pFile;
	unsigned char *sequence = NULL;
	IDnum sequenceCount, sequenceIndex;
	ReadSet *reads;
	CnyUnifiedSeqFileHeader header;
	CnyUnifiedSeqReadInfo readInfo;
	memset(&readInfo, 0, sizeof(readInfo));

	pFile = openCnySeqForRead(filename, &header);
	readInfo.m_numCategories = header.m_numCategories;
	readInfo.m_minSeqLen = header.m_minSeqLen;
	readInfo.m_maxSeqLen = header.m_maxSeqLen;
	readInfo.m_bIsRef = false;

	if (pFile != NULL)
		velvetLog("Reading CNY read set file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	readInfo.m_pReadBuffer = mallocOrExit(USF_READ_BUF_SIZE, sizeof(*readInfo.m_pReadBuffer));
	readInfo.m_pCurrentReadPtr = readInfo.m_pReadBufEnd = 0;

	reads = newReadSet();

	resetCnySeqCurrentRead(pFile, &header, &readInfo);
	sequenceCount = header.m_sequenceCnt;

	velvetLog("%li sequences found\n", (long) sequenceCount);

	reads->readCount = sequenceCount;

	if (reads->readCount == 0) {
		reads->sequences = NULL;
		reads->categories = NULL;
	        fclose(pFile);
		return reads;	
	}

	reads->sequences = callocOrExit(sequenceCount, char *);
	reads->categories = callocOrExit(sequenceCount, Category);

	// read all sequence and category info in one pass
	for (sequenceIndex = 0; sequenceIndex < sequenceCount; sequenceIndex += 1) {
		reads->categories[sequenceIndex] = readInfo.m_currCategory;
		if (sizeof(ShortLength) == sizeof(int16_t) && readInfo.m_currentReadLength > SHRT_MAX) {
			velvetLog("Read %li of length %lli, longer than limit %i\n",
					(long) sequenceIndex + 1, (long long) readInfo.m_currentReadLength, SHRT_MAX);
			velvetLog("You should modify ShortLength to int32_t in globals.h (around) line 84, recompile and re-run\n");
			exit(1);
		}
		reads->sequences[sequenceIndex] = mallocOrExit(readInfo.m_currentReadLength + 1, char);
		sequence = (unsigned char *) reads->sequences[sequenceIndex];
		getCnySeqNucl(pFile, &header, &readInfo, sequence);
		if (readInfo.m_bIsRef) {
		    readInfo.m_referenceID = readCnySeqUint32(pFile, &header, &readInfo);
		    readInfo.m_pos = readCnySeqUint32(pFile, &header, &readInfo);
		}
		sequence[readInfo.m_currentReadLength] = '\0';
		if (sequenceIndex < sequenceCount) {
			advanceCnySeqCurrentRead(pFile, &header, &readInfo);
		}
	}

	fclose(pFile);

	velvetLog("Done\n");
	return reads;

}

// write routines
#define ADENINE		0
#define CYTOSINE	1
#define GUANINE		2
#define THYMINE		3
#define INVALID         5

static void cnySeqHostBufferFull(BinarySequencesWriter *cnySeqWriteInfo)
{
    // The current Host buffer is full
    switch (cnySeqWriteInfo->m_hostBuffersInUse) {
	case 1:	// buf[0] is full, start using buf[1]
	    cnySeqWriteInfo->m_pHostBufPtr = cnySeqWriteInfo->m_pWriteBuffer[1];
	    cnySeqWriteInfo->m_pHostBufPtrMax = cnySeqWriteInfo->m_pHostBufPtr + WRITE_BUF_SIZE;
	    cnySeqWriteInfo->m_hostBufferFilePos[1] = cnySeqWriteInfo->m_hostBufferFilePos[0] + WRITE_BUF_SIZE;
	    cnySeqWriteInfo->m_hostBuffersInUse = 2;
	    break;
	case 2: // buf[0] and buf[1] are full, start using buf[2]
	    cnySeqWriteInfo->m_pHostBufPtr = cnySeqWriteInfo->m_pWriteBuffer[2];
	    cnySeqWriteInfo->m_pHostBufPtrMax = cnySeqWriteInfo->m_pHostBufPtr + WRITE_BUF_SIZE;
	    cnySeqWriteInfo->m_hostBufferFilePos[2] = cnySeqWriteInfo->m_hostBufferFilePos[1] + WRITE_BUF_SIZE;
	    cnySeqWriteInfo->m_hostBuffersInUse = 3;
	    break;
	case 3: // all three buffers are full, write out buf[2] and reuse
	    if (fseek(cnySeqWriteInfo->m_pFile, cnySeqWriteInfo->m_hostBufferFilePos[2], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }

	    if (fwrite(cnySeqWriteInfo->m_pWriteBuffer[2], WRITE_BUF_SIZE, 1, cnySeqWriteInfo->m_pFile) != 1) {
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }

	    // also need to copy to coprocessor memory

	    cnySeqWriteInfo->m_pHostBufPtr = cnySeqWriteInfo->m_pWriteBuffer[2];
	    cnySeqWriteInfo->m_pHostBufPtrMax = cnySeqWriteInfo->m_pHostBufPtr + WRITE_BUF_SIZE;
	    cnySeqWriteInfo->m_hostBufferFilePos[2] = cnySeqWriteInfo->m_hostBufferFilePos[2] + WRITE_BUF_SIZE;
	    break;
	default:
	    velvetLog("Unknown CnySeq host buffer state %d\n", cnySeqWriteInfo->m_hostBuffersInUse);
	    exit(1);
	    break;
    }
}

static void moveCnySeqNucleotides(BinarySequencesWriter *cnySeqWriteInfo)
{
    // move nucleotides in buffer to allow a four byte length value
    // the current sequence may span two buffers

    uint64_t bufIdx = (cnySeqWriteInfo->m_hostBuffersInUse == 2) ? (cnySeqWriteInfo->m_pHostBufPtr - cnySeqWriteInfo->m_pWriteBuffer[1] + WRITE_BUF_SIZE) : (cnySeqWriteInfo->m_pHostBufPtr - cnySeqWriteInfo->m_pWriteBuffer[0]);
    if (bufIdx + 4 > 2 * WRITE_BUF_SIZE) {
	velvetLog("CnySeq bufIdx %ld too large\n", bufIdx);
	exit(1);
    }

    if (bufIdx + 4 >= WRITE_BUF_SIZE) {
	// continue writing to buf[1]
	cnySeqWriteInfo->m_pHostBufPtr = cnySeqWriteInfo->m_pWriteBuffer[1] + (bufIdx + 4 - WRITE_BUF_SIZE);
	cnySeqWriteInfo->m_pHostBufPtrMax = cnySeqWriteInfo->m_pWriteBuffer[1] + WRITE_BUF_SIZE;
	cnySeqWriteInfo->m_hostBufferFilePos[1] = cnySeqWriteInfo->m_hostBufferFilePos[0] + WRITE_BUF_SIZE;
	cnySeqWriteInfo->m_hostBuffersInUse = 2;
    } else
	cnySeqWriteInfo->m_pHostBufPtr += 4;

    cnySeqWriteInfo->m_insertCurrentIndex += 16;

    uint64_t cnt;
    for (cnt = (cnySeqWriteInfo->m_insertLength+3)>>2; cnt > 0; cnt -= 1) {
	cnySeqWriteInfo->m_pWriteBuffer[(bufIdx+4) >> WRITE_BUF_SHFT][(bufIdx+4) & WRITE_BUF_MASK] = cnySeqWriteInfo->m_pWriteBuffer[bufIdx >> WRITE_BUF_SHFT][bufIdx & WRITE_BUF_MASK];
	bufIdx -= 1;
    }
}

static void writeCnySeqNucleotide(uint8_t nucleotide, BinarySequencesWriter *cnySeqWriteInfo)
{
    if (cnySeqWriteInfo->m_insertLength == SHORT_NUCL_LENGTH-1) {
	moveCnySeqNucleotides(cnySeqWriteInfo);
    }
    if ((cnySeqWriteInfo->m_insertCurrentIndex & 0x3) == 0)
	*cnySeqWriteInfo->m_pHostBufPtr = 0;

    *cnySeqWriteInfo->m_pHostBufPtr = *cnySeqWriteInfo->m_pHostBufPtr | (nucleotide << ((cnySeqWriteInfo->m_insertCurrentIndex & 0x3) * 2));

    cnySeqWriteInfo->m_insertLength += 1;
    cnySeqWriteInfo->m_insertCurrentIndex += 1;

    if ((cnySeqWriteInfo->m_insertCurrentIndex & 0x3) == 0) {
	cnySeqWriteInfo->m_pHostBufPtr += 1;

	if (cnySeqWriteInfo->m_pHostBufPtr == cnySeqWriteInfo->m_pHostBufPtrMax)
	    cnySeqHostBufferFull(cnySeqWriteInfo);
    }
}

void cnySeqInsertNucleotideString(const char *pReadBuf, BinarySequencesWriter *cnySeqWriteInfo) {
    uint8_t		nucleotide;

    static boolean bInit = false;
    static uint8_t charMap[256];

    if (!bInit) {
	bInit = true;
	memset(charMap, INVALID, 256);
	charMap[(int)'A'] = charMap[(int)'a'] = ADENINE;
	charMap[(int)'C'] = charMap[(int)'c'] = CYTOSINE;
	charMap[(int)'G'] = charMap[(int)'g'] = GUANINE;
	charMap[(int)'T'] = charMap[(int)'t'] = THYMINE;
	charMap[(int)'\0'] = 4;
    }

    for (;;) {
	nucleotide = charMap[(int)*pReadBuf];
	if (nucleotide < 4) {
	    writeCnySeqNucleotide(nucleotide, cnySeqWriteInfo);
	    pReadBuf += 1;
	    continue;
	} else if (nucleotide == 4) {
	    return;
	} else {
	    velvetLog("CnySeq unexpected char %d\n", (int) *pReadBuf);
	    exit(1);
	}
    }
}

BinarySequencesWriter * openCnySeqForWrite(const char *unifiedSeqFileName)
{
    BinarySequencesWriter *cnySeqWriteInfo = callocOrExit(1, BinarySequencesWriter);
    cnySeqWriteInfo->m_pWriteBuffer[0] = NULL;
    cnySeqWriteInfo->m_pWriteBuffer[1] = NULL;
    cnySeqWriteInfo->m_pWriteBuffer[2] = NULL;

#ifdef COLOR
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_bColor = true;
#else
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_bColor = false;
#endif

    if ((cnySeqWriteInfo->m_pFile = fopen(unifiedSeqFileName, "wb")) == 0) {
	velvetLog("Unable to open %s for writing\n", unifiedSeqFileName);
	exit(1);
    }

    memcpy(&cnySeqWriteInfo->m_unifiedSeqFileHeader.m_magic, "CSQ0", 4);
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_timeStamp = time(0);
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_bFileWriteCompleted = false;

    if (fwrite(&cnySeqWriteInfo->m_unifiedSeqFileHeader, sizeof(CnyUnifiedSeqFileHeader), 1, cnySeqWriteInfo->m_pFile) != 1) {
	velvetLog("Unable to write file %s\n", unifiedSeqFileName);
	exit(1);
    }

    cnySeqWriteInfo->m_insertCurrentIndex = 0;
    cnySeqWriteInfo->m_pWriteBuffer[0] = mallocOrExit(WRITE_BUF_SIZE, uint8_t);
    cnySeqWriteInfo->m_pWriteBuffer[1] = mallocOrExit(WRITE_BUF_SIZE, uint8_t);
    cnySeqWriteInfo->m_pWriteBuffer[2] = mallocOrExit(WRITE_BUF_SIZE, uint8_t);

    cnySeqWriteInfo->m_hostBufferFilePos[0] = sizeof(CnyUnifiedSeqFileHeader);

    cnySeqWriteInfo->m_pHostBufPtr = cnySeqWriteInfo->m_pWriteBuffer[0];
    cnySeqWriteInfo->m_pHostBufPtrMax = cnySeqWriteInfo->m_pWriteBuffer[0] + WRITE_BUF_SIZE;
    cnySeqWriteInfo->m_hostBuffersInUse = 1;
    cnySeqWriteInfo->m_fileSegmentWriteIdx = 0;	// file segment currently being written
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_sequenceCnt = 0;
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_minSeqLen = ~0LL;
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_maxSeqLen = 0;
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_totalSeqLen = 0;
    return cnySeqWriteInfo;
}

static void alignCnySeqToNextByteBoundary(BinarySequencesWriter *cnySeqWriteInfo)
{
    if ((cnySeqWriteInfo->m_insertCurrentIndex & 0x3) != 0)
	cnySeqWriteInfo->m_pHostBufPtr += 1;

    cnySeqWriteInfo->m_insertCurrentIndex = (cnySeqWriteInfo->m_insertCurrentIndex + 3) & ~0x3LL;

    if (cnySeqWriteInfo->m_pHostBufPtr == cnySeqWriteInfo->m_pHostBufPtrMax) {
	cnySeqHostBufferFull(cnySeqWriteInfo);
    }
}

static void writeCnySeqUint8(uint8_t uint8, BinarySequencesWriter *cnySeqWriteInfo)
{
    *cnySeqWriteInfo->m_pHostBufPtr++ = uint8;
    cnySeqWriteInfo->m_insertCurrentIndex += 4;

    if (cnySeqWriteInfo->m_pHostBufPtr == cnySeqWriteInfo->m_pHostBufPtrMax) {
	cnySeqHostBufferFull(cnySeqWriteInfo);
    }
}

static void writeCnySeqUint32(uint32_t uint32, BinarySequencesWriter *cnySeqWriteInfo)
{
    int i;
    for (i = 0; i < 4; i += 1)
	writeCnySeqUint8((uint8_t)(uint32 >> (i*8)), cnySeqWriteInfo);
}

void inputCnySeqFileStart(Category category, BinarySequencesWriter *cnySeqWriteInfo)
{
    if (category > REFERENCE) {
	velvetLog("Found category %d beyond max of %d\n", category, REFERENCE);
	exit(1);
    }

    alignCnySeqToNextByteBoundary(cnySeqWriteInfo);
    writeCnySeqUint8(0xc0, cnySeqWriteInfo);
    writeCnySeqUint32(category, cnySeqWriteInfo);
}

void cnySeqInsertStart(BinarySequencesWriter *cnySeqWriteInfo)
{
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_sequenceCnt += 1;

    alignCnySeqToNextByteBoundary(cnySeqWriteInfo);

    cnySeqWriteInfo->m_insertLength = 0;
    cnySeqWriteInfo->m_pHostLengthBufPtr = cnySeqWriteInfo->m_pHostBufPtr;
    cnySeqWriteInfo->m_pHostLengthBufPtrMax = cnySeqWriteInfo->m_pHostBufPtrMax;
    cnySeqWriteInfo->m_pHostBufPtr += 1;
    cnySeqWriteInfo->m_insertLengthIndex = cnySeqWriteInfo->m_insertCurrentIndex >> 2;	// byte index
    cnySeqWriteInfo->m_insertCurrentIndex += 4; // allow for single byte header
    cnySeqWriteInfo->m_insertStartIndex = cnySeqWriteInfo->m_insertCurrentIndex;

    if (cnySeqWriteInfo->m_pHostBufPtr == cnySeqWriteInfo->m_pHostBufPtrMax)
    {
	cnySeqHostBufferFull(cnySeqWriteInfo);
    }

}

void cnySeqInsertEnd(BinarySequencesWriter *cnySeqWriteInfo)
{
    uint8_t *tmp;

    if (cnySeqWriteInfo->m_bIsRef && cnySeqWriteInfo->m_insertLength < SHORT_NUCL_LENGTH) {
	moveCnySeqNucleotides(cnySeqWriteInfo);
    }

    // fill last few empty nucleotides with a fixed pattern for consistency checking
    if ((cnySeqWriteInfo->m_insertCurrentIndex & 0x3) != 0) {
	*cnySeqWriteInfo->m_pHostBufPtr |= 0xAA << ((cnySeqWriteInfo->m_insertCurrentIndex & 0x3)*2);
    }

    // collect read length statistics
    if (cnySeqWriteInfo->m_unifiedSeqFileHeader.m_minSeqLen > cnySeqWriteInfo->m_insertLength)
	cnySeqWriteInfo->m_unifiedSeqFileHeader.m_minSeqLen = cnySeqWriteInfo->m_insertLength;
    if (cnySeqWriteInfo->m_unifiedSeqFileHeader.m_maxSeqLen < cnySeqWriteInfo->m_insertLength)
	cnySeqWriteInfo->m_unifiedSeqFileHeader.m_maxSeqLen = cnySeqWriteInfo->m_insertLength;

    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_totalSeqLen += cnySeqWriteInfo->m_insertLength;

    if (cnySeqWriteInfo->m_insertLength >= SHORT_NUCL_LENGTH || cnySeqWriteInfo->m_bIsRef) {
	if (cnySeqWriteInfo->m_bIsRef) {
	    // set ref bit
	    *(cnySeqWriteInfo->m_pHostLengthBufPtr++) = 0xa0 | ((cnySeqWriteInfo->m_insertLength >> 32) & 0x1f);
	} else {
	    // one byte control, four byte length
	    *(cnySeqWriteInfo->m_pHostLengthBufPtr++) = 0x80 | ((cnySeqWriteInfo->m_insertLength >> 32) & 0x1f);
	}
	int idx;
	for (idx = 0; idx < 4; idx += 1) {
	    if (cnySeqWriteInfo->m_pHostLengthBufPtr == cnySeqWriteInfo->m_pHostLengthBufPtrMax) {
		if (cnySeqWriteInfo->m_hostBuffersInUse < 2) {
		    velvetLog("CnySeq m_hostBuffersInUse %d\n", cnySeqWriteInfo->m_hostBuffersInUse);
		    exit(1);
		}
		cnySeqWriteInfo->m_pHostLengthBufPtr = cnySeqWriteInfo->m_pWriteBuffer[1];
	    }
	    *cnySeqWriteInfo->m_pHostLengthBufPtr++ = (cnySeqWriteInfo->m_insertLength >> (idx*8)) & 0xff;
	}
	if (cnySeqWriteInfo->m_bIsRef) {
	    // write out map info

	    alignCnySeqToNextByteBoundary(cnySeqWriteInfo);

	    int idx;
	    for (idx = 0; idx < 4; idx += 1) {
		if (cnySeqWriteInfo->m_pHostBufPtr == cnySeqWriteInfo->m_pHostBufPtrMax)
		    cnySeqHostBufferFull(cnySeqWriteInfo);
		*cnySeqWriteInfo->m_pHostBufPtr++ = (cnySeqWriteInfo->m_referenceID >> (idx*8)) & 0xff;
		cnySeqWriteInfo->m_insertCurrentIndex += 4; // single byte
	    }
	    for (idx = 0; idx < 4; idx += 1) {
		if (cnySeqWriteInfo->m_pHostBufPtr == cnySeqWriteInfo->m_pHostBufPtrMax)
		    cnySeqHostBufferFull(cnySeqWriteInfo);
		*cnySeqWriteInfo->m_pHostBufPtr++ = (cnySeqWriteInfo->m_pos >> (idx*8)) & 0xff;
		cnySeqWriteInfo->m_insertCurrentIndex += 4; // single byte
	    }
	    cnySeqWriteInfo->m_bIsRef = false;

	}

    } else {
	// one byte length;
	*cnySeqWriteInfo->m_pHostLengthBufPtr = (uint8_t)cnySeqWriteInfo->m_insertLength;
    }

    alignCnySeqToNextByteBoundary(cnySeqWriteInfo);

    switch (cnySeqWriteInfo->m_hostBuffersInUse) {
	case 1:	// buf[0] is being written
	    break;
	case 2: // buf[0] and buf[1] are being written, write buf[0] to disk
	    if (fseek(cnySeqWriteInfo->m_pFile, cnySeqWriteInfo->m_hostBufferFilePos[0], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }

	    if (fwrite(cnySeqWriteInfo->m_pWriteBuffer[0], WRITE_BUF_SIZE, 1, cnySeqWriteInfo->m_pFile) != 1) {
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }

	    // swap buf[0] and buf[1]
	    tmp = cnySeqWriteInfo->m_pWriteBuffer[0];
	    cnySeqWriteInfo->m_pWriteBuffer[0] = cnySeqWriteInfo->m_pWriteBuffer[1];
	    cnySeqWriteInfo->m_pWriteBuffer[1] = tmp;

	    cnySeqWriteInfo->m_hostBufferFilePos[0] = cnySeqWriteInfo->m_hostBufferFilePos[1];

	    cnySeqWriteInfo->m_hostBuffersInUse = 1;
	    break;
	case 3: // buf[0], [1] and [2] are in use, write buf[0] and [1] to disk
	    if (fseek(cnySeqWriteInfo->m_pFile, cnySeqWriteInfo->m_hostBufferFilePos[0], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }


	    if (fwrite(cnySeqWriteInfo->m_pWriteBuffer[0], WRITE_BUF_SIZE, 1, cnySeqWriteInfo->m_pFile) != 1){
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }


	    if (fseek(cnySeqWriteInfo->m_pFile, cnySeqWriteInfo->m_hostBufferFilePos[1], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }


	    if (fwrite(cnySeqWriteInfo->m_pWriteBuffer[1], WRITE_BUF_SIZE, 1, cnySeqWriteInfo->m_pFile) != 1){
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }

	    // swap buf[0] and buf[2]
	    tmp = cnySeqWriteInfo->m_pWriteBuffer[0];
	    cnySeqWriteInfo->m_pWriteBuffer[0] = cnySeqWriteInfo->m_pWriteBuffer[2];
	    cnySeqWriteInfo->m_pWriteBuffer[2] = tmp;

	    cnySeqWriteInfo->m_hostBufferFilePos[0] = cnySeqWriteInfo->m_hostBufferFilePos[2];

	    cnySeqWriteInfo->m_hostBuffersInUse = 1;
	    break;
    }
}

void closeCnySeqForWrite(BinarySequencesWriter *cnySeqWriteInfo)
{
    if (cnySeqWriteInfo->m_pWriteBuffer[0])
	free(cnySeqWriteInfo->m_pWriteBuffer[0]);
    if (cnySeqWriteInfo->m_pWriteBuffer[1])
	free(cnySeqWriteInfo->m_pWriteBuffer[1]);
    if (cnySeqWriteInfo->m_pWriteBuffer[2])
	free(cnySeqWriteInfo->m_pWriteBuffer[2]);

    // should be only one buffer in use
    if (cnySeqWriteInfo->m_hostBuffersInUse != 1) {
	velvetLog("CnySeq host buffers in use %d\n", cnySeqWriteInfo->m_hostBuffersInUse);
	exit(1);
    }

    if (fseek(cnySeqWriteInfo->m_pFile, cnySeqWriteInfo->m_hostBufferFilePos[0], SEEK_SET) < 0) {
	velvetLog("Unable to seek CnySeq\n");
	exit(1);
    }

    if (fwrite(cnySeqWriteInfo->m_pWriteBuffer[0], (uint32_t)(cnySeqWriteInfo->m_pHostBufPtr - cnySeqWriteInfo->m_pWriteBuffer[0]), 1, cnySeqWriteInfo->m_pFile) != 1) {
	velvetLog("Unable to write CnySeq\n");
	exit(1);
    }

    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_bFileWriteCompleted = true;
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_seqNuclStoreSize = cnySeqWriteInfo->m_insertCurrentIndex >> 2;
    cnySeqWriteInfo->m_unifiedSeqFileHeader.m_numCategories = CATEGORIES;

    if (fseek(cnySeqWriteInfo->m_pFile, 0, SEEK_SET) < 0) {
	velvetLog("Unable to seek CnySeq\n");
	exit(1);
    }

    if (fwrite(&cnySeqWriteInfo->m_unifiedSeqFileHeader, sizeof(CnyUnifiedSeqFileHeader), 1, cnySeqWriteInfo->m_pFile) != 1) {
	velvetLog("Unable to write CnySeq\n");
	exit(1);
    }

    if (fclose(cnySeqWriteInfo->m_pFile) < 0) {
	velvetLog("Unable to close CnySeq\n");
	exit(1);
    }
}
