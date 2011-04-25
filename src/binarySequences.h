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
#ifndef _BINARYSEQ_H_
#define _BINARYSEQ_H_

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

// Reading
ReadSet *importCnyReadSet(char *filename);

// Writing
struct binarySequencesWriter_st {
        FILE *		m_pFile;
	FILE *		m_nameFile;
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
typedef struct binarySequencesWriter_st BinarySequencesWriter;
BinarySequencesWriter * openCnySeqForWrite(const char *unifiedSeqFileName);
void cnySeqInsertSequenceName(const char *name, IDnum readID, BinarySequencesWriter *cnySeqWriteInfo); 
void cnySeqInsertNucleotideString(const char *pReadBuf, BinarySequencesWriter *cnySeqWriteInfo);
void inputCnySeqFileStart(Category category, BinarySequencesWriter *cnySeqWriteInfo);
void cnySeqInsertStart(BinarySequencesWriter *cnySeqWriteInfo);
void cnySeqInsertEnd(BinarySequencesWriter *cnySeqWriteInfo);
void closeCnySeqForWrite(BinarySequencesWriter *cnySeqWriteInfo);
#endif
