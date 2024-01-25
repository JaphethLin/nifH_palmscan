#ifndef sfasta_h
#define sfasta_h

#include "myutils.h"
#include <map>

class SeqInfo;

// Sequential reader for FASTA file format.
// Serves sequences in file order to save memory.
// Caches biggish chunks to compromise memory vs. speed.

// FASTA �ļ���ʽ��˳���ȡ��
// ���ļ�˳���ṩ�����Խ�ʡ�ڴ档
// ����ϴ�Ŀ��Խ����ڴ���ٶȡ�

class SFasta
	{
public:
	string m_FileName;		// �ļ���
	FILE *m_File;			// �ļ����
	bool m_AllowGaps;		// ��ȡ�ļ��������Ƿ��� Gap

	uint64 m_FileSize;		// �ļ���С

// Position to start next read
	uint64 m_FilePos;		// ��ʼ��һ�ζ�ȡ��λ��

// Cached data.
	char *m_Buffer;			// ���ݻ����� Buffer

// Bytes allocated to m_Buffer
	uint32 m_BufferSize;	// ����� Buffer ���ֽ���

// Current position in buffer, normally points to '>'
	uint32 m_BufferOffset;	// �������еĵ�ǰλ�ã�ͨ��ָ��'>'

// File data in buffer <= m_BufferSize
	uint32 m_BufferBytes;

// Current label
// Points into m_Buffer, not a separate buffer.
	char *m_Label;

// Current sequence length
	unsigned m_SeqLength;

// Current seq index
	unsigned m_SeqIndex;

	unsigned m_LongestLength;
	unsigned m_TooLongCount;

protected:
	bool m_IsNucleoSet;
	bool m_IsNucleo;
	map<char, unsigned> m_BadCharToCount;

public:
	SFasta();
	~SFasta();

	void Clear();
	void Open(const string &FileName);
	void Rewind();
	bool GetIsNucleo();

// Get next sequence. ��ȡ��һ������
// Returns zero on end-of-file ��ȡ���ļ���βʱ���� 0 
	const char *GetNextSeq();

// Length of most recent sequence returned by GetNextSeq().
	unsigned GetSeqLength() const { return m_SeqLength; }

// Label of most recent sequence returned by GetNextSeq().
	const char *GetLabel() const { return m_Label; }

// Index of most recent sequence returned by GetNextSeq().
	unsigned GetSeqIndex() const { return m_SeqIndex; }

	bool GetNextSI(SeqInfo &SI);

	unsigned GetPctDoneX10() const;
	double GetPctDone() const;

	void LogMe() const;

protected:
	void FillCache();
	const char *GetNextSeqLo();
	};

#endif // sfasta_h
