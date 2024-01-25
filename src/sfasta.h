#ifndef sfasta_h
#define sfasta_h

#include "myutils.h"
#include <map>

class SeqInfo;

// Sequential reader for FASTA file format.
// Serves sequences in file order to save memory.
// Caches biggish chunks to compromise memory vs. speed.

// FASTA 文件格式的顺序读取器
// 按文件顺序提供序列以节省内存。
// 缓存较大的块以降低内存和速度。

class SFasta
	{
public:
	string m_FileName;		// 文件名
	FILE *m_File;			// 文件句柄
	bool m_AllowGaps;		// 读取文件的序列是否含有 Gap

	uint64 m_FileSize;		// 文件大小

// Position to start next read
	uint64 m_FilePos;		// 开始下一次读取的位置

// Cached data.
	char *m_Buffer;			// 数据缓存区 Buffer

// Bytes allocated to m_Buffer
	uint32 m_BufferSize;	// 分配给 Buffer 的字节数

// Current position in buffer, normally points to '>'
	uint32 m_BufferOffset;	// 缓冲区中的当前位置，通常指向'>'

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

// Get next sequence. 获取下一条序列
// Returns zero on end-of-file 读取到文件结尾时返回 0 
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
