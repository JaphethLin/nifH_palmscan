#ifndef seqdb_h
#define seqdb_h

// SeqDB 类的定义 (从fasta文件中读取标签和序列)
class SeqDB
	{
public:
	bool m_IsAligned;	// 是否对齐
	bool m_IsNucleo;	// 是否为核苷酸序列 or 氨基酸序列
	bool m_IsNucleoSet;	// 是否为核苷酸集合
	unsigned m_ColCount;// 列数，针对对齐的序列集
	vector<string> m_Labels;	// 序列标题 字符串向量
	vector<string> m_Seqs;		// 序列 字符串向量

public:
	SeqDB()
		{
		m_IsAligned = false;
		m_IsNucleo = false;
		m_IsNucleoSet = false;
		m_ColCount = UINT_MAX;
		}

	// 向 SeqDB 对象中添加一条 fasta 序列
	unsigned AddSeq(const string &Label, const string &Seq);
	

	const string &GetSeq(unsigned SeqIndex) const;
	const string &GetLabel(unsigned SeqIndex) const;
	unsigned GetSeqLength(unsigned SeqIndex) const;

	// 获取 m_IsAligned
	bool IsAligned() const;

	// 获取 m_ColCount
	unsigned GetColCount() const;
	
	// 获取 m_IsNucleo
	bool GetIsNucleo();

	// 获取序列数量
	unsigned GetSeqCount() const { return SIZE(m_Seqs); }
	
	// 从 fasta 文件中读取序列标签和序列
	void FromFasta(const string &FileName);

	void WritePretty(FILE *f) const;
	void WriteMSAPretty(FILE *f) const;
	void LogMe() const;

private:
	void SetIsNucleo();
	};

#endif // seqdb_h
