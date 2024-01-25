#ifndef seqdb_h
#define seqdb_h

// SeqDB ��Ķ��� (��fasta�ļ��ж�ȡ��ǩ������)
class SeqDB
	{
public:
	bool m_IsAligned;	// �Ƿ����
	bool m_IsNucleo;	// �Ƿ�Ϊ���������� or ����������
	bool m_IsNucleoSet;	// �Ƿ�Ϊ�����Ἧ��
	unsigned m_ColCount;// ��������Զ�������м�
	vector<string> m_Labels;	// ���б��� �ַ�������
	vector<string> m_Seqs;		// ���� �ַ�������

public:
	SeqDB()
		{
		m_IsAligned = false;
		m_IsNucleo = false;
		m_IsNucleoSet = false;
		m_ColCount = UINT_MAX;
		}

	// �� SeqDB ���������һ�� fasta ����
	unsigned AddSeq(const string &Label, const string &Seq);
	

	const string &GetSeq(unsigned SeqIndex) const;
	const string &GetLabel(unsigned SeqIndex) const;
	unsigned GetSeqLength(unsigned SeqIndex) const;

	// ��ȡ m_IsAligned
	bool IsAligned() const;

	// ��ȡ m_ColCount
	unsigned GetColCount() const;
	
	// ��ȡ m_IsNucleo
	bool GetIsNucleo();

	// ��ȡ��������
	unsigned GetSeqCount() const { return SIZE(m_Seqs); }
	
	// �� fasta �ļ��ж�ȡ���б�ǩ������
	void FromFasta(const string &FileName);

	void WritePretty(FILE *f) const;
	void WriteMSAPretty(FILE *f) const;
	void LogMe() const;

private:
	void SetIsNucleo();
	};

#endif // seqdb_h
