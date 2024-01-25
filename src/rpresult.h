#pragma once

class RPResult
	{
public:
	string m_QueryLabel;			// ��ѯ������
	uint m_QL;						// ��ѯ���г���
	uint m_QLnt;					// ���г��ȣ������ᣩ
	string m_Gene;					// "RdRP","RT"
	bool m_HiConf;					// �Ƿ�Ϊ�����Ŷ�����
	string m_Comments;				// �����÷ֹ����еı��
	string m_XXX;					// Motif ˳��"ABC","CAB",""��
	int m_Frame;
	float m_PSSMTotalScore;			// ��� Group ���� Motif PSSM �ĵ÷�֮��
	float m_PSSMMinScore;			// ��� Group ����͵� Motif PSSM �÷�
	uint m_StartPos;				// �����ϱȶԵ���ʼλ��
	uint m_SegLength;				// ���бȶԵĿ��
	uint m_StartPosNt;				
	uint m_SegLengthNt;
	int m_Dist12;					// motif1 �� motif2 �ľ���
	int m_Dist23;					// motif2 �� motif3 �ľ���
	double m_FinalScore;			// �������յ÷�
	string m_TopGroupName;			// ��� Group ����
	string m_TrimmedSeq;			// �޼���Ĳ�ѯ���� (��һ�� motif �ĵ� 1 λ������ -> ���һ�� motif �ĵ���� 1 λ������)
	string m_TrimmedSeqNt;
	string m_MotifsSeq;
	string m_MotifsSeq2;
	string m_SuperMotif;			// SuperMotif
	vector<string> m_Aln;			// 1.[TopLine] 2.[QLine] 3.[ALine] 4.[PLine]		/// m_Aln:[TopLine][QLine][ALine][PLine]
	vector<string> m_FullAln;		// empty											///	[0] '    'A:4-19(2.4)'  ''   'B:27-47(1.2)'      ''   'B:56-67(9.0)'
																						/// [1] '    'ABACCDJSHMALDSP <8> LKSPAIEMNDSJKASJJDNF <9> IFDUBAIFBIA''' '[63]
																						/// [2] '    '|||||+||||+||.|'   '|||||||+++|||.||+|||'   '|||||||||||''
																						/// [3] '    'ABACCGJSHMARDSP'   'LKSPAIEMFASJKASJTTTF'   'IFDUBAIFBIA''
											


public:
	void Clear()
		{
		m_QueryLabel.clear();
		m_QL = UINT_MAX;
		m_QLnt = UINT_MAX;
		m_Gene.clear();
		m_HiConf = false;
		m_Comments.clear();
		m_XXX.clear();

		m_StartPos = UINT_MAX;
		m_StartPosNt = UINT_MAX;
		m_SegLength = UINT_MAX;
		m_SegLengthNt = UINT_MAX;
		m_Frame = -999;

		m_PSSMTotalScore = -999;
		m_PSSMMinScore = -999;
		m_Dist12 = -999;
		m_Dist23 = -999;
		m_SegLength = UINT_MAX;
		m_TopGroupName.clear();
		m_TrimmedSeq.clear();
		m_TrimmedSeqNt.clear();
		m_MotifsSeq.clear();
		m_MotifsSeq2.clear();
		m_SuperMotif.clear();
		m_Aln.clear();
		m_FullAln.clear();
		m_FinalScore = -999;
		}

	void ToFEVStr(string &s) const;
	void GetCat(string &Cat) const
		{
		if (m_Gene == "unclassified")
			Cat = "unclassified";
		else
			{
			if (m_HiConf)
				Cat = "high-confidence-" + m_Gene;
			else
				Cat = "low-confidence-" + m_Gene;
			}		
		}
	};
