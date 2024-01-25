#pragma once

class RPResult
	{
public:
	string m_QueryLabel;			// 查询序列名
	uint m_QL;						// 查询序列长度
	uint m_QLnt;					// 序列长度（核苷酸）
	string m_Gene;					// "RdRP","RT"
	bool m_HiConf;					// 是否为高置信度序列
	string m_Comments;				// 奖罚得分过程中的标记
	string m_XXX;					// Motif 顺序（"ABC","CAB",""）
	int m_Frame;
	float m_PSSMTotalScore;			// 最佳 Group 所有 Motif PSSM 的得分之和
	float m_PSSMMinScore;			// 最佳 Group 中最低的 Motif PSSM 得分
	uint m_StartPos;				// 序列上比对的起始位置
	uint m_SegLength;				// 序列比对的跨度
	uint m_StartPosNt;				
	uint m_SegLengthNt;
	int m_Dist12;					// motif1 和 motif2 的距离
	int m_Dist23;					// motif2 和 motif3 的距离
	double m_FinalScore;			// 序列最终得分
	string m_TopGroupName;			// 最佳 Group 名称
	string m_TrimmedSeq;			// 修剪后的查询序列 (第一个 motif 的第 1 位氨基酸 -> 最后一个 motif 的第最后 1 位氨基酸)
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
