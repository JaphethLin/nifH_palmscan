#include "myutils.h"
#include "rdrpmodel.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "alnparams.h"
#include "pathinfo.h"

void SeqToFasta(FILE *f, const char *Label, const char *Seq, unsigned L);

/***
Typical spec line:

A	Duplorna	e:/res/.../model1.pssm
B	Burksadf	/res/.../model2.pssm


Fields:
	1. MotifLetter (X)
	2. Short name for reports (.X appended)
	3. Path to model file
***/

void RdRpModel::FromSpecFile(const string &FileName)
	{
	Clear();

	m_MotifLetters.push_back('A');
	m_MotifLetters.push_back('B');
	m_MotifLetters.push_back('C');

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	uint PSSMCount = 0;
	while (ReadLineStdioFile(f, Line))	// 单行读取文件，若行为空或以 ‘#’ 开头，则跳过
		{
		if (StartsWith(Line, "#"))	
			continue;
		if (Line.empty())
			continue;

		Split(Line, Fields, '\t');		// 以 \t 分割每行字段
		asserta(SIZE(Fields) == 3);		
		asserta(SIZE(Fields[0]) == 1);
		char MotifLetter = toupper(Fields[0][0]);	// 读取 Motif 字段
		const string &GroupName = Fields[1];		// 读取 group 字段
		const string &FastaFileName = Fields[2];	// 读取 文件路径 字段

		bool Found = false;							// 判定该行指定的 PSSM 所属的 Motif 是否已被存入，若没有被发现则存入，若已经被存入则跳过
		for (uint i = 0; i < SIZE(m_GroupNames); ++i)
			{
			if (m_GroupNames[i] == GroupName)
				{
				Found = true;
				break;
				}
			}
		if (!Found)
			m_GroupNames.push_back(GroupName);			

		m_PSSMMotifLetters.push_back(MotifLetter);		// 存入该行指定的 PSSM 对应 Motif 的字符
		m_PSSMGroupNames.push_back(GroupName);			// 存入该行指定的 PSSM 对应 Group 的字符
		m_PSSMs.resize(PSSMCount+1);				
		m_PSSMs[PSSMCount].FromFasta(FastaFileName);	// 从先前计算好的 PSSM 文件中读取模型信息，存入 PSSMs 向量中
		++PSSMCount;
		}
	CloseStdioFile(f);

	SetPIV();
	}

void RdRpModel::SetPIV()	// 索引 PSSM 模型和 Motif 与 Group 之间的关系
	{
	const uint MotifCount = SIZE(m_MotifLetters);
	const uint GroupCount = SIZE(m_GroupNames);
	m_PIV.resize(MotifCount);
	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		{
		m_PIV[MotifIndex].resize(GroupCount);
		char MotifLetter = m_MotifLetters[MotifIndex];
		for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
			{
			const string &GroupName = m_GroupNames[GroupIndex];
			uint PSSMIndex = GetPSSMIndex(MotifLetter, GroupName);
			m_PIV[MotifIndex][GroupIndex] = PSSMIndex;
			}
		}
	}

void RdRpModel::SearchAA1(uint PSSMIndex, const string &Seq,
  RPHit &Hit)
	{
	CheckThread();
	if (m_Gapped)
		SearchAA1_Gapped(PSSMIndex, Seq, Hit);			// Gapped 检索
	else
		SearchAA1_Ungapped(PSSMIndex, Seq, Hit);		// Ungapped 检索
	
	extern FILE *g_fPSSMAln;
	if (g_fPSSMAln == 0)	// 若输出文件找不到，则不输出
		return;
	if (Hit.m_Score < 0)	// 若比对得分小于 0 则不输出当前比对
		return;

	FILE *f = g_fPSSMAln;	// 获得文件句柄

	const string &GroupName = m_PSSMGroupNames[PSSMIndex].c_str();	// 获得当前 PSSM 的 Group 名称
	char MotifLetter = m_PSSMMotifLetters[PSSMIndex];				// 获得当前 PSSM 的 Motif 标签

	uint GroupIndex = UINT_MAX;
	uint MotifIndex = UINT_MAX;
	for (uint g = 0; g < SIZE(m_GroupNames); ++g)					// 获取当前 PSSMIndex 对应的 GroupIndex 和 MotifIndex
		{
		for (uint m = 0; m < SIZE(m_MotifLetters); ++m)
			{
			uint PSSMIndex2 = GetPSSMIndex(m, g);
			if (PSSMIndex2 == PSSMIndex)
				{
				GroupIndex = g;
				MotifIndex = m;
				break;
				}
			}
		}
	asserta(GroupIndex != UINT_MAX);
	asserta(MotifIndex != UINT_MAX);
	
	string Q;
	string P;
	string A;
	uint QLo;
	uint QHi;

	///
	/// \n
	/// PSSM Group1(A) 2.3
	///
	fprintf(f, "\n");
	fprintf(f, "PSSM %s(%c) %.1f\n",
	  GroupName.c_str(), MotifLetter, Hit.m_Score);
	GetAln(MotifIndex, GroupIndex, Q, P, A, QLo, QHi);		
	/// Q = Query 序列命中片段的氨基酸序列
	/// P = 该 PSSM 模型的保守序列
	/// A = 比对的 Score 序列 '|'(>=1)  '+'(>=0.5)  '.'(>0)  " "(<=0)	
	/// QLo = Query 序列命中片段的第一个氨基酸的位置数  
	/// QLo = Query 序列命中片段的最后一个氨基酸的位置数
	fprintf(f, "%s\n", Q.c_str());
	fprintf(f, "%s\n", A.c_str());
	fprintf(f, "%s\n", P.c_str());
	}

// 无 Gap 氨基酸序列检索
void RdRpModel::SearchAA1_Ungapped(uint PSSMIndex, const string &Seq,
  RPHit &Hit)
	{
	CheckThread();
	const uint L = SIZE(Seq);						// L -> 序列长度
	const PSSM &P = GetPSSM(PSSMIndex);				// P -> PSSM
	PSSMHitAA HitA;
	GetTopHitAA(P, -9, "_notused_", Seq.c_str(), L, HitA);	 // 获得片段最佳得分和起始位置
	Hit.SetUngapped(P, HitA.AAPos, HitA.Score);		// 通过 HitA (PSSMHitAA 结构) 给 Hit (RPHit 类) 赋值
	}

// 有 Gap 氨基酸序列检索
void RdRpModel::SearchAA1_Gapped(uint PSSMIndex, const string &Seq,
  RPHit &Hit)
	{
	CheckThread();
	// 定义 TrimQueryTermGaps 函数
	void TrimQueryTermGaps(const char *Path,
	  unsigned &QLo, unsigned &QHi,
	  unsigned &ColLo, unsigned &ColHi);

	// 定义 Viterbi_PSSM 函数
	float Viterbi_PSSM(XDPMem &Mem, const char *A, unsigned LA,
	  const PSSM &P, const AlnParams &AP, PathInfo &PI);

	const uint L = SIZE(Seq);
	const PSSM &P = GetPSSM(PSSMIndex);
	PathInfo *PI = m_OM.GetPathInfo();
	float Score = Viterbi_PSSM(m_DPMem, Seq.c_str(), L, P, m_AP, *PI);

	uint QLo, QHi, ColLo, ColHi;
	TrimQueryTermGaps(PI->GetPath(), QLo, QHi, ColLo, ColHi);
	uint PL = P.GetColCount();
	Hit.SetGapped(P, QLo, Score, PI, ColLo, ColHi);
	
	m_OM.Down(PI);
	}

// 获取 Motif2 比对首位到 Motif1 比对末尾的距离
int RdRpModel::GetDist(char MotifLetter1, char MotifLetter2) const
	{
	const RPHit *Hit1 = GetHit(MotifLetter1);
	const RPHit *Hit2 = GetHit(MotifLetter2);
	int Dist = GetDist(Hit1, Hit2);
	return Dist;
	}

int RdRpModel::GetDist(const RPHit *Hit1, const RPHit *Hit2) const
	{
	if (Hit1 == 0 || Hit2 == 0)
		return -999;
//	uint PL1 = Hit1->m_PSSM->GetColCount();		// Lim Deleted
	uint PL1 = Hit1->GetQuerySegLength();		// Lim Edited
	uint QPos1 = Hit1->m_QPos;
	uint QPos2 = Hit2->m_QPos;
	int Dist = int(QPos2) - int(QPos1+PL1);
	return Dist;
	}

// 获取 Motif1 比对首位到 Motif2/3 比对末尾的距离
uint RdRpModel::GetSpan(const RPHit *Hit1, const RPHit *Hit2) const
	{
	if (Hit1 == 0 || Hit2 == 0)
		return -999;
	uint PL2 = Hit2->m_PSSM->GetColCount();
	uint QStart1 = Hit1->m_QPos;
	uint QStart2 = Hit2->m_QPos;
	if (QStart2 <= QStart1)
		return 0;
	uint QEnd2 = QStart2 + Hit2->GetQuerySegLength();
	uint Span = QEnd2 - QStart1;
	return Span;
	}

// 获得该 Motif 在 Query 序列上找到的片段的氨基酸序列片段 (根据 RPHit 对象 Path 中的 M 和 D )
void RdRpModel::GetXSeq(char MotifLetter, string &XSeq) const
	{
	XSeq = "";
	uint GroupIndex = GetTopGroup();	// 最佳 Group
	if (GroupIndex == UINT_MAX)
		{
		XSeq = "";
		return;
		}

	const uint QL = SIZE(m_QuerySeq);
	const uint PL = GetPSSMLength(MotifLetter, GroupIndex);
	const RPHit *Hit = GetHit(MotifLetter, GroupIndex);
	if (Hit == 0)				// Hit 对象不为空
		{
		XSeq = "";
		return;
		}
	uint QPos = Hit->m_QPos;	
	const string &Path = Hit->m_Path;
	for (uint Col = 0; Col < SIZE(Path); ++Col)
		{
		char Op = Path[Col];
		switch (Op)
			{
		case 'M':
			{
			asserta(QPos < QL);
			char q = m_QuerySeq[QPos];
			XSeq += q;
			++QPos;
			break;
			}

		case 'D':
			{
			asserta(QPos < QL);
			char q = m_QuerySeq[QPos++];
			XSeq += q;
			break;
			}

		case 'I':
			{
//			XSeq += '-';
			break;
			}

		default:
			asserta(false);
			}
		}
	}

void RdRpModel::GetAln(uint MotifIndex, uint GroupIndex,
  string &QRow, string &PRow, string &AnnotRow,
  uint &QLo, uint &QHi) const
	{
	CheckThread();

	QRow.clear();
	PRow.clear();
	AnnotRow.clear();

	const uint QL = SIZE(m_QuerySeq);							// QL：查询序列的长度
	const uint PL = GetPSSMLength(MotifIndex, GroupIndex);		// PL：PSSM 长度
	const RPHit *Hit = GetHit(MotifIndex, GroupIndex);			// 获得 RPHit 对象
	const PSSM *P = GetPSSM(MotifIndex, GroupIndex);			// PSSM 对象
	if (Hit == 0)
		return;
	uint QPos = Hit->m_QPos;		// 序列命中起始位置
	QLo = QPos + 1;
	uint PPos = 0;					// PSSM 第一个位置

	const string &Path = Hit->m_Path;
	for (uint Col = 0; Col < SIZE(Path); ++Col)
		{
		char Op = Path[Col];
		switch (Op)
			{
		case 'M':	// Match //
			{
			asserta(QPos < QL);
			asserta(PPos < PL);
			char q = m_QuerySeq[QPos];				// 查询序列命中片段中 QPos 位置的氨基酸
			float Score = P->GetScore1(PPos, q);	// QPos 位置的氨基酸在 PSSM 中 PPos 位置上的 Score
			char p = P->GetConsChar(PPos);			// PSSM 中 PPos 位置上的保守氨基酸 

			if (Score >= 1)							// AnnotRow : " ||++.++|||+|| " 
				AnnotRow += '|';					// [ >=1 : | ]  [ >=0.5 : + ] [ >0 : . ] [ <=0 : " " ]
			else if (Score >= 0.5)
				AnnotRow += '+';
			else if (Score > 0)
				AnnotRow += ".";
			else
				AnnotRow += " ";

			QRow += q;		// QRow : Query 序列氨基酸序列
			PRow += p;		// PRow : PSSM 保守氨基酸序列

			++QPos;			// QPos 后移一位
			++PPos;			// PPos 后移一位
			break;
			}

		case 'D':	// Delete //	PSSM 相对与 Query 序列，在该位置上有一个删除
			{
			asserta(QPos < QL);
			char q = m_QuerySeq[QPos++];

			QRow += q;
			PRow += '-';
			AnnotRow += " ";
			break;
			}

		case 'I':	// Insert //	PSSM 相对与 Query 序列，在该位置上有一个插入
			{
			asserta(PPos < PL);
			char p = P->GetConsChar(PPos++);

			QRow += '-';
			PRow += p;
			AnnotRow += " ";
			break;
			}

		default:
			asserta(false);
			}
		}
	QHi = QPos;
	}

void RdRpModel::LogAlnViterbi(uint PSSMIndex, const string &QSeq, const RPHit &Hit) const
	{
	const uint QL = SIZE(QSeq);

	const PSSM &P = GetPSSM(PSSMIndex);
	const uint PL = P.GetColCount();

	string QRow;
	string PRow;
	string AnnotRow;

	uint QPos = Hit.m_QPos;
	uint PPos = 0;
	const string &Path = Hit.m_Path;
	for (uint Col = 0; Col < SIZE(Path); ++Col)
		{
		char Op = Path[Col];
		switch (Op)
			{
		case 'M':
			{
			asserta(QPos < QL);
			asserta(PPos < PL);
			char q = QSeq[QPos];
			float Score = P.GetScore1(PPos, q);
			char p = P.GetConsChar(PPos);

			if (Score >= 1)
				AnnotRow += '|';
			else if (Score >= 0.5)
				AnnotRow += '+';
			else if (Score > 0)
				AnnotRow += ".";
			else
				AnnotRow += " ";

			QRow += q;
			PRow += p;

			++QPos;
			++PPos;
			break;
			}

		case 'D':
			{
			asserta(QPos < QL);
			char q = QSeq[QPos++];

			QRow += q;
			PRow += '-';
			break;
			}

		case 'I':
			{
			asserta(PPos < PL);
			char p = P.GetConsChar(PPos++);

			QRow += '-';
			PRow += p;
			break;
			}

		default:
			asserta(false);
			}
		}

	char MotifLetter = m_PSSMMotifLetters[PSSMIndex];
	Log("\n");
	Log("Q  %s\n", QRow.c_str());
	Log("   %s\n", AnnotRow.c_str());
	Log("%c  %s\n", MotifLetter, PRow.c_str());
	Log("\n");
	}

// 检索氨基酸序列
void RdRpModel::SearchAA(const string &QueryLabel, const string &QuerySeq)
	{
	CheckThread();

	m_QueryLabel = QueryLabel;		// 待检索序列的标签
	m_QuerySeq = QuerySeq;			// 待检索的序列

	const uint PSSMCount = GetPSSMCount();		// PSSM 模型的数量
	m_HitsU.resize(PSSMCount);					// 重置无Gap命中对象为模型数量对应大小
	m_HitsG.resize(PSSMCount);					// 重置有Gap命中对象为模型数量对应大小

	for (uint PSSMIndex = 0; PSSMIndex < PSSMCount; ++PSSMIndex)			// 循环使用各个 PSSM 模型检索序列
		{
		RPHit &Hit = (m_Gapped ? m_HitsG[PSSMIndex] : m_HitsU[PSSMIndex]);	// 判定是 Gapped 还是 Ungapped 检索
		SearchAA1(PSSMIndex, m_QuerySeq, Hit);								// 获取 Query 序列与各个 PSSM 模型比对的命中结果 Hit
		}

	CheckThread();
	SetResult();		// 核心计算函数，计算 m_Result（RPResult 类） 
	}

uint RdRpModel::GetPSSMIndex(char MotifLetter, const string &GroupName) const
	{
	const uint PSSMCount = GetPSSMCount();
	for (uint i = 0; i < PSSMCount; ++i)
		{
		if (m_PSSMMotifLetters[i] == MotifLetter &&
		  m_PSSMGroupNames[i] == GroupName)
			return i;
		}
	return UINT_MAX;
	}

void RdRpModel::LogHitTable() const
	{
	Log("\n");
	Log(">%s\n", m_QueryLabel.c_str());
	const uint MotifCount = SIZE(m_MotifLetters);
	const uint GroupCount = SIZE(m_GroupNames);
	const uint PSSMCount = SIZE(m_PSSMs);

	Log("%10.10s  ", "");
	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		{
		char MotifLetter = m_MotifLetters[MotifIndex];
		Log("  %12.12s%c", "", MotifLetter);
		if (MotifIndex + 1 != MotifCount)
			Log("  %5.5s", "gap");
		}
	Log("  %5.5s - %5.5s (%s)", "Start", "End", "Len");
	Log("\n");

	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		const string &GroupName = m_GroupNames[GroupIndex];
		Log("  %10.10s", GroupName.c_str());
		uint AStart = UINT_MAX;
		uint CEnd = UINT_MAX;
		for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
			{
			char MotifLetter = m_MotifLetters[MotifIndex];
//			uint PSSMIndex = GetPSSMIndex(MotifLetter, GroupName);
			uint PSSMIndex = GetPSSMIndex(MotifIndex, GroupIndex);
			uint PL = GetPSSMLength(MotifIndex, GroupIndex);
			if (PSSMIndex == UINT_MAX)
				Log("  %13.13s", "<missing>");
			else
				{
				const RPHit &HitG = m_HitsG[PSSMIndex];
				uint QPosG = HitG.m_QPos;
				if (MotifIndex == 0)
					AStart = QPosG;
				else if (MotifIndex + 1 == MotifCount)
					CEnd = QPosG + PL - 1;
				float Score = HitG.m_Score;
				Log("  %5u(%6.2f)", QPosG, Score);

				if (MotifIndex + 1 < MotifCount)
					{
//					uint NextPSSMIndex = GetPSSMIndex(MotifLetter+1, GroupName);
					uint NextPSSMIndex = GetPSSMIndex(MotifIndex+1, GroupIndex);
					const RPHit &NextHitG = m_HitsG[NextPSSMIndex];
					uint NextPos = NextHitG.m_QPos;
					int Gap = int(NextPos) - int(QPosG) - int(PL);
					Log("  %5d", Gap);
					}
				}
			if (AStart != UINT_MAX && CEnd != UINT_MAX)
				{
				uint ACLength = CEnd - AStart + 1;
				Log("  %5u - %5u (%u)", AStart, CEnd, ACLength);
				}
			}
		Log("\n");
		}
	}

uint RdRpModel::GetPSSMLength(uint MotifIndex, uint GroupIndex) const
	{
	const PSSM *P = GetPSSM(MotifIndex, GroupIndex);
	if (P == 0)
		return 0;
	uint L = P->GetColCount();
	return L;
	}

// 获得某 Group 的总得分，即该 Group 下所有 Motif 的 PSSM 得分之和
float RdRpModel::GetTotalScore(uint GroupIndex) const
	{
	float Total = 0;
	for (uint MotifIndex = 0; MotifIndex < GetMotifCount(); ++MotifIndex)
		{
		const RPHit *HitG = GetHit(MotifIndex, GroupIndex);
		if (HitG == 0)
			return 0;
		if (HitG->m_Score < 0)
			return 0;
		Total += HitG->m_Score;
		}
	return Total;
	}

// 获得最佳 Group 中 Motif C 的比对 Score
float RdRpModel::GetCScore() const
	{
	uint GroupIndex = GetTopGroup();	// 获得最佳 Group
	if (GroupIndex == UINT_MAX)			// 若最佳 Group 不存在(均小于0)，则返回 -999
		return -999;
	float Min = 999;
	uint MotifIndex = GetMotifIndex('C');
	const RPHit *HitG = GetHit(MotifIndex, GroupIndex);		// 获得 Motif C 的比对结果
	if (HitG == 0)						// 若 Motif C 的比对结果不存在，则返回 -999
		return -999;					
	float Score = HitG->m_Score;		// 获得 Motif C 的比对 Score
	return Score;
	}

// Lim edited //
// 获得最佳 Group 中 Motif B 的比对 Score
float RdRpModel::GetBScore() const
{
	uint GroupIndex = GetTopGroup();	// 获得最佳 Group
	if (GroupIndex == UINT_MAX)			// 若最佳 Group 不存在(均小于0)，则返回 -999
		return -999;
	float Min = 999;
	uint MotifIndex = GetMotifIndex('B');
	const RPHit* HitG = GetHit(MotifIndex, GroupIndex);		// 获得 Motif B 的比对结果
	if (HitG == 0)						// 若 Motif B 的比对结果不存在，则返回 -999
		return -999;
	float Score = HitG->m_Score;		// 获得 Motif B 的比对 Score
	return Score;
}

// Lim edited //
// 获得最佳 Group 中 Motif A 的比对 Score
float RdRpModel::GetAScore() const
{
	uint GroupIndex = GetTopGroup();	// 获得最佳 Group
	if (GroupIndex == UINT_MAX)			// 若最佳 Group 不存在(均小于0)，则返回 -999
		return -999;
	float Min = 999;
	uint MotifIndex = GetMotifIndex('A');
	const RPHit* HitG = GetHit(MotifIndex, GroupIndex);		// 获得 Motif A 的比对结果
	if (HitG == 0)						// 若 Motif A 的比对结果不存在，则返回 -999
		return -999;
	float Score = HitG->m_Score;		// 获得 Motif A 的比对 Score
	return Score;
}


// 获得最佳 Group 中 Motif 比对的最小 Score
float RdRpModel::GetMinScore() const
	{
	uint GroupIndex = GetTopGroup();
	if (GroupIndex == UINT_MAX)			// 若最佳 Group 不存在(均小于0)，则返回 -999
		return -999;
	float Min = 999;
	for (uint MotifIndex = 0; MotifIndex < GetMotifCount(); ++MotifIndex)
		{
		const RPHit *HitG = GetHit(MotifIndex, GroupIndex);
		if (HitG == 0)					// 若有一种 Motif 的比对结果不存在，则返回 -1
			return -1;
		Min = min(Min, HitG->m_Score);
		}
	return Min;
	}

uint RdRpModel::GetTotalGaps() const
	{
	uint GroupIndex = GetTopGroup();
	if (GroupIndex == UINT_MAX)
		return 0;
	uint Total = 0;
	for (uint MotifIndex = 0; MotifIndex < GetMotifCount(); ++MotifIndex)
		{
		const RPHit *HitG = GetHit(MotifIndex, GroupIndex);
		if (HitG == 0)
			return 0;
		Total += HitG->GetGapCount();
		}
	return Total;
	}

// 返回 m_HitsU[PSSMIndex] 对象（RPHit 类）
const RPHit *RdRpModel::GetHit(uint MotifIndex, uint GroupIndex) const	
	{
	uint PSSMIndex = GetPSSMIndex(MotifIndex, GroupIndex);
	if (PSSMIndex == UINT_MAX)
		return 0;
	asserta(PSSMIndex < SIZE(m_HitsG));
	if (m_Gapped)
		return &m_HitsG[PSSMIndex];
	return &m_HitsU[PSSMIndex];
	}

// 获得最高得分的 Group 的 Index 和 Score，若总 Score <= 0 则返回 UINT_MAX
void RdRpModel::GetTopGroup(uint &TopIx, float &TopScore) const
	{
	TopIx = UINT_MAX;
	TopScore = 0;

	for (uint GroupIndex = 0; GroupIndex < GetGroupCount(); ++GroupIndex)
		{
		float Score = GetTotalScore(GroupIndex);	// 计算 group 总得分
		if (Score > TopScore)
			{
			TopIx = GroupIndex;			
			TopScore = Score;
			}
		}
	}

bool RdRpModel::GetABCRange(uint &Start, uint &End) const
	{
	Start = UINT_MAX;
	End = UINT_MAX;
	uint GroupIndex = GetTopGroup();
	if (GroupIndex == UINT_MAX)
		return false;

	const RPHit *HitA = GetHit('A', GroupIndex);
	const RPHit *HitB = GetHit('B', GroupIndex);
	const RPHit *HitC = GetHit('C', GroupIndex);
	if (HitA == 0 || HitB == 0 || HitC == 0)
		{
		Start = 0;
		End = 0;
		return false;
		}

	uint PosA = HitA->m_QPos;
	uint PosB = HitB->m_QPos;
	uint PosC = HitC->m_QPos;

	uint QAL = HitA->GetQuerySegLength();  // Query 序列中该 Hit 的长度
	uint QBL = HitB->GetQuerySegLength();
	uint QCL = HitC->GetQuerySegLength();

	uint HiA = PosA + QAL - 1;
	uint HiB = PosB + QBL - 1;
	uint HiC = PosC + QCL - 1;

	Start = min(min(PosA, PosB), PosC);		// 三个 motif 在序列上的最小起始位置
	End = max(max(HiA, HiB), HiC);			// 三个 motif 在序列上的最大结束位置

	const uint QL = SIZE(m_QuerySeq);		// Query 序列的总长度
	if (End >= QL)
		{
		if (QL - End <= 2)
			End = QL - 1;
		else
			{
			Start = 0;
			End = 0;
			return false;
			}
		}
	return true;	// 命中的结束位置 < 序列最后一个氨基酸位置 返回 true
	}

// 获得 ABC motif 的顺序字符串
void RdRpModel::GetABCOrder(string &XXX) const
	{
	const RPHit *ptrHit1;
	const RPHit *ptrHit2;
	const RPHit *ptrHit3;
	GetOrderedHits(XXX, ptrHit1, ptrHit2, ptrHit3);
	}

void RdRpModel::GetOrderedHits(string &XXX, const RPHit *&Hit1,
  const RPHit *&Hit2, const RPHit *&Hit3) const
	{
	XXX.clear();
	Hit1 = 0;
	Hit2 = 0;
	Hit3 = 0;

	uint GroupIndex = GetTopGroup();
	if (GroupIndex == UINT_MAX)
		return;

	const RPHit *HitA = GetHit('A', GroupIndex);
	const RPHit *HitB = GetHit('B', GroupIndex);
	const RPHit *HitC = GetHit('C', GroupIndex);
	if (HitA == 0 || HitB == 0 || HitC == 0)  // 未找到 PSSMIndex 对应的 RPHit 对象
		return;

	uint PosA = HitA->m_QPos;
	uint PosB = HitB->m_QPos;
	uint PosC = HitC->m_QPos;

/** Lim deleted

	if (PosA == PosB || PosB == PosC || PosA == PosC)
		return;

	if (PosA < PosB && PosB < PosC)
		{
		XXX = "ABC";
		Hit1 = HitA;
		Hit2 = HitB;
		Hit3 = HitC;
		}
	else if (PosC < PosA && PosA < PosB)
		{
		XXX = "CAB";
		Hit1 = HitC;
		Hit2 = HitA;
		Hit3 = HitB;
		}

**/

	if (PosA >= PosB || PosB >= PosC || PosA >= PosC)
		return;
	else
		{
			XXX = "ABC";
			Hit1 = HitA;
			Hit2 = HitB;
			Hit3 = HitC;
		}
	}

// 通过 Motif 的标签获得 Motif 的索引，若未找到，则返回 UINT_MAX
uint RdRpModel::GetMotifIndex(char MotifLetter) const
	{
	for (uint MotifIndex = 0; MotifIndex < GetMotifCount(); ++MotifIndex)
		{
		if (m_MotifLetters[MotifIndex] == MotifLetter)
			return MotifIndex;
		}
	return UINT_MAX;
	}

// 通过 Group 的名称获得 Group 的索引，若未找到，则返回 0
uint RdRpModel::GetGroupIndex(const string &GroupName) const
	{
	for (uint GroupIndex = 0; GroupIndex < GetGroupCount(); ++GroupIndex)
		{
		if (m_GroupNames[GroupIndex] == GroupName)
			return GroupIndex;
		}
	asserta(false);
	return 0;
	}

// 根据第一个 Motif 比对的起始与最后一个 Motif 比对的结束位置 裁剪 Query 序列
void RdRpModel::GetTrimmedSeq(string &Seq) const
	{
	Seq.clear();

	uint TopGroup = GetTopGroup();		// 最佳 Group
	if (TopGroup == UINT_MAX)
		return;

	uint Start, End;
	bool Ok = GetABCRange(Start, End);	// 获取 Query 序列上的起始和结束位置
	if (!Ok)
		return;

	const uint QL = SIZE(m_QuerySeq);	// Query 序列长度
	if (End >= QL)
		{
		Die("End %u >= QL %u", End, QL);
		return;
		}

	for (uint Pos = Start; Pos <= End; ++Pos)	// 裁剪 Query 序列
		Seq += m_QuerySeq[Pos];
	}

/// MotifA 匹配到的片段 + "xxx" + MotifB 匹配到的片段 + "xxx" + MotifC 匹配到的片段
void RdRpModel::GetMotifsSeq(string &s) const
	{
	s.clear();

	string A, B, C;
	GetXSeq('A', A);								// 找到 "A" 表示的 Motif 的在 Query 序列上匹配到的氨基酸片段
	GetXSeq('B', B);
	GetXSeq('C', C);
	if (A.empty() || B.empty() || C.empty())		// 任一 Motif 所属的 Group 或 Hit 找不到 
		return;

	s = A + "xxx" + B + "xxx" + C;					// 整合三个 Motif 匹配到的片段
	}

/// 
/// 若 GroupName 为 Duplorna、Kitrino、Lenua 或 Pisu，则对 Motif A 比对得到的 Query 序列片段，在第 8 位插入一个 "-", 并删除第 12 位
/// 
void RdRpModel::GetMotifsSeq2(string &s) const
	{
	s.clear();
	
	uint TopGroup;
	float TopScore;
	GetTopGroup(TopGroup, TopScore);
	if (TopGroup == UINT_MAX)
		return;

	string GroupName;
	GetGroupName(TopGroup, GroupName);

	bool AGap = (GroupName == "Duplorna" ||
				 GroupName == "Kitrino" ||
				 GroupName == "Lenua" ||
				 GroupName == "Pisu");

	string A, B, C;
	GetXSeq('A', A);
	GetXSeq('B', B);
	GetXSeq('C', C);
	if (A.empty() || B.empty() || C.empty())		// 任一 Motif 所属的 Group 或 Hit 找不到 
		return;

	// asserta(SIZE(A) == 12);		// 若比对上的 Motif A 属于以上四种 Group ，那比对得到的 Query 序列片段长度必为 12
	if (AGap)
		{
		string A2;
		for (uint i = 0; i < 12; ++i)
			{
			if (i < 7)
				A2 += A[i];
			else if (i == 7)
				A2 += '-';
			else
				A2 += A[i-1];
			}
		asserta(SIZE(A2) == 12);
		A = A2;
		}

	s = A + B + C;
	}

void RdRpModel::GetFullAln(vector<string> &Rows) const
	{
	Die("GetFullAln not implemented");
	//Rows.clear();
	//Rows.resize(4);

	//string &TopLine = Rows[0];
	//string &QLine = Rows[1];
	//string &ALine = Rows[2];
	//string &PLine = Rows[3];

	//if (m_Result.m_StartPos == UINT_MAX)
	//	return;
	//
	//const uint QL = SIZE(m_QuerySeq);
	//asserta(QL ==  m_Result.m_QL);
	//asserta(m_Result.m_StartPos < m_Result.m_QL);
	//for (uint i = 0; i < QL; ++i)
	//	QLine[i] = m_QuerySeq[i];
	}

/// Rows: [1: TopLine] [2: QLine] [3: ALine] [4: PLine]
void RdRpModel::GetAlnRows(vector<string> &Rows) const
	{
	Rows.clear();
	Rows.resize(4);				// 包含四个字符串

	string &TopLine = Rows[0];	
	string &QLine = Rows[1];
	string &ALine = Rows[2];
	string &PLine = Rows[3];

	uint TopGroup;								// 最佳 Group
	float TopScore;								// 最佳 Group 的 Score
	GetTopGroup(TopGroup, TopScore);
	if (TopGroup == UINT_MAX)		// 若 TopGroup == UINT_MAX 则表示没有匹配到总 Score > 0 的 group
		{
		Rows.clear();
		return;
		}

	string TopGroupName;						// 最佳 Group 的名称
	GetGroupName(TopGroup, TopGroupName);	

	string XXX;									// 序列中 Motif 的顺序 
	GetABCOrder(XXX);							// 若顺序正确为 "ABC" ，否则为 ""

	bool AnyMotifs = false;
	char PrevMotifLetter = 0;
	for (uint k = 0; k < 3; ++k)
		{
		char MotifLetter = XXX[k];
		uint MotifIndex = GetMotifIndex(MotifLetter);	// 根据 MotifLetter 找到对应的索引，若顺序不正确 MotifLetter 为 "" , MotifIndex 为 UINT_MAX
		if (MotifIndex == UINT_MAX)						// 未找到
			continue;
		AnyMotifs = true;

		uint QLo, QHi;							// Query 序列命中的 起始位置（QLo）和 结束位置（QHi）
		string QRow, PRow, ARow;				// Query 序列的氨基酸序列（QRow），PSSM 序列的氨基酸序列（PRow），Query 序列的 Score
		GetAln(MotifIndex, TopGroup, QRow, PRow, ARow, QLo, QHi);	

		uint n = SIZE(QRow);
		string s;
		string t;
		Psa(t, "%c:%u-%u", MotifLetter, QLo, QHi);
		const RPHit *HitG = GetHit(MotifIndex, TopGroup);
		if (HitG != 0)
			Psa(t, "(%.1f)", HitG->m_Score);
		uint m = SIZE(t);
		
		if (PrevMotifLetter != 0)	
			{
			int Dist = GetDist(PrevMotifLetter, MotifLetter);
			Psa(s, " <%d> ", Dist);
			QLine += s;
			for (uint k = 0; k < SIZE(s); ++k)
				{
				TopLine += ' ';
				PLine += ' ';
				ALine += ' ';
				}
			}
		else			// 处理第一个出现的 Motif 
			{
			TopLine += "   ";
			QLine += "   ";
			PLine += "   ";
			ALine += "   ";
			}

		TopLine += t;
		for (uint k = m; k < n; ++k)
			TopLine += ' ';

		Psa(QLine, "%s", QRow.c_str());
		Psa(PLine, "%s", PRow.c_str());
		Psa(ALine, "%s", ARow.c_str());

		if (k < 2)
			{
			for (uint kk = n; kk < m; ++kk)
				{
				QLine += ' ';
				PLine += ' ';
				ALine += ' ';
				}
			}

		PrevMotifLetter = MotifLetter;
		}

	///	'    'A:4-19(2.4)'  ''   'B:27-47(1.2)'      ''   'B:56-67(9.0)'
	/// '    'ABACCDJSHMALDSP <8> LKSPAIEMNDSJKASJJDNF <9> IFDUBAIFBIA''
	/// '    '|||||+||||+||.|'   '|||||||+++|||.||+|||'   '|||||||||||''
	/// '    'ABACCGJSHMARDSP'   'LKSPAIEMFASJKASJTTTF'   'IFDUBAIFBIA''


	if (!AnyMotifs)		// 若找不到任何 Motif 则清空 Rows 并结束当前函数
		{
		Rows.clear();
		return;
		}

	uint Start;			// 第一个 Motif 起始位置
	uint End;			// 最后一个 Motif 结束位置
	bool RangeOk = GetABCRange(Start, End);
	if (RangeOk)
		Psa(QLine, "  [%d]", int(End) - int(Start) + 1);

	///	'    'A:4-19(2.4)'  ''   'B:27-47(1.2)'      ''   'B:56-67(9.0)'
	/// '    'ABACCDJSHMALDSP <8> LKSPAIEMNDSJKASJJDNF <9> IFDUBAIFBIA''' '[63]
	/// '    '|||||+||||+||.|'   '|||||||+++|||.||+|||'   '|||||||||||''
	/// '    'ABACCGJSHMARDSP'   'LKSPAIEMFASJKASJTTTF'   'IFDUBAIFBIA''
	}


/// 从已经比对到的 Query 序列片段中找 SuperMotif 片段
void RdRpModel::GetSuperMotif(const string &MotifsSeq, string &s) const
	{
	s.clear();
	uint n = SIZE(MotifsSeq);
	if (n == 0)
		return;
	if (n != 40)
		return;

	s += MotifsSeq[3];	// D
	s += MotifsSeq[8];	// d
	s += MotifsSeq[16];	// G
	s += MotifsSeq[34];	// G/S (MAY is RT)
	s += MotifsSeq[35];	// D
	s += MotifsSeq[36];	// D
	}
