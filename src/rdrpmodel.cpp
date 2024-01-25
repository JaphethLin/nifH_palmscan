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
	while (ReadLineStdioFile(f, Line))	// ���ж�ȡ�ļ�������Ϊ�ջ��� ��#�� ��ͷ��������
		{
		if (StartsWith(Line, "#"))	
			continue;
		if (Line.empty())
			continue;

		Split(Line, Fields, '\t');		// �� \t �ָ�ÿ���ֶ�
		asserta(SIZE(Fields) == 3);		
		asserta(SIZE(Fields[0]) == 1);
		char MotifLetter = toupper(Fields[0][0]);	// ��ȡ Motif �ֶ�
		const string &GroupName = Fields[1];		// ��ȡ group �ֶ�
		const string &FastaFileName = Fields[2];	// ��ȡ �ļ�·�� �ֶ�

		bool Found = false;							// �ж�����ָ���� PSSM ������ Motif �Ƿ��ѱ����룬��û�б���������룬���Ѿ�������������
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

		m_PSSMMotifLetters.push_back(MotifLetter);		// �������ָ���� PSSM ��Ӧ Motif ���ַ�
		m_PSSMGroupNames.push_back(GroupName);			// �������ָ���� PSSM ��Ӧ Group ���ַ�
		m_PSSMs.resize(PSSMCount+1);				
		m_PSSMs[PSSMCount].FromFasta(FastaFileName);	// ����ǰ����õ� PSSM �ļ��ж�ȡģ����Ϣ������ PSSMs ������
		++PSSMCount;
		}
	CloseStdioFile(f);

	SetPIV();
	}

void RdRpModel::SetPIV()	// ���� PSSM ģ�ͺ� Motif �� Group ֮��Ĺ�ϵ
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
		SearchAA1_Gapped(PSSMIndex, Seq, Hit);			// Gapped ����
	else
		SearchAA1_Ungapped(PSSMIndex, Seq, Hit);		// Ungapped ����
	
	extern FILE *g_fPSSMAln;
	if (g_fPSSMAln == 0)	// ������ļ��Ҳ����������
		return;
	if (Hit.m_Score < 0)	// ���ȶԵ÷�С�� 0 �������ǰ�ȶ�
		return;

	FILE *f = g_fPSSMAln;	// ����ļ����

	const string &GroupName = m_PSSMGroupNames[PSSMIndex].c_str();	// ��õ�ǰ PSSM �� Group ����
	char MotifLetter = m_PSSMMotifLetters[PSSMIndex];				// ��õ�ǰ PSSM �� Motif ��ǩ

	uint GroupIndex = UINT_MAX;
	uint MotifIndex = UINT_MAX;
	for (uint g = 0; g < SIZE(m_GroupNames); ++g)					// ��ȡ��ǰ PSSMIndex ��Ӧ�� GroupIndex �� MotifIndex
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
	/// Q = Query ��������Ƭ�εİ���������
	/// P = �� PSSM ģ�͵ı�������
	/// A = �ȶԵ� Score ���� '|'(>=1)  '+'(>=0.5)  '.'(>0)  " "(<=0)	
	/// QLo = Query ��������Ƭ�εĵ�һ���������λ����  
	/// QLo = Query ��������Ƭ�ε����һ���������λ����
	fprintf(f, "%s\n", Q.c_str());
	fprintf(f, "%s\n", A.c_str());
	fprintf(f, "%s\n", P.c_str());
	}

// �� Gap ���������м���
void RdRpModel::SearchAA1_Ungapped(uint PSSMIndex, const string &Seq,
  RPHit &Hit)
	{
	CheckThread();
	const uint L = SIZE(Seq);						// L -> ���г���
	const PSSM &P = GetPSSM(PSSMIndex);				// P -> PSSM
	PSSMHitAA HitA;
	GetTopHitAA(P, -9, "_notused_", Seq.c_str(), L, HitA);	 // ���Ƭ����ѵ÷ֺ���ʼλ��
	Hit.SetUngapped(P, HitA.AAPos, HitA.Score);		// ͨ�� HitA (PSSMHitAA �ṹ) �� Hit (RPHit ��) ��ֵ
	}

// �� Gap ���������м���
void RdRpModel::SearchAA1_Gapped(uint PSSMIndex, const string &Seq,
  RPHit &Hit)
	{
	CheckThread();
	// ���� TrimQueryTermGaps ����
	void TrimQueryTermGaps(const char *Path,
	  unsigned &QLo, unsigned &QHi,
	  unsigned &ColLo, unsigned &ColHi);

	// ���� Viterbi_PSSM ����
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

// ��ȡ Motif2 �ȶ���λ�� Motif1 �ȶ�ĩβ�ľ���
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

// ��ȡ Motif1 �ȶ���λ�� Motif2/3 �ȶ�ĩβ�ľ���
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

// ��ø� Motif �� Query �������ҵ���Ƭ�εİ���������Ƭ�� (���� RPHit ���� Path �е� M �� D )
void RdRpModel::GetXSeq(char MotifLetter, string &XSeq) const
	{
	XSeq = "";
	uint GroupIndex = GetTopGroup();	// ��� Group
	if (GroupIndex == UINT_MAX)
		{
		XSeq = "";
		return;
		}

	const uint QL = SIZE(m_QuerySeq);
	const uint PL = GetPSSMLength(MotifLetter, GroupIndex);
	const RPHit *Hit = GetHit(MotifLetter, GroupIndex);
	if (Hit == 0)				// Hit ����Ϊ��
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

	const uint QL = SIZE(m_QuerySeq);							// QL����ѯ���еĳ���
	const uint PL = GetPSSMLength(MotifIndex, GroupIndex);		// PL��PSSM ����
	const RPHit *Hit = GetHit(MotifIndex, GroupIndex);			// ��� RPHit ����
	const PSSM *P = GetPSSM(MotifIndex, GroupIndex);			// PSSM ����
	if (Hit == 0)
		return;
	uint QPos = Hit->m_QPos;		// ����������ʼλ��
	QLo = QPos + 1;
	uint PPos = 0;					// PSSM ��һ��λ��

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
			char q = m_QuerySeq[QPos];				// ��ѯ��������Ƭ���� QPos λ�õİ�����
			float Score = P->GetScore1(PPos, q);	// QPos λ�õİ������� PSSM �� PPos λ���ϵ� Score
			char p = P->GetConsChar(PPos);			// PSSM �� PPos λ���ϵı��ذ����� 

			if (Score >= 1)							// AnnotRow : " ||++.++|||+|| " 
				AnnotRow += '|';					// [ >=1 : | ]  [ >=0.5 : + ] [ >0 : . ] [ <=0 : " " ]
			else if (Score >= 0.5)
				AnnotRow += '+';
			else if (Score > 0)
				AnnotRow += ".";
			else
				AnnotRow += " ";

			QRow += q;		// QRow : Query ���а���������
			PRow += p;		// PRow : PSSM ���ذ���������

			++QPos;			// QPos ����һλ
			++PPos;			// PPos ����һλ
			break;
			}

		case 'D':	// Delete //	PSSM ����� Query ���У��ڸ�λ������һ��ɾ��
			{
			asserta(QPos < QL);
			char q = m_QuerySeq[QPos++];

			QRow += q;
			PRow += '-';
			AnnotRow += " ";
			break;
			}

		case 'I':	// Insert //	PSSM ����� Query ���У��ڸ�λ������һ������
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

// ��������������
void RdRpModel::SearchAA(const string &QueryLabel, const string &QuerySeq)
	{
	CheckThread();

	m_QueryLabel = QueryLabel;		// ���������еı�ǩ
	m_QuerySeq = QuerySeq;			// ������������

	const uint PSSMCount = GetPSSMCount();		// PSSM ģ�͵�����
	m_HitsU.resize(PSSMCount);					// ������Gap���ж���Ϊģ��������Ӧ��С
	m_HitsG.resize(PSSMCount);					// ������Gap���ж���Ϊģ��������Ӧ��С

	for (uint PSSMIndex = 0; PSSMIndex < PSSMCount; ++PSSMIndex)			// ѭ��ʹ�ø��� PSSM ģ�ͼ�������
		{
		RPHit &Hit = (m_Gapped ? m_HitsG[PSSMIndex] : m_HitsU[PSSMIndex]);	// �ж��� Gapped ���� Ungapped ����
		SearchAA1(PSSMIndex, m_QuerySeq, Hit);								// ��ȡ Query ��������� PSSM ģ�ͱȶԵ����н�� Hit
		}

	CheckThread();
	SetResult();		// ���ļ��㺯�������� m_Result��RPResult �ࣩ 
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

// ���ĳ Group ���ܵ÷֣����� Group ������ Motif �� PSSM �÷�֮��
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

// ������ Group �� Motif C �ıȶ� Score
float RdRpModel::GetCScore() const
	{
	uint GroupIndex = GetTopGroup();	// ������ Group
	if (GroupIndex == UINT_MAX)			// ����� Group ������(��С��0)���򷵻� -999
		return -999;
	float Min = 999;
	uint MotifIndex = GetMotifIndex('C');
	const RPHit *HitG = GetHit(MotifIndex, GroupIndex);		// ��� Motif C �ıȶԽ��
	if (HitG == 0)						// �� Motif C �ıȶԽ�������ڣ��򷵻� -999
		return -999;					
	float Score = HitG->m_Score;		// ��� Motif C �ıȶ� Score
	return Score;
	}

// Lim edited //
// ������ Group �� Motif B �ıȶ� Score
float RdRpModel::GetBScore() const
{
	uint GroupIndex = GetTopGroup();	// ������ Group
	if (GroupIndex == UINT_MAX)			// ����� Group ������(��С��0)���򷵻� -999
		return -999;
	float Min = 999;
	uint MotifIndex = GetMotifIndex('B');
	const RPHit* HitG = GetHit(MotifIndex, GroupIndex);		// ��� Motif B �ıȶԽ��
	if (HitG == 0)						// �� Motif B �ıȶԽ�������ڣ��򷵻� -999
		return -999;
	float Score = HitG->m_Score;		// ��� Motif B �ıȶ� Score
	return Score;
}

// Lim edited //
// ������ Group �� Motif A �ıȶ� Score
float RdRpModel::GetAScore() const
{
	uint GroupIndex = GetTopGroup();	// ������ Group
	if (GroupIndex == UINT_MAX)			// ����� Group ������(��С��0)���򷵻� -999
		return -999;
	float Min = 999;
	uint MotifIndex = GetMotifIndex('A');
	const RPHit* HitG = GetHit(MotifIndex, GroupIndex);		// ��� Motif A �ıȶԽ��
	if (HitG == 0)						// �� Motif A �ıȶԽ�������ڣ��򷵻� -999
		return -999;
	float Score = HitG->m_Score;		// ��� Motif A �ıȶ� Score
	return Score;
}


// ������ Group �� Motif �ȶԵ���С Score
float RdRpModel::GetMinScore() const
	{
	uint GroupIndex = GetTopGroup();
	if (GroupIndex == UINT_MAX)			// ����� Group ������(��С��0)���򷵻� -999
		return -999;
	float Min = 999;
	for (uint MotifIndex = 0; MotifIndex < GetMotifCount(); ++MotifIndex)
		{
		const RPHit *HitG = GetHit(MotifIndex, GroupIndex);
		if (HitG == 0)					// ����һ�� Motif �ıȶԽ�������ڣ��򷵻� -1
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

// ���� m_HitsU[PSSMIndex] ����RPHit �ࣩ
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

// �����ߵ÷ֵ� Group �� Index �� Score������ Score <= 0 �򷵻� UINT_MAX
void RdRpModel::GetTopGroup(uint &TopIx, float &TopScore) const
	{
	TopIx = UINT_MAX;
	TopScore = 0;

	for (uint GroupIndex = 0; GroupIndex < GetGroupCount(); ++GroupIndex)
		{
		float Score = GetTotalScore(GroupIndex);	// ���� group �ܵ÷�
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

	uint QAL = HitA->GetQuerySegLength();  // Query �����и� Hit �ĳ���
	uint QBL = HitB->GetQuerySegLength();
	uint QCL = HitC->GetQuerySegLength();

	uint HiA = PosA + QAL - 1;
	uint HiB = PosB + QBL - 1;
	uint HiC = PosC + QCL - 1;

	Start = min(min(PosA, PosB), PosC);		// ���� motif �������ϵ���С��ʼλ��
	End = max(max(HiA, HiB), HiC);			// ���� motif �������ϵ�������λ��

	const uint QL = SIZE(m_QuerySeq);		// Query ���е��ܳ���
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
	return true;	// ���еĽ���λ�� < �������һ��������λ�� ���� true
	}

// ��� ABC motif ��˳���ַ���
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
	if (HitA == 0 || HitB == 0 || HitC == 0)  // δ�ҵ� PSSMIndex ��Ӧ�� RPHit ����
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

// ͨ�� Motif �ı�ǩ��� Motif ����������δ�ҵ����򷵻� UINT_MAX
uint RdRpModel::GetMotifIndex(char MotifLetter) const
	{
	for (uint MotifIndex = 0; MotifIndex < GetMotifCount(); ++MotifIndex)
		{
		if (m_MotifLetters[MotifIndex] == MotifLetter)
			return MotifIndex;
		}
	return UINT_MAX;
	}

// ͨ�� Group �����ƻ�� Group ����������δ�ҵ����򷵻� 0
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

// ���ݵ�һ�� Motif �ȶԵ���ʼ�����һ�� Motif �ȶԵĽ���λ�� �ü� Query ����
void RdRpModel::GetTrimmedSeq(string &Seq) const
	{
	Seq.clear();

	uint TopGroup = GetTopGroup();		// ��� Group
	if (TopGroup == UINT_MAX)
		return;

	uint Start, End;
	bool Ok = GetABCRange(Start, End);	// ��ȡ Query �����ϵ���ʼ�ͽ���λ��
	if (!Ok)
		return;

	const uint QL = SIZE(m_QuerySeq);	// Query ���г���
	if (End >= QL)
		{
		Die("End %u >= QL %u", End, QL);
		return;
		}

	for (uint Pos = Start; Pos <= End; ++Pos)	// �ü� Query ����
		Seq += m_QuerySeq[Pos];
	}

/// MotifA ƥ�䵽��Ƭ�� + "xxx" + MotifB ƥ�䵽��Ƭ�� + "xxx" + MotifC ƥ�䵽��Ƭ��
void RdRpModel::GetMotifsSeq(string &s) const
	{
	s.clear();

	string A, B, C;
	GetXSeq('A', A);								// �ҵ� "A" ��ʾ�� Motif ���� Query ������ƥ�䵽�İ�����Ƭ��
	GetXSeq('B', B);
	GetXSeq('C', C);
	if (A.empty() || B.empty() || C.empty())		// ��һ Motif ������ Group �� Hit �Ҳ��� 
		return;

	s = A + "xxx" + B + "xxx" + C;					// �������� Motif ƥ�䵽��Ƭ��
	}

/// 
/// �� GroupName Ϊ Duplorna��Kitrino��Lenua �� Pisu����� Motif A �ȶԵõ��� Query ����Ƭ�Σ��ڵ� 8 λ����һ�� "-", ��ɾ���� 12 λ
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
	if (A.empty() || B.empty() || C.empty())		// ��һ Motif ������ Group �� Hit �Ҳ��� 
		return;

	// asserta(SIZE(A) == 12);		// ���ȶ��ϵ� Motif A ������������ Group ���ǱȶԵõ��� Query ����Ƭ�γ��ȱ�Ϊ 12
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
	Rows.resize(4);				// �����ĸ��ַ���

	string &TopLine = Rows[0];	
	string &QLine = Rows[1];
	string &ALine = Rows[2];
	string &PLine = Rows[3];

	uint TopGroup;								// ��� Group
	float TopScore;								// ��� Group �� Score
	GetTopGroup(TopGroup, TopScore);
	if (TopGroup == UINT_MAX)		// �� TopGroup == UINT_MAX ���ʾû��ƥ�䵽�� Score > 0 �� group
		{
		Rows.clear();
		return;
		}

	string TopGroupName;						// ��� Group ������
	GetGroupName(TopGroup, TopGroupName);	

	string XXX;									// ������ Motif ��˳�� 
	GetABCOrder(XXX);							// ��˳����ȷΪ "ABC" ������Ϊ ""

	bool AnyMotifs = false;
	char PrevMotifLetter = 0;
	for (uint k = 0; k < 3; ++k)
		{
		char MotifLetter = XXX[k];
		uint MotifIndex = GetMotifIndex(MotifLetter);	// ���� MotifLetter �ҵ���Ӧ����������˳����ȷ MotifLetter Ϊ "" , MotifIndex Ϊ UINT_MAX
		if (MotifIndex == UINT_MAX)						// δ�ҵ�
			continue;
		AnyMotifs = true;

		uint QLo, QHi;							// Query �������е� ��ʼλ�ã�QLo���� ����λ�ã�QHi��
		string QRow, PRow, ARow;				// Query ���еİ��������У�QRow����PSSM ���еİ��������У�PRow����Query ���е� Score
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
		else			// �����һ�����ֵ� Motif 
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


	if (!AnyMotifs)		// ���Ҳ����κ� Motif ����� Rows ��������ǰ����
		{
		Rows.clear();
		return;
		}

	uint Start;			// ��һ�� Motif ��ʼλ��
	uint End;			// ���һ�� Motif ����λ��
	bool RangeOk = GetABCRange(Start, End);
	if (RangeOk)
		Psa(QLine, "  [%d]", int(End) - int(Start) + 1);

	///	'    'A:4-19(2.4)'  ''   'B:27-47(1.2)'      ''   'B:56-67(9.0)'
	/// '    'ABACCDJSHMALDSP <8> LKSPAIEMNDSJKASJJDNF <9> IFDUBAIFBIA''' '[63]
	/// '    '|||||+||||+||.|'   '|||||||+++|||.||+|||'   '|||||||||||''
	/// '    'ABACCGJSHMARDSP'   'LKSPAIEMFASJKASJTTTF'   'IFDUBAIFBIA''
	}


/// ���Ѿ��ȶԵ��� Query ����Ƭ������ SuperMotif Ƭ��
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
