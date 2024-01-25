#include "myutils.h"
#include "rdrpmodel.h"
#include "heuristics.h"

static bool MatchMotif1(char q, char t)
	{
	if (q == t || t == 'x' || t == '.')
		return true;
	return false;
	}

// Motif �� Pattern ��һһ��Ӧ�򷵻� FALSE��Pattern �е� 'x' �� '.' ���⣩ 
static bool MatchMotif(const string &Motif, const string &Pattern)
	{
	const uint L = SIZE(Motif);
	asserta(SIZE(Pattern) == L);
	for (uint i = 0; i < L; ++i)
		if (!MatchMotif1(Motif[i], Pattern[i]))
			return false;
	return true;
	}

void RdRpModel::SetResult()
	{
	CheckThread();
	m_Result.Clear();
	m_Result.m_QueryLabel = m_QueryLabel;
	m_Result.m_Comments.clear();
	m_Result.m_Gene = "unclassified";
	m_Result.m_HiConf = false;
	m_Result.m_FinalScore = 0;

	GetAlnRows(m_Result.m_Aln);
	/// m_Aln:[TopLine][QLine][ALine][PLine]
	///	[0] '    'A:4-19(2.4)'  ''   'B:27-47(1.2)'      ''   'B:56-67(9.0)'
	/// [1] '    'ABACCDJSHMALDSP <8> LKSPAIEMNDSJKASJJDNF <9> IFDUBAIFBIA''' '[63]
	/// [2] '    '|||||+||||+||.|'   '|||||||+++|||.||+|||'   '|||||||||||''
	/// [3] '    'ABACCGJSHMARDSP'   'LKSPAIEMFASJKASJTTTF'   'IFDUBAIFBIA''
	
	//GetFullAln(m_Result.m_FullAln);
	GetMotifsSeq(m_Result.m_MotifsSeq);
	/// MotifA ƥ�䵽��Ƭ�� + "xxx" + MotifB ƥ�䵽��Ƭ�� + "xxx" + MotifC ƥ�䵽��Ƭ��
	
	GetMotifsSeq2(m_Result.m_MotifsSeq2);
	/// �� GroupName Ϊ Duplorna��Kitrino��Lenua �� Pisu����� Motif A �ȶԵõ��� Query ����Ƭ�Σ��ڵ� 8 λ����һ�� "-", ��ɾ���� 12 λ
	
	GetSuperMotif(m_Result.m_MotifsSeq, m_Result.m_SuperMotif);
	/// ���Ѿ��ȶԵ��� Query ����Ƭ������ SuperMotif Ƭ��
	
	GetTrimmedSeq(m_Result.m_TrimmedSeq);
	/// ���ݵ�һ�� Motif �ȶԵ���ʼ�����һ�� Motif �ȶԵĽ���λ�� �ü� Query ����

	bool AllowHigh = true;	// ���ܵ÷�����ͨ��
	
	// ���� Query �����е���ֹ�Ӻ�X����
	uint XCount = 0;
	for (uint i = 0; i < SIZE(m_Result.m_TrimmedSeq); ++i)
		{
		char c = m_Result.m_TrimmedSeq[i];
		if (c == '*')								// �� Query ���бȶ��ϵ�Ƭ���а��� * ����ζ�ź�����ֹ�� //
			{
			m_Result.m_Comments += "stop-codon.";	// ������� stop-codon.
			m_Result.m_FinalScore = 0.0;			// ���յ÷�Ϊ 0
			return;
			}

		if (c == 'X')								// ��¼ Query ���бȶ��ϵ�Ƭ�� X ������ (δ֪�����ⰱ����)
			++XCount;
		}

	if (XCount > m_MaxX)							// �� Query ���бȶ��ϵ�Ƭ���а��� > m_MaxX �� X ���� too-many-Xs. ͬʱ���յ÷�Ϊ 0 �����������ν������
		{
		m_Result.m_Comments += "too-many-Xs.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	// ���� Query �����ܳ���
	m_Result.m_QL = SIZE(m_QuerySeq);				// m_QL����ȡ Query ���г���
	if (m_Result.m_QL < MIN_SEG)					// �� Query ���еĳ��� < MIN_SEG �����ж�Ϊ���г��ȹ��� ����� query-too-short. ͬʱ���յ÷�Ϊ 0 �����������ν������
		{
		m_Result.m_Comments += "query-too-short.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
	if (m_Result.m_QL < SHORT_QUERY)				// �� Query ���еĳ��� < SHORT_QUERY �����ж�Ϊ���г��Ƚ϶� ����� short-query. ��������
		m_Result.m_Comments += "short-query.";

	// ���� Query �����Ƿ�ȶԲ��� motifs
	uint TopGroup;
	float TopScore;
	GetTopGroup(TopGroup, TopScore);				
	if (TopGroup == UINT_MAX)						// �� Query ������� Motifs �ıȶԾ������� 0 ����� motifs-not-found. ͬʱ���յ÷�Ϊ 0 �����������ν������
		{
		m_Result.m_Comments += "motifs-not-found.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	// ���� Query ���бȶ� Motifs ���ܷ��Ƿ���ϸ�����Ҫ��
	GetGroupName(TopGroup, m_Result.m_TopGroupName);	//  m_TopGroupName <- ��� Group ������
	m_Result.m_PSSMTotalScore = TopScore;				//	m_PSSMTotalScore <- ��� Group �� Score
	m_Result.m_FinalScore = TopScore;					//	m_FinalScore <- ��� Group �� Score
	if (TopScore >= MIN_HICONF)						// ����� Group ���ܷ� >= MIN_HICONF, ��� high-PSSM-score.
		m_Result.m_Comments += "good-PSSM-score.";

	if (StartsWith(m_Result.m_TopGroupName, "nifH"))
		m_Result.m_Gene = "nifH";

//	if (StartsWith(m_Result.m_TopGroupName, "RT"))		// �Ƚ� m_TopGroupName ��ͷ�ַ��Ƿ�� "RT" �ַ������, ������� m_Gene Ϊ RT������Ϊ RdRP
//		m_Result.m_Gene = "RT";
//	else
//		m_Result.m_Gene = "RdRP";						// m_Gene <- RdRP or RT



	const RPHit *ptrHit1 = 0;
	const RPHit *ptrHit2 = 0;
	const RPHit *ptrHit3 = 0;
	// Motifs ˳�򲻶ԣ�0��
	GetOrderedHits(m_Result.m_XXX, ptrHit1, ptrHit2, ptrHit3);	// m_XXX �ַ���Ϊ Motifs ��˳�� , 'ABC' or '' or 'CAB' (Lim �޸�ǰ)		// m_XXX �ַ���Ϊ Motifs ��˳�� , 'ABC' (Lim �޸ĺ�)
	if (m_Result.m_XXX.empty())						// �� m_XXX �ַ���Ϊ '' , ��˵�� Motif �ȶ��ϵ�˳��Ϊ 'ABC' �� 'CAB' (Lim �޸�ǰ)	// �� m_XXX �ַ���Ϊ '' , ��˵�� Motif �ȶ��ϵ�˳��Ϊ 'ABC' (Lim �޸ĺ�)
		{											//		��� motifs-out-of-order. , ͬʱ���յ÷�Ϊ 0 , ���������ν������
		AllowHigh = false;
		m_Result.m_Comments += "motifs-out-of-order.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

/*	Lim deleted
	// Motifs �ȶԵ� RT , ��˳��Ϊ CAB , 0��
	bool Permuted = (m_Result.m_XXX != "ABC");		// �� m_XXX �ַ���Ϊ 'CAB' �� m_Gene Ϊ 'RT' , 
	if (Permuted && m_Result.m_Gene == "RT")		//		��� reject-permuted-RT. , ͬʱ���յ÷�Ϊ 0 , ���������ν������
		{
		m_Result.m_Comments += "reject-permuted-RT.";		// �ܾ� RT
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/

/*	Lim deleted
	// �ܾ��ȶ�˳��Ϊ CAB ������, 0��
	bool Permuted = (m_Result.m_XXX != "ABC");
	if (Permuted)
		{
		m_Result.m_Comments += "reject-permuted-nonABC.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/

/*	Lim deleted
	// Motif C �÷ֹ��ͣ�0��
	float CScore = GetCScore();						// �� Motif C �� Score < MIN_C_SCORE��
	if (CScore < MIN_C_SCORE)						//		��� reject-low-C-score. , ͬʱ���յ÷�Ϊ 0 , ���������ν������
		{
		AllowHigh = false;
		m_Result.m_Comments += "reject-low-C-score.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/

/*	Lim deleted
	// Motif ����͵÷ֹ�����˳��Ϊ ABC��0��
	bool Permuted = (m_Result.m_XXX != "ABC");
	float MinScore = GetMinScore();					// ����� Group �з�ֵ��͵� Motif Score < MIN_PSSM_SCORE_PERMUTED��
	m_Result.m_PSSMMinScore = MinScore;				//		��� permuted-low-PSSM-score. , ͬʱ���յ÷�Ϊ 0 , ���������ν������
	if (Permuted && MinScore < MIN_PSSM_SCORE_PERMUTED)
		{
		AllowHigh = false;
		m_Result.m_Comments += "permuted-low-PSSM-score.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/


	// Motif ����͵÷ֹ��ͣ��۷�
	float MinScore = GetMinScore();
	m_Result.m_PSSMMinScore = MinScore;
	if (MinScore < MIN_PSSM_SCORE)					// ����� Group �з�ֵ��͵� Motif Score < MIN_PSSM_SCORE��
		{											//		��� penalty-low-PSSM-score. , ͬʱ���յ÷ּ�ȥ MIN_PSSM_PENALTY , ���������ν������
		AllowHigh = false;	// ���ܵ÷�����ͨ��
		m_Result.m_Comments += "penalty-low-PSSM-score.";
		m_Result.m_FinalScore -= MIN_PSSM_PENALTY;
		}

	asserta(ptrHit1 != 0 && ptrHit2 != 0 && ptrHit3 != 0);	// ���˶��ԣ����� Motif �ıȶԽ��������
	m_Result.m_StartPos = ptrHit1->m_QPos;					// m_StartPos <- ��һ���ȶ��ϵ� Motif �� Query ���е���ʼλ��

	// �� Motif �ص���0��
	m_Result.m_Dist12 = GetDist(ptrHit1, ptrHit2);			// ��ȡ Motif2 �ȶ���λ�� Motif1 �ȶ�ĩβ�ľ���
	m_Result.m_Dist23 = GetDist(ptrHit2, ptrHit3);			// ��ȡ Motif3 �ȶ���λ�� Motif2 �ȶ�ĩβ�ľ���
	if (m_Result.m_Dist12 <= 0 || m_Result.m_Dist23 <= 0)	// �������� Motif ���� <= 0 , ���ص�  
		{													//		��� reject-overlapping-motifs. , ͬʱ���յ÷�Ϊ 0 , ���������ν������
		AllowHigh = false;
		m_Result.m_Comments += "reject-overlapping-motifs.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	// ���� Motifs �Ŀ��̫��0��
	uint Span = GetSpan(ptrHit1, ptrHit3);				// ��ȡ Motif1 �ȶ���λ�� Motif3 �ȶ�ĩβ�Ŀ��
	m_Result.m_SegLength = Span;						// m_SegLength <- Motif1 �ȶ���λ�� Motif3 �ȶ�ĩβ�Ŀ��
	if (Span < MIN_SEG || Span > MAX_SEG)			// �� Motifs ��� < MIN_SEG �� ��� > MAX_SEG ,
		{											// ��� bad-segment-length. , ͬʱ���յ÷ּ�ȥ MIN_PSSM_PENALTY , ���������ν������
		AllowHigh = false;
		m_Result.m_Comments += "bad-segment-length.";
		m_Result.m_FinalScore = 0.0;
		return;
		}


	if (m_Result.m_Gene == "nifH")						// �� m_Gene Ϊ nifH Ϊǰ��
		{

/*	Lim deleted
		const string &Super = m_Result.m_SuperMotif;	
		if (m_Result.m_XXX == "ABC" && !Super.empty())	// Motifs ˳��Ϊ ABC , �ҿ����ҵ� SuperMotif ����
			{
			const string &GDD = Super.substr(3, 3);		// GDD <- SuperMotif ���еĵ� 4,5,6 λ������
			if (GDD != "GDD" && GDD != "SDD" && GDD != "GDN")	// �� SuperMotif ���еĵ� 4,5,6 λ�����᲻Ϊ "GDD" / "SDD" / "GDN" ,
				{												// ���� penalty-missing-GDD/SDD/GDN. , ���յ÷ּ�ȥ PENALTY_NOT_GDD_SDD_GDN
				m_Result.m_Comments += "penalty-missing-GDD/SDD/GDN.";
				m_Result.m_FinalScore -= PENALTY_NOT_GDD_SDD_GDN;
				}
			else if (MatchMotif(Super, "DDGGDD"))				// �� SuperMotif ����Ϊ DDGGDD,
				{												// ���� reward-DDGGDD. , ���յ÷ּ��� REWARD_DDGGDD
				m_Result.m_Comments += "reward-DDGGDD.";
				m_Result.m_FinalScore += REWARD_DDGGDD;
				}
			else if (MatchMotif(Super, "DDGSDD"))				// �� SuperMotif ����Ϊ DDGSDD,
				{												// ���� reward-DDGSDD. , ���յ÷ּ��� REWARD_DDGSDD
				m_Result.m_Comments += "reward-DDGSDD.";
				m_Result.m_FinalScore += REWARD_DDGSDD;
				}
			else if (MatchMotif(Super, "DNMSDD"))				// �� SuperMotif ����Ϊ DNMSDD,
			{													// ���� reward-DNMSDD. , ���յ÷ּ��� REWARD_DNMSDD
				m_Result.m_Comments += "reward-DNMSDD.";
				m_Result.m_FinalScore += REWARD_DNMSDD;
				}
			else if (MatchMotif(Super, "DNxGDN"))				// �� SuperMotif ����Ϊ DNxGDN,
				{												// ���� reward-DNxGDN. , ���յ÷ּ��� REWARD_DNxGDN
				m_Result.m_Comments += "reward-DNxGDN.";
				m_Result.m_FinalScore += REWARD_DNxGDN;
				}
			}
*/																// [ ���� 1 ] ( MIN_SEG1 [ ���� 2 ] ( MIN_SEG2 [ ���� 3 ] MAX_SEG2 ) [ ���� 2 ] MAX_SEG1) [ ���� 1 ]
		if (Span < MIN_SEG1 || Span > MAX_SEG1)					// �� Motifs �Ŀ���� ( MIN_SEG1 , MAX_SEG1 ) ֮�� (����1)
			{													// ���� segment-length-penalty1. , ���յ÷ּ�ȥ SEG_LENGTH_PENALTY1
			m_Result.m_Comments += "segment-length-penalty1.";
			m_Result.m_FinalScore -= SEG_LENGTH_PENALTY1;
			}
		else if (Span < MIN_SEG2 || Span > MAX_SEG2)			// �� Motifs �Ŀ���� ( MIN_SEG2 , MAX_SEG2 ) ֮�� (����2)
			{													// ���� segment-length-penalty2. , ���յ÷ּ�ȥ SEG_LENGTH_PENALTY2
			m_Result.m_Comments += "segment-length-penalty2.";
			m_Result.m_FinalScore -= SEG_LENGTH_PENALTY2;
			}
		else													// �� Motifs �Ŀ�� >= MIN_SEG2 �� <= MAX_SEG2 (����3)
			m_Result.m_Comments += "good-segment-length.";		// ���� good-segment-length. 

/*	Lim deleted
		if (m_Result.m_XXX == "CAB")							// ���ȶ�Ϊ RdRP , �� Motif ˳��Ϊ CAB ,
			{													// ���� permuted-penalty. , ���յ÷ּ�ȥ PERMUTED_PENALTY
			m_Result.m_Comments += "permuted-penalty.";
			m_Result.m_FinalScore -= PERMUTED_PENALTY;
			}
*/		
		if (m_Result.m_XXX == "ABC")	//	Lim added
			{
			int AB = m_Result.m_Dist12;
			int BC = m_Result.m_Dist23;
			if (AB < MIN_AB1 || AB > MAX_AB1)					// �� Motif 1 �� 2 ֮��ľ��� < MIN_AB1 �� > MAX_AB1
				{												// ���� AB-dist-penalty. , ���յ÷ּ�ȥ AB_PENALTY1
				m_Result.m_Comments += "AB-dist-penalty.";
				m_Result.m_FinalScore -= AB_PENALTY1;
				}
			if (BC < MIN_BC1 || BC > MAX_BC1)					// �� Motif 2 �� 3 ֮��ľ��� < MIN_BC1 �� > MAX_BC1
				{												// ���� BC-dist-penalty. , ���յ÷ּ�ȥ BC_PENALTY1
				m_Result.m_Comments += "BC-dist-penalty.";
				m_Result.m_FinalScore -= BC_PENALTY1;
				}
			}
		} // end if RdRP

	// ���յ÷ֹ��� (<0) , 0��
	if (m_Result.m_FinalScore < 0)								// �� m_FinalScore ���յ÷� < 0 (��Ҫ���ڷ��� , ��ʼ��ֵΪ 0)
		{														// �� m_FinalScore <- 0 , ��� m_Gene Ϊ unclassified , ��� m_Comments Ϊ negative-score.
		m_Result.m_FinalScore = 0;								// ��������
		m_Result.m_Gene = "unclassified";
		m_Result.m_Comments += "negative-score.";
		return;
		}

	if (!AllowHigh && m_Result.m_FinalScore >= MIN_HICONF)		// ��Ҫ������ Motif �÷ֹ��� , ���ܵ÷ָߵ���� 
		m_Result.m_FinalScore = MIN_LOCONF;						// ���÷ָ�ֵΪ MIN_LOCONF

	m_Result.m_HiConf = (m_Result.m_FinalScore >= myHICONF);	// �� AllowHight ��Ϊ True �� m_FinalScore >= Highly Confirm Score ,
	}															// m_HiConf <- True , ����Ϊ False

///
/// score=3.1
/// query=Label 
/// gene=RdRP
/// confidence=0.9
/// comments="good-segment-length.BC-dist-penalty."
/// order=ABC
/// qlen=219
/// pp_start=12
/// pp_end=201
/// pp_length=189
/// pssm_total_score=1.3
/// pssm_min_score=0.3
/// order=ABC
/// v1_length=34
/// v2_length=28
/// group=Lena
/// super=GDNAGDMMDG
/// motifs=
/// 
void RPResult::ToFEVStr(string &s) const
	{
	Ps(s, "score=%.1f", m_FinalScore);
	Psa(s, "\tquery=%s", m_QueryLabel.c_str());
	Psa(s, "\tgene=%s", m_Gene.c_str());
	// Lim edited
	if (m_Gene == "nifH" || m_Gene == "RdRP" || m_Gene == "RT")
		Psa(s, "\tconfidence=%s", m_HiConf ? "high" : "low");
	Psa(s, "\tcomments=%s", m_Comments.c_str());
	if (!m_XXX.empty())
		Psa(s, "\torder=%s", m_XXX.c_str());
	if (m_QLnt != 0 && m_QLnt != UINT_MAX)
		Psa(s, "\tqlen=%u", m_QLnt);
	else if (m_QL != 0 && m_QL != UINT_MAX)
		Psa(s, "\tqlen=%u", m_QL);
	if (m_StartPosNt != UINT_MAX)
		{
		Psa(s, "\txlat=yes");
		Psa(s, "\tpp_start=%u", m_StartPosNt + 1);
		Psa(s, "\tpp_end=%u", m_StartPosNt + 3*m_SegLength);
		Psa(s, "\tpp_length=%u", 3*m_SegLength);
		}
	else if (m_StartPos != UINT_MAX)
		{
		Psa(s, "\tpp_start=%u", m_StartPos + 1);
		Psa(s, "\tpp_end=%u", m_StartPos + m_SegLength);
		Psa(s, "\tpp_length=%u", m_SegLength);
		}
	if (m_Frame != 0 && m_Frame != -999)
		Psa(s, "\tframe=%+d", m_Frame);
	if (m_PSSMTotalScore > -10)
		Psa(s, "\tpssm_total_score=%.1f", m_PSSMTotalScore);
	if (m_PSSMMinScore > -10)
		Psa(s, "\tpssm_min_score=%.1f", m_PSSMMinScore);
	if (m_XXX == "ABC" || m_XXX == "CAB")
		{
		Psa(s, "\torder=%s", m_XXX.c_str());
		Psa(s, "\tv1_length=%d", m_Dist12);
		Psa(s, "\tv2_length=%d", m_Dist23);
		}
	if (!m_TopGroupName.empty())
		Psa(s, "\tgroup=%s", m_TopGroupName.c_str());
	if (!m_SuperMotif.empty())
		Psa(s, "\tsuper=%s", m_SuperMotif.c_str());
	if (!m_MotifsSeq.empty())
		Psa(s, "\tmotifs=%s", m_MotifsSeq.c_str());
	}

// δ���н���ĳ�ʼ��
void RdRpModel::SetResult_NoNtHit(const string &QueryLabel, uint QLnt,
  RPResult &Result)
	{
	Result.Clear();
	Result.m_QueryLabel = QueryLabel;
	Result.m_QL = 0;
	Result.m_QLnt = QLnt;
	Result.m_Comments.clear();
	Result.m_Gene = "unclassified";
	Result.m_HiConf = false;
	Result.m_Comments = "six-frame-none-found.";
	Result.m_FinalScore = 0;
	}

// δ���н���ĳ�ʼ��
void RdRpModel::SetResult_NoAaHit(const string &QueryLabel, uint QLaa,
  RPResult &Result)
	{
	Result.Clear();
	Result.m_QueryLabel = QueryLabel;
	Result.m_QL = QLaa;
	Result.m_QLnt = UINT_MAX;
	Result.m_Comments.clear();
	Result.m_Gene = "unclassified";
	Result.m_HiConf = false;
	Result.m_Comments = "aa-none-found.";
	Result.m_FinalScore = 0;
	}
