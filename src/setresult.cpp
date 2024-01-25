#include "myutils.h"
#include "rdrpmodel.h"
#include "heuristics.h"

static bool MatchMotif1(char q, char t)
	{
	if (q == t || t == 'x' || t == '.')
		return true;
	return false;
	}

// Motif 和 Pattern 不一一对应则返回 FALSE（Pattern 中的 'x' 和 '.' 除外） 
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
	/// MotifA 匹配到的片段 + "xxx" + MotifB 匹配到的片段 + "xxx" + MotifC 匹配到的片段
	
	GetMotifsSeq2(m_Result.m_MotifsSeq2);
	/// 若 GroupName 为 Duplorna、Kitrino、Lenua 或 Pisu，则对 Motif A 比对得到的 Query 序列片段，在第 8 位插入一个 "-", 并删除第 12 位
	
	GetSuperMotif(m_Result.m_MotifsSeq, m_Result.m_SuperMotif);
	/// 从已经比对到的 Query 序列片段中找 SuperMotif 片段
	
	GetTrimmedSeq(m_Result.m_TrimmedSeq);
	/// 根据第一个 Motif 比对的起始与最后一个 Motif 比对的结束位置 裁剪 Query 序列

	bool AllowHigh = true;	// 高总得分允许通过
	
	// 检验 Query 序列中的终止子和X数量
	uint XCount = 0;
	for (uint i = 0; i < SIZE(m_Result.m_TrimmedSeq); ++i)
		{
		char c = m_Result.m_TrimmedSeq[i];
		if (c == '*')								// 若 Query 序列比对上的片段中包含 * ，意味着含有终止子 //
			{
			m_Result.m_Comments += "stop-codon.";	// 添加评述 stop-codon.
			m_Result.m_FinalScore = 0.0;			// 最终得分为 0
			return;
			}

		if (c == 'X')								// 记录 Query 序列比对上的片段 X 的数量 (未知或任意氨基酸)
			++XCount;
		}

	if (XCount > m_MaxX)							// 若 Query 序列比对上的片段中包含 > m_MaxX 个 X 则标记 too-many-Xs. 同时最终得分为 0 ，并结束本次结果计算
		{
		m_Result.m_Comments += "too-many-Xs.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	// 检验 Query 序列总长度
	m_Result.m_QL = SIZE(m_QuerySeq);				// m_QL：获取 Query 序列长度
	if (m_Result.m_QL < MIN_SEG)					// 若 Query 序列的长度 < MIN_SEG ，则判定为序列长度过短 ，标记 query-too-short. 同时最终得分为 0 ，并结束本次结果计算
		{
		m_Result.m_Comments += "query-too-short.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
	if (m_Result.m_QL < SHORT_QUERY)				// 若 Query 序列的长度 < SHORT_QUERY ，则判定为序列长度较短 ，标记 short-query. 继续计算
		m_Result.m_Comments += "short-query.";

	// 检验 Query 序列是否比对不到 motifs
	uint TopGroup;
	float TopScore;
	GetTopGroup(TopGroup, TopScore);				
	if (TopGroup == UINT_MAX)						// 若 Query 序列与各 Motifs 的比对均不大于 0 ，标记 motifs-not-found. 同时最终得分为 0 ，并结束本次结果计算
		{
		m_Result.m_Comments += "motifs-not-found.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	// 检验 Query 序列比对 Motifs 的总分是否符合高置信要求
	GetGroupName(TopGroup, m_Result.m_TopGroupName);	//  m_TopGroupName <- 最佳 Group 的名称
	m_Result.m_PSSMTotalScore = TopScore;				//	m_PSSMTotalScore <- 最佳 Group 的 Score
	m_Result.m_FinalScore = TopScore;					//	m_FinalScore <- 最佳 Group 的 Score
	if (TopScore >= MIN_HICONF)						// 若最佳 Group 的总分 >= MIN_HICONF, 标记 high-PSSM-score.
		m_Result.m_Comments += "good-PSSM-score.";

	if (StartsWith(m_Result.m_TopGroupName, "nifH"))
		m_Result.m_Gene = "nifH";

//	if (StartsWith(m_Result.m_TopGroupName, "RT"))		// 比较 m_TopGroupName 开头字符是否和 "RT" 字符串相等, 若相等则 m_Gene 为 RT，否则为 RdRP
//		m_Result.m_Gene = "RT";
//	else
//		m_Result.m_Gene = "RdRP";						// m_Gene <- RdRP or RT



	const RPHit *ptrHit1 = 0;
	const RPHit *ptrHit2 = 0;
	const RPHit *ptrHit3 = 0;
	// Motifs 顺序不对，0分
	GetOrderedHits(m_Result.m_XXX, ptrHit1, ptrHit2, ptrHit3);	// m_XXX 字符串为 Motifs 的顺序 , 'ABC' or '' or 'CAB' (Lim 修改前)		// m_XXX 字符串为 Motifs 的顺序 , 'ABC' (Lim 修改后)
	if (m_Result.m_XXX.empty())						// 若 m_XXX 字符串为 '' , 则说明 Motif 比对上的顺序不为 'ABC' 或 'CAB' (Lim 修改前)	// 若 m_XXX 字符串为 '' , 则说明 Motif 比对上的顺序不为 'ABC' (Lim 修改后)
		{											//		标记 motifs-out-of-order. , 同时最终得分为 0 , 并结束本次结果计算
		AllowHigh = false;
		m_Result.m_Comments += "motifs-out-of-order.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

/*	Lim deleted
	// Motifs 比对到 RT , 但顺序为 CAB , 0分
	bool Permuted = (m_Result.m_XXX != "ABC");		// 若 m_XXX 字符串为 'CAB' 且 m_Gene 为 'RT' , 
	if (Permuted && m_Result.m_Gene == "RT")		//		标记 reject-permuted-RT. , 同时最终得分为 0 , 并结束本次结果计算
		{
		m_Result.m_Comments += "reject-permuted-RT.";		// 拒绝 RT
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/

/*	Lim deleted
	// 拒绝比对顺序为 CAB 的序列, 0分
	bool Permuted = (m_Result.m_XXX != "ABC");
	if (Permuted)
		{
		m_Result.m_Comments += "reject-permuted-nonABC.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/

/*	Lim deleted
	// Motif C 得分过低，0分
	float CScore = GetCScore();						// 若 Motif C 的 Score < MIN_C_SCORE，
	if (CScore < MIN_C_SCORE)						//		标记 reject-low-C-score. , 同时最终得分为 0 , 并结束本次结果计算
		{
		AllowHigh = false;
		m_Result.m_Comments += "reject-low-C-score.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/

/*	Lim deleted
	// Motif 的最低得分过低且顺序不为 ABC，0分
	bool Permuted = (m_Result.m_XXX != "ABC");
	float MinScore = GetMinScore();					// 若最佳 Group 中分值最低的 Motif Score < MIN_PSSM_SCORE_PERMUTED，
	m_Result.m_PSSMMinScore = MinScore;				//		标记 permuted-low-PSSM-score. , 同时最终得分为 0 , 并结束本次结果计算
	if (Permuted && MinScore < MIN_PSSM_SCORE_PERMUTED)
		{
		AllowHigh = false;
		m_Result.m_Comments += "permuted-low-PSSM-score.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
*/


	// Motif 的最低得分过低，扣分
	float MinScore = GetMinScore();
	m_Result.m_PSSMMinScore = MinScore;
	if (MinScore < MIN_PSSM_SCORE)					// 若最佳 Group 中分值最低的 Motif Score < MIN_PSSM_SCORE，
		{											//		标记 penalty-low-PSSM-score. , 同时最终得分减去 MIN_PSSM_PENALTY , 并结束本次结果计算
		AllowHigh = false;	// 高总得分允许通过
		m_Result.m_Comments += "penalty-low-PSSM-score.";
		m_Result.m_FinalScore -= MIN_PSSM_PENALTY;
		}

	asserta(ptrHit1 != 0 && ptrHit2 != 0 && ptrHit3 != 0);	// 到此断言：三个 Motif 的比对结果均存在
	m_Result.m_StartPos = ptrHit1->m_QPos;					// m_StartPos <- 第一个比对上的 Motif 在 Query 序列的起始位置

	// 两 Motif 重叠，0分
	m_Result.m_Dist12 = GetDist(ptrHit1, ptrHit2);			// 获取 Motif2 比对首位到 Motif1 比对末尾的距离
	m_Result.m_Dist23 = GetDist(ptrHit2, ptrHit3);			// 获取 Motif3 比对首位到 Motif2 比对末尾的距离
	if (m_Result.m_Dist12 <= 0 || m_Result.m_Dist23 <= 0)	// 若任意两 Motif 距离 <= 0 , 即重叠  
		{													//		标记 reject-overlapping-motifs. , 同时最终得分为 0 , 并结束本次结果计算
		AllowHigh = false;
		m_Result.m_Comments += "reject-overlapping-motifs.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	// 三个 Motifs 的跨度太大，0分
	uint Span = GetSpan(ptrHit1, ptrHit3);				// 获取 Motif1 比对首位到 Motif3 比对末尾的跨度
	m_Result.m_SegLength = Span;						// m_SegLength <- Motif1 比对首位到 Motif3 比对末尾的跨度
	if (Span < MIN_SEG || Span > MAX_SEG)			// 若 Motifs 跨度 < MIN_SEG 或 跨度 > MAX_SEG ,
		{											// 标记 bad-segment-length. , 同时最终得分减去 MIN_PSSM_PENALTY , 并结束本次结果计算
		AllowHigh = false;
		m_Result.m_Comments += "bad-segment-length.";
		m_Result.m_FinalScore = 0.0;
		return;
		}


	if (m_Result.m_Gene == "nifH")						// 以 m_Gene 为 nifH 为前提
		{

/*	Lim deleted
		const string &Super = m_Result.m_SuperMotif;	
		if (m_Result.m_XXX == "ABC" && !Super.empty())	// Motifs 顺序为 ABC , 且可以找到 SuperMotif 序列
			{
			const string &GDD = Super.substr(3, 3);		// GDD <- SuperMotif 序列的第 4,5,6 位氨基酸
			if (GDD != "GDD" && GDD != "SDD" && GDD != "GDN")	// 若 SuperMotif 序列的第 4,5,6 位氨基酸不为 "GDD" / "SDD" / "GDN" ,
				{												// 则标记 penalty-missing-GDD/SDD/GDN. , 最终得分减去 PENALTY_NOT_GDD_SDD_GDN
				m_Result.m_Comments += "penalty-missing-GDD/SDD/GDN.";
				m_Result.m_FinalScore -= PENALTY_NOT_GDD_SDD_GDN;
				}
			else if (MatchMotif(Super, "DDGGDD"))				// 若 SuperMotif 序列为 DDGGDD,
				{												// 则标记 reward-DDGGDD. , 最终得分加上 REWARD_DDGGDD
				m_Result.m_Comments += "reward-DDGGDD.";
				m_Result.m_FinalScore += REWARD_DDGGDD;
				}
			else if (MatchMotif(Super, "DDGSDD"))				// 若 SuperMotif 序列为 DDGSDD,
				{												// 则标记 reward-DDGSDD. , 最终得分加上 REWARD_DDGSDD
				m_Result.m_Comments += "reward-DDGSDD.";
				m_Result.m_FinalScore += REWARD_DDGSDD;
				}
			else if (MatchMotif(Super, "DNMSDD"))				// 若 SuperMotif 序列为 DNMSDD,
			{													// 则标记 reward-DNMSDD. , 最终得分加上 REWARD_DNMSDD
				m_Result.m_Comments += "reward-DNMSDD.";
				m_Result.m_FinalScore += REWARD_DNMSDD;
				}
			else if (MatchMotif(Super, "DNxGDN"))				// 若 SuperMotif 序列为 DNxGDN,
				{												// 则标记 reward-DNxGDN. , 最终得分加上 REWARD_DNxGDN
				m_Result.m_Comments += "reward-DNxGDN.";
				m_Result.m_FinalScore += REWARD_DNxGDN;
				}
			}
*/																// [ 区间 1 ] ( MIN_SEG1 [ 区间 2 ] ( MIN_SEG2 [ 区间 3 ] MAX_SEG2 ) [ 区间 2 ] MAX_SEG1) [ 区间 1 ]
		if (Span < MIN_SEG1 || Span > MAX_SEG1)					// 若 Motifs 的跨度在 ( MIN_SEG1 , MAX_SEG1 ) 之外 (区间1)
			{													// 则标记 segment-length-penalty1. , 最终得分减去 SEG_LENGTH_PENALTY1
			m_Result.m_Comments += "segment-length-penalty1.";
			m_Result.m_FinalScore -= SEG_LENGTH_PENALTY1;
			}
		else if (Span < MIN_SEG2 || Span > MAX_SEG2)			// 若 Motifs 的跨度在 ( MIN_SEG2 , MAX_SEG2 ) 之外 (区间2)
			{													// 则标记 segment-length-penalty2. , 最终得分减去 SEG_LENGTH_PENALTY2
			m_Result.m_Comments += "segment-length-penalty2.";
			m_Result.m_FinalScore -= SEG_LENGTH_PENALTY2;
			}
		else													// 若 Motifs 的跨度 >= MIN_SEG2 或 <= MAX_SEG2 (区间3)
			m_Result.m_Comments += "good-segment-length.";		// 则标记 good-segment-length. 

/*	Lim deleted
		if (m_Result.m_XXX == "CAB")							// 若比对为 RdRP , 但 Motif 顺序为 CAB ,
			{													// 则标记 permuted-penalty. , 最终得分减去 PERMUTED_PENALTY
			m_Result.m_Comments += "permuted-penalty.";
			m_Result.m_FinalScore -= PERMUTED_PENALTY;
			}
*/		
		if (m_Result.m_XXX == "ABC")	//	Lim added
			{
			int AB = m_Result.m_Dist12;
			int BC = m_Result.m_Dist23;
			if (AB < MIN_AB1 || AB > MAX_AB1)					// 若 Motif 1 和 2 之间的距离 < MIN_AB1 或 > MAX_AB1
				{												// 则标记 AB-dist-penalty. , 最终得分减去 AB_PENALTY1
				m_Result.m_Comments += "AB-dist-penalty.";
				m_Result.m_FinalScore -= AB_PENALTY1;
				}
			if (BC < MIN_BC1 || BC > MAX_BC1)					// 若 Motif 2 和 3 之间的距离 < MIN_BC1 或 > MAX_BC1
				{												// 则标记 BC-dist-penalty. , 最终得分减去 BC_PENALTY1
				m_Result.m_Comments += "BC-dist-penalty.";
				m_Result.m_FinalScore -= BC_PENALTY1;
				}
			}
		} // end if RdRP

	// 最终得分过低 (<0) , 0分
	if (m_Result.m_FinalScore < 0)								// 若 m_FinalScore 最终得分 < 0 (主要由于罚分 , 初始分值为 0)
		{														// 则 m_FinalScore <- 0 , 标记 m_Gene 为 unclassified , 标记 m_Comments 为 negative-score.
		m_Result.m_FinalScore = 0;								// 结束计算
		m_Result.m_Gene = "unclassified";
		m_Result.m_Comments += "negative-score.";
		return;
		}

	if (!AllowHigh && m_Result.m_FinalScore >= MIN_HICONF)		// 主要针对最低 Motif 得分过低 , 但总得分高的情况 
		m_Result.m_FinalScore = MIN_LOCONF;						// 将得分赋值为 MIN_LOCONF

	m_Result.m_HiConf = (m_Result.m_FinalScore >= myHICONF);	// 若 AllowHight 仍为 True 且 m_FinalScore >= Highly Confirm Score ,
	}															// m_HiConf <- True , 否则为 False

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

// 未命中结果的初始化
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

// 未命中结果的初始化
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
