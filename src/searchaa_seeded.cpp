#include "myutils.h"
#include "rdrpmodel.h"

// 扫描序列，若包含子序列 GDD , SDD 或 GDN，则返回子序列中第一个氨基酸的位置，若未找到则返回 UINT_MAX
uint RdRpModel::GetNextSeedPos(const string &QuerySeq, uint SeedPos)
	{
	const uint QL = SIZE(QuerySeq);
	const char *Q = QuerySeq.c_str();
	while (SeedPos + 2 < QL)
		{
		char c1 = Q[SeedPos];
		char c2 = Q[SeedPos+1];
		char c3 = Q[SeedPos+2];
#define	eq(x, y, z)	if (c1 == x && c2 == y && c3 == z) return SeedPos;
		eq('C', 'G', 'G')
#undef eq
		++SeedPos;
		}
	return UINT_MAX;
	}

static void CorrectResultSeeded(RPResult &Res, uint LoPos)
	{
	Res.m_StartPos += LoPos;
	}

void RdRpModel::SearchAA_Seeded(const string &QueryLabel, const string &QuerySeq)
	{
	RPResult BestResult;
	BestResult.Clear();
	const uint QL = SIZE(QuerySeq);
	SetResult_NoAaHit(QueryLabel, QL, m_Result);			// 以未找到的形式初始化比对结果
	SetResult_NoAaHit(QueryLabel, QL, BestResult);
	uint SeedPos = 0;
	for (;;)
		{
		SeedPos = GetNextSeedPos(QuerySeq, SeedPos);		// 查找子序列 GDD , SDD 或 GDN 的起始位置，若未找到则返回 UINT_MAX
		if (SeedPos == UINT_MAX)
			break;

		int LoPos = int(SeedPos) - 350;						// 裁剪包含 GDD , SDD 或 GDN 的序列区间（长度 < 509）
		int HiPos = int(SeedPos) + 159;
		if (LoPos < 0)
			LoPos = 0;
		if (HiPos >= int(QL))
			HiPos = QL- 1;
		int ChunkL = HiPos - LoPos + 1;
		if (ChunkL < int(MIN_SEG))
			break;
		string Chunk;
		for (int i = LoPos; i <= HiPos; ++i)				// 获得剪裁后的序列
			Chunk += QuerySeq[i];

		SearchAA(QueryLabel, Chunk);						// 检索剪裁后的序列
		if (m_Result.m_FinalScore > BestResult.m_FinalScore)
			{
			BestResult = m_Result;							// 获得多次剪裁过程中最佳得分的片段序列
			asserta(LoPos >= 0);
			CorrectResultSeeded(BestResult, uint(LoPos));	// 矫正序列的起始位置（因为从 LoPos 位点开始剪切）
			}

		SeedPos += 50;										// 位置后移 50 pb，继续检索片段
		}
	m_Result = BestResult;
	}
