#include "myutils.h"
#include "sfasta.h"
#include "pssm.h"
#include "pssmsearch.h"

#define TRACE	0


///
/// 获得最佳得分的氨基酸序列片段（掌纹），信息储存于 PSSMHitAA 结构体
/// 
/// P：用于评估的 PSSM 模型
/// MinScore：\
/// Label：氨基酸序列的标签
/// AASeq：氨基酸序列
/// L：氨基酸序列的长度
/// Hit：PSSMHitAA 结构体

void GetTopHitAA(const PSSM &P, float MinScore, 
  const string &Label, const char *AASeq, unsigned L,
  PSSMHitAA &Hit)
	{
	Hit.Label = &Label;
	Hit.AASeq = AASeq;
	Hit.L = L;
	Hit.Score = -9e9f;
	Hit.AAPos = UINT_MAX;
	asserta(!P.m_IsNucleo);
	const unsigned M = P.GetColCount();
	if (M > L)		// Query 序列的长度必须大于模型长度
		return;
	float BestScore = -9e9f;		// 初始化最佳得分
	unsigned BestPos = UINT_MAX;	// 初始化最佳起始位置
	for (unsigned Pos2 = 0; Pos2 <= L - M; ++Pos2)
		{
		float Score2 = P.GetScore(AASeq + Pos2);	// 地址类型，指针不断后移，从 0 加到 L-M
		if (Score2 > BestScore)						// 获得得分最高的序列片段
			{
			BestScore = Score2;
			BestPos = Pos2;
			}
		}
	Hit.Score = BestScore;
	Hit.AAPos = BestPos;
#if	TRACE
	Log("GetTopHitAA Score %.2f, MinScore %.2f, Pos %u\n",
	  BestScore, MinScore, BestPos);
#endif
	}


// 获得与 Query 序列比对得分最高的 Group ，以及与该 group 中 PSSMs 的比对信息
void SearchAATop()
	{
	asserta(optset_psm);
	const string InputFileName = string(opt_search_aa_top);		// 输入文件名
	const string OutputFileName = string(opt_output);			// 输出文件名 
	asserta(InputFileName != "");
	asserta(OutputFileName != "");

	FILE *fOut = CreateStdioFile(OutputFileName);

	PSSM P;
	P.FromFile(opt_psm);						// 从文件读取 PSSM 矩阵
	unsigned M = P.GetColCount();				// 获得 PSSM 序列长度

	SFasta SF;					
	SF.Open(InputFileName);						// SFasta 类用来读取 fatsa 文件

	ProgressStep(0, 1002, "Search aa %s", SF.m_FileName.c_str());
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();		// SF 存储 fasta 文件的各类信息，Seq 为指向序列开头的指针
		if (Seq == 0)
			break;
		ProgressStep(SF.GetPctDoneX10(), 1002, "Search aa %s", SF.m_FileName.c_str());

		const string Label = SF.GetLabel();
		const unsigned L = SF.GetSeqLength();

		PSSMHitAA Hit;
		GetTopHitAA(P, (float) opt_minscore, Label, (const char *) Seq, L, Hit);	// Hit 对象获得了该 Query 序列与最高得分 Group 的 PSSMs 的比对信息
		unsigned Pos = Hit.AAPos;
		if (Pos == UINT_MAX)
			continue;

		fprintf(fOut, "%s", Label.c_str());
		fprintf(fOut, "\t%u", Pos + 1);
		fprintf(fOut, "\t%.3f", Hit.Score);
		fprintf(fOut, "\t%*.*s", M, M, Seq + Pos);
		fprintf(fOut, "\n");
		}
	ProgressStep(1001, 1002, "Search aa %s", SF.m_FileName.c_str());

	CloseStdioFile(fOut);
	}
