#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "rdrpmodel.h"
#include "pssmsearch.h"
#include "seqinfo.h"
#include "fastaseqsource.h"
#include "randomseqsource.h"
#include <time.h>

void Translate(const string &NtSeq, int Frame, string &AASeq);
bool GetIsNucleo(const string &Seq);
void SeqToUpper(string &Seq);

static FILE *g_fFev;		// .fev 文件句柄
static FILE *g_fBed;		// 文件句柄
static FILE *g_fTrimFa;		// 文件储存修剪的序列
static FILE *g_fTrimFaNt;	
static FILE *g_fMotFa;
static FILE *g_fMotFa2;
static FILE *g_fNotFa;		// 未确定的序列
FILE *g_fPSSMAln;			// PSSM 比对结果文件
FILE *g_fFullAln;			// 全部比对结果文件
static omp_lock_t ProgressLock;
static omp_lock_t OutputLock;
static uint g_FoundCount;
static uint g_QueryCount;
static bool g_Calibrate = false;
static uint64 g_TotalLetters = 0;
static uint32 *g_ScoreToCount[100];

static void OutputBed(const RdRpModel &Mod, const string &Label)
	{
	if (g_fBed == 0)
		return;

	const RPResult &Res = Mod.m_Result;

/***
1. Chrom.
2. chromStart 0-based first position.
3. chromEnd 1-based last position.
4. Name (high-confidence-RdRP)
5. Score (0 .. 1000)
6. Strand.	
***/

	asserta(Res.m_FinalScore >= 0);
	uint Score = uint((Res.m_FinalScore/60)*1000);
	if (Score > 1000)
		Score = 1000;

	string Feature = "pp-" + Res.m_Gene;
	fprintf(g_fBed, "%s", Label.c_str()); // 1
	if (Res.m_StartPosNt != UINT_MAX)
		{
		fprintf(g_fBed, "\t%u", Res.m_StartPosNt); // 2
		fprintf(g_fBed, "\t%u", Res.m_StartPosNt + Res.m_SegLengthNt); // 3
		fprintf(g_fBed, "\t%s", Feature.c_str()); // 4
		fprintf(g_fBed, "\t%u", Score); // 5
		fprintf(g_fBed, "\t%c", pom(Res.m_Frame > 0));
		}
	else
		{
		fprintf(g_fBed, "\t%u", Res.m_StartPos); // 2
		fprintf(g_fBed, "\t%u", Res.m_StartPos + Res.m_SegLength); // 3
		fprintf(g_fBed, "\t%s", Feature.c_str()); // 4
		fprintf(g_fBed, "\t%u", Score); // 5
		fprintf(g_fBed, "\t.");
		}
	fprintf(g_fBed, "\n");
	}

static void OutputFev(const RPResult &Res)
	{
	if (g_fFev == 0)
		return;

	string FEV;
	Res.ToFEVStr(FEV);
	fprintf(g_fFev, "%s\n", FEV.c_str());
	}

/// 输出比对数据到 .txt 文件
static void OutputAln(const RPResult &Res)
	{
	if (g_fRep == 0)
		return;

	const vector<string> &Aln = Res.m_Aln;
	if (Aln.empty())
		return;
	asserta(SIZE(Aln) == 4);
	string Cat;
	Res.GetCat(Cat);
	/// unclassified or high/low-confidence-RdRP/RT 

	fprintf(g_fRep, "\n");
	fprintf(g_fRep, ">%s\n", Res.m_QueryLabel.c_str());
	for (int i = 0; i < 4; ++i)
	  fprintf(g_fRep, "%s\n", Aln[i].c_str());
	fprintf(g_fRep, "Score %.1f, %s: %s\n",
	  Res.m_FinalScore,
	  Cat.c_str(),
	  Res.m_Comments.c_str());

	/// example:
	///
	///	>RDRP_ROTHD / 100 - 748
	///		A:238 - 249(14.8)      B : 312 - 325(19.0)      C : 349 - 356(11.5)
	///		LYTDVSQWDSSQ    <62> SGEKQTKAANSIAN  <23> VDGDDNYA[119]
	///		| ++ | . | +||||          |||+. | . + || +. | +.|||. | +
	///		lmlDysdFNsqH         SGeRaTtfiNsvlN       hvGDDvll
	///		Score 55.3, high - confidence - RdRP : high - PSSM - score.reward - DDGGDD.good - segment - length.
	}

// 同上
static void OutputFullAln(const RPResult &Res)
	{
	if (g_fFullAln == 0)
		return;

	const vector<string> &Aln = Res.m_FullAln;
	if (Aln.empty())
		return;
	asserta(SIZE(Aln) == 4);

	string Cat;
	Res.GetCat(Cat);

	fprintf(g_fFullAln, "\n");
	fprintf(g_fFullAln, ">%s\n", Res.m_QueryLabel.c_str());
	for (int i = 0; i < 4; ++i)
	  fprintf(g_fFullAln, "%s\n", Aln[i].c_str());
	fprintf(g_fFullAln, "Score %.1f, %s: %s\n",
	  Res.m_FinalScore,
	  Cat.c_str(),
	  Res.m_Comments.c_str());
	}

static void Search1_AA(RdRpModel &Mod, const string &Label, string &Seq)
	{
	uint QL = SIZE(Seq);
//	if (QL < 500)
		Mod.SearchAA(Label, Seq);			// 序列长度 < 500 用 Mod.SearchAA
//	else
//		Mod.SearchAA_Seeded(Label, Seq);	// 序列长度 >= 500 用 Mod.SearchAA_Seeded
	}

static void CorrectResultTranslated(const string &NtSeq, const string &AASeq,
  RPResult &Res, int Frame)
	{
	if (Res.m_FinalScore <= 0 || Res.m_Gene == "unclassified")
		return;

	uint QLnt = SIZE(NtSeq);
	uint QLaa = SIZE(AASeq);
	int SegLengthNt = 3*Res.m_SegLength;

	uint StartPosAA = Res.m_StartPos;
	uint EndPosAA = StartPosAA + Res.m_SegLength;
	asserta(EndPosAA <= SIZE(AASeq));

	asserta(Frame != 0);
	int Off = abs(Frame) - 1;
	if (Off < 0 || Off > 2)
		Die("Frame=%d, Off=%d", Frame, Off);

	int StartPosNt;
	int EndPosNt;
	if (Frame > 0)
		{
		StartPosNt = 3*int(Res.m_StartPos) + Off;
		EndPosNt = StartPosNt + SegLengthNt - 1;
		}
	else
		{
		EndPosNt = QLnt - int(3*Res.m_StartPos) - Off - 1;
		StartPosNt = EndPosNt - SegLengthNt + 1;
		}

	int SegLengthNt2 = EndPosNt - StartPosNt + 1;
	asserta(SegLengthNt2 == SegLengthNt);
	asserta(SegLengthNt > 1);
	asserta(StartPosNt >= 0);
	if (EndPosNt >= int(QLnt))
		Die("EndPosNt = %d, QLnt = %d", EndPosNt, QLnt);

	Res.m_Frame = Frame;
	Res.m_StartPosNt = uint(StartPosNt);
	Res.m_SegLengthNt = uint(SegLengthNt);

	string &t = Res.m_TrimmedSeqNt;
	t.clear();
	for (int Pos = StartPosNt; Pos <= EndPosNt; ++Pos)
		t += NtSeq[Pos];
	if (Frame < 0)
		RevCompSeq(t);
	}

static void Search1_Nt(RdRpModel &Mod, const string &Label, string &NtSeq)
	{
	const uint QLnt = SIZE(NtSeq);
	string AASeq;
	RPResult BestResult;
	BestResult.Clear();
	RdRpModel::SetResult_NoNtHit(Label, QLnt, BestResult);
	for (int Frame = 3; Frame >= -3; --Frame)
		{
		if (Frame == 0)
			continue;
		Translate(NtSeq, Frame, AASeq);
		const uint AAL = SIZE(AASeq);
		Search1_AA(Mod, Label, AASeq);
		if (Mod.m_Result.m_FinalScore >= BestResult.m_FinalScore)
			{
			BestResult = Mod.m_Result;
			CorrectResultTranslated(NtSeq, AASeq, BestResult, Frame);
			}
		}
	if (BestResult.m_FinalScore == -999)
		RdRpModel::SetResult_NoNtHit(Label, QLnt, Mod.m_Result);
	else
		{
		Mod.m_Result = BestResult;
		Mod.m_Result.m_QL = 0;
		Mod.m_Result.m_QLnt = QLnt;
		}
	}

static void MakeNewLabel(const RdRpModel &Mod, const string &Label,
  bool IsNucleo, string &NewLabel)
	{
	if (!opt_coords)
		{
		NewLabel = Label;
		return;
		}

	const RPResult &Res = Mod.m_Result;
	if (IsNucleo)
		{
		uint Start = Res.m_StartPosNt + 1;
		uint End = Start + Res.m_SegLengthNt;
		asserta(Res.m_Frame != 0);
		Ps(NewLabel, "%s %u-%u(%+d)",
		  Label.c_str(), Start, End, Res.m_Frame);
		}
	else
		{
		uint Start = Res.m_StartPos + 1;
		uint End = Start + Res.m_SegLength;
		Ps(NewLabel, "%s %u-%u",
		  Label.c_str(), Start, End);
		}
	}

/// 向 g_fTrimFa 文件输出修剪后的 Query 序列 (单行最大80个字符)
/// >Old_label startPos-endPos(m_frame)
/// ABDIYQOSNIFBVIQBIXNQBDFIQBQS...
static void OutputTrimFasta(RdRpModel &Mod, const string &Label)
	{
	if (g_fTrimFa == 0)
		return;

	const string &TrimSeq = Mod.m_Result.m_TrimmedSeq;
	if (TrimSeq.empty())
		return;

	string NewLabel;
	MakeNewLabel(Mod, Label, false, NewLabel);
	SeqToFasta(g_fTrimFa, NewLabel.c_str(), TrimSeq.c_str(), SIZE(TrimSeq));
	}

static void OutputTrimFastaNt(RdRpModel &Mod, const string &Label)
	{
	if (g_fTrimFaNt == 0)
		return;

	const string &TrimSeqNt = Mod.m_Result.m_TrimmedSeqNt;
	if (TrimSeqNt.empty())
		return;

	string NewLabel;
	MakeNewLabel(Mod, Label, true, NewLabel);
	SeqToFasta(g_fTrimFaNt, NewLabel.c_str(), TrimSeqNt.c_str(), SIZE(TrimSeqNt));
	}

static void OutputMotifsFasta(RdRpModel &Mod, const string &Label)
	{
	if (g_fMotFa != 0)
		{
		const string &MotifsSeq = Mod.m_Result.m_MotifsSeq;
		SeqToFasta(g_fMotFa, Label.c_str(), MotifsSeq.c_str(), SIZE(MotifsSeq));
		}
	if (g_fMotFa2 != 0)
		{
		const string &MotifsSeq2 = Mod.m_Result.m_MotifsSeq2;
		SeqToFasta(g_fMotFa2, Label.c_str(), MotifsSeq2.c_str(), SIZE(MotifsSeq2));
		}
	}

// opt_loconf, opt_rt, opt_rdrp
static void Output(RdRpModel &Mod, const string &Label, string &Seq)
	{
	omp_set_lock(&OutputLock);		// 获得互斥器	// Output() 中的内容由一个线程执行
	++g_QueryCount;					// Query 序列计数 +1
	g_TotalLetters += SIZE(Seq);	// 总氨基酸计数 +size(Seq)

	const RPResult &Res = Mod.m_Result;		// 当前序列读取结果

// Confidence
	bool ShowConf = false;								// 是否展示置信度
	if (opt_all)										// opt_all == 1 时 , 该序列展示置信度
		ShowConf = true;
	if (Res.m_HiConf)									// m_HiConf == True 时 , 该序列展示置信度 
		ShowConf = true;
	if (opt_loconf && Res.m_Gene != "unclassified")		// opt_loconf == 1 且序列最终得分 > 0 时 , 该序列展示置信度 
		ShowConf = true;

// Gene
	bool ShowGene = false;								// 是否展示基因
	if (Res.m_Gene == "RT" && optset_rt)				// 该序列被判定为 RT , 且 optset_rt == 1 , 展示该序列
		ShowGene = true;
	if (Res.m_Gene == "RdRP" && optset_rdrp)			// 该序列被判定为 RdRP , 且 optset_rdrp == 1 , 展示该序列
		ShowGene = true;
	if (Res.m_Gene == "nifH" && optset_nifH)			// 该序列被判定为 nifH , 且 optset_nifH == 1 , 展示该序列
		ShowGene = true;
//	if (!optset_rt && !optset_rdrp)						// optset_rdrp == 0 且 optset_rt == 0 , 则提示必须指定一种序列展示类型
//		Die("Must specify one or both of -rt and -rdrp");
	if (!optset_nifH && !optset_rdrp && !optset_rt)
		Die("Must specify at least one gene.");

	const bool Show = (ShowConf && ShowGene);			// 置信度和是否展示序列均通过 , 则进行展示

	if (Show)
		{
		uint IntScore = uint(Mod.m_Result.m_FinalScore);	// 1. 获得最终得分
		if (IntScore >= 100)								// 2. 若最终得分大于 99 ，则赋值为 99
			IntScore = 99;
		++(g_ScoreToCount[IntScore]);						// 3. Score 分布统计 , 该 Score 的位置 +1
		++g_FoundCount;										// 4. 发现序列总计数 +1
		
											// 若文件句柄非空：
		OutputFev(Mod.m_Result);						// 向 g_fFev 文件输出比对总结文件, example：score=35.9	query=Y2R2_DROME/498-741	gene=RT	confidence=high	comments=high-PSSM-score.	order=ABC	qlen=244	pp_start=78	pp_end=172	pp_length=95	pssm_total_score=35.9	pssm_min_score=10.1	order=ABC	v1_length=46	v2_length=14	group=RT	motifs=TFVDFKGAFDNVxxxPQGSISGPFIWDILMxxxAYADDLLL
		OutputAln(Mod.m_Result);						// 向 g_fRep 文件输出比对数据，详见 https://github.com/rcedgar/palmscan/blob/main/test/results/PF02123_RdRP_4.txt
		OutputFullAln(Mod.m_Result);					// 向 g_fFullAln 文件输出比对数据，同上
		OutputTrimFasta(Mod, Label);					// 向 g_fTrimFa 文件输出修剪后的 Query 序列 (单行最大80个字符)
		OutputTrimFastaNt(Mod, Label);					// 向 g_fTrimFaNt 文件输出修剪后的 Nt 序列 (若 m_TrimmedSeqNt 非空) 
		OutputMotifsFasta(Mod, Label);					// 向 g_fMotFa / g_fMotFa2 文件输出 m_MotifsSeq >Label\n( A(Q) + "xxx" + B(Q) + "xxx" + C(Q) )
		OutputBed(Mod, Label);							// 向 g_fBed 文件输出 Bed 格式, example：Label    m_StartPos    m_EndPos    pp-RdRP    finalscore/60*1000    .\n
		}
	else
		SeqToFasta(g_fNotFa, Label.c_str(), Seq.c_str(), SIZE(Seq));

	omp_unset_lock(&OutputLock);	// 释放互斥器
	}

static void Search1(RdRpModel &Mod, const string &Label, string &Seq)	// 包装函数，判断是氨基酸查找还是核苷酸查找
	{
	SeqToUpper(Seq);
	bool IsNucleo = GetIsNucleo(Seq);
	if (IsNucleo)
		Search1_Nt(Mod, Label, Seq);
	else
		Search1_AA(Mod, Label, Seq);
	}

void SearchPP()
	{
	const string &QueryFileName = opt_search_pp;	// opt_search_pp 为查询的 Query 序列集文件句柄
	const string &ModelFileName = opt_model;		// opt_model 为 pssm model 文件句柄

	if (!opt_notrunclabels)
		opt_trunclabels = true;

	g_fFev = CreateStdioFile(opt_fevout);
	g_fRep = CreateStdioFile(opt_report);
	g_fTrimFa = CreateStdioFile(opt_ppout);
	g_fTrimFaNt = CreateStdioFile(opt_ppout_nt);
	g_fMotFa = CreateStdioFile(opt_motifs_fastaout);
	g_fMotFa2 = CreateStdioFile(opt_motifs_fastaout2);
	g_fBed = CreateStdioFile(opt_bedout);
	g_fNotFa = CreateStdioFile(opt_nohit_fastaout);
	g_fFullAln = CreateStdioFile(opt_alnout);
	g_fPSSMAln = CreateStdioFile(opt_pssm_alnout);

	omp_init_lock(&ProgressLock);
	omp_init_lock(&OutputLock);							// 初始化互斥器

	const uint ThreadCount = GetRequestedThreadCount();
	vector<ObjMgr *> OMs;
	vector<RdRpModel *> Mods;
	for (int i = 0; i < int(ThreadCount); ++i)			// 根据线程数量创建对象 
		{
		ObjMgr *OM = new ObjMgr;						// 创建对象管理器
		OMs.push_back(OM);
		RdRpModel *Mod = new RdRpModel;					// 创建 RdRP 模型
		Mod->m_Thread = i;								// 设置线程编号
		if (optset_model)
			Mod->FromModelFile(ModelFileName);			// 读取模型文件，构建 Mod 对象（RdRpModel 类）
		else
			{
			extern vector<string> g_ModelStrings;
			Mod->FromStrings(g_ModelStrings);
			}

		Mods.push_back(Mod);
		}

	uint MaxSecs = 0;
	SeqSource *SS;
	if (QueryFileName == ".calibrate.nucleo.")
		{
		g_Calibrate = true;
		RandomSeqSource *RSS = new RandomSeqSource;
		RSS->m_Nucleo = true;
		RSS->m_SeqCount = UINT_MAX-1;
		SS = RSS;
		MaxSecs = opt_secs;
		}
	else if (QueryFileName == ".calibrate.amino.")	// 校准氨基酸
		{
		g_Calibrate = true;
		RandomSeqSource *RSS = new RandomSeqSource;
		RSS->m_Nucleo = false;
		RSS->m_SeqCount = UINT_MAX-1;
		SS = RSS;
		MaxSecs = opt_secs;
		}
	else					// 正常运行时的步骤
		{	
		g_Calibrate = false;
		FASTASeqSource *FSS = new FASTASeqSource;
		FSS->Open(QueryFileName);
		SS = FSS;
		MaxSecs = 0;
		}
	bool Stop = false;
	uint LastElapsedSecs = 0;
	uint CurrElapsedSecs = 0;
	ProgressStep(0, 1000, "Searching");

#pragma omp parallel num_threads(ThreadCount)		// 并行编程，用 ThreadCount 个线程并行 //
	{
	int ThreadIndex = omp_get_thread_num();			// 获取当前线程索引
	RdRpModel &Mod = *Mods[ThreadIndex];			// 为当前线程分配 Mod 对象
	ObjMgr *OM = OMs[ThreadIndex];					// 为当前线程分配 OM 对象（对象管理器，引用数加减，对象的销毁）
	for (;;)										// 循环执行序列检索
		{
		if (Stop)
			break;

		SeqInfo *QSI = OM->GetSeqInfo();			// 创建一个 QSI（SeqInfo）对象，并在调度到 busy 链表中

		if (g_QueryCount%100 == 0)
			CurrElapsedSecs = GetElapsedSecs();		// 每查询 100 条序列，计算一次总运行时间 GetElapsedSecs()

		if (g_Calibrate)	// 若在校准模式下，若运行时间超过 MaxSecs 则停止
			{
			if (CurrElapsedSecs > MaxSecs)
				{
				Stop = true;
				break;
				}
			}

		if (ThreadIndex == 0 && CurrElapsedSecs > LastElapsedSecs)			// 编号为 0 的线程中 , 若当前运行时间超过上一次查询的总运行时间 , 则打印文件读取量和发现序列占比
			{
			if (g_Calibrate)
				Progress("Searching %u/%u hits\r", g_FoundCount, g_QueryCount);
			else
				{
				uint Pct10 = SS->GetPctDoneX10();							// 获得当前已读取的文件量 (单位为千分之一)
				double HitPct = GetPct(g_FoundCount, g_QueryCount);			// 获得当前已找到的序列占总序列的比例 (单位为 %)
				ProgressStep(Pct10, 1000, "Searching %u/%u hits (%.1f%%)",
					g_FoundCount, g_QueryCount, HitPct);
				}
			LastElapsedSecs = CurrElapsedSecs;								// 记录当前查询的总运行时间
			}

		bool Ok = SS->GetNext(QSI);			// 用 SS（SeqSource）获得下一条序列信息存入 QSI（SeqInfo）
		if (!Ok)							// 若获取不到下一条则停止
			break;

		const string Label = string(QSI->m_Label);	// 获取序列的标签信息
		string Seq;
		for (uint i = 0; i < QSI->m_L; ++i)			// 获取序列
			Seq += char(QSI->m_Seq[i]);

		Search1(Mod, Label, Seq);			// 主要函数 基于 Mod (RdRpModel 类)，对标签为 Label 的查询序列 Seq 进行检索
		Output(Mod, Label, Seq);			// 主要函数 制作输出结果
		OM->Down(QSI);
		}
	}
	double HitPct = GetPct(g_FoundCount, g_QueryCount);		// 获得百分比 g_FoundCount / g_QueryCount *100（发现序列所占的比例）
	ProgressStep(999, 1000, "Searching %u/%u hits (%.1f%%)",
	  g_FoundCount, g_QueryCount, HitPct);
	
	if (g_Calibrate)				// 打印 Score 从 0 到 >99 的分布信息
		{
		ProgressLog("\n");
		ProgressLog("Total letters %s\n", Int64ToStr(g_TotalLetters));
		Log("Score	Count\n");
		for (uint IntScore = 0; IntScore < 100; ++IntScore)
			Log("%u	%u\n", IntScore, g_ScoreToCount[IntScore]);
		}

	CloseStdioFile(g_fFev);
	CloseStdioFile(g_fRep);
	CloseStdioFile(g_fTrimFa);
	CloseStdioFile(g_fTrimFaNt);
	CloseStdioFile(g_fMotFa);
	CloseStdioFile(g_fMotFa2);
	CloseStdioFile(g_fBed);
	CloseStdioFile(g_fNotFa);
	CloseStdioFile(g_fFullAln);
	CloseStdioFile(g_fPSSMAln);
	}
