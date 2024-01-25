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

static FILE *g_fFev;		// .fev �ļ����
static FILE *g_fBed;		// �ļ����
static FILE *g_fTrimFa;		// �ļ������޼�������
static FILE *g_fTrimFaNt;	
static FILE *g_fMotFa;
static FILE *g_fMotFa2;
static FILE *g_fNotFa;		// δȷ��������
FILE *g_fPSSMAln;			// PSSM �ȶԽ���ļ�
FILE *g_fFullAln;			// ȫ���ȶԽ���ļ�
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

/// ����ȶ����ݵ� .txt �ļ�
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

// ͬ��
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
		Mod.SearchAA(Label, Seq);			// ���г��� < 500 �� Mod.SearchAA
//	else
//		Mod.SearchAA_Seeded(Label, Seq);	// ���г��� >= 500 �� Mod.SearchAA_Seeded
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

/// �� g_fTrimFa �ļ�����޼���� Query ���� (�������80���ַ�)
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
	omp_set_lock(&OutputLock);		// ��û�����	// Output() �е�������һ���߳�ִ��
	++g_QueryCount;					// Query ���м��� +1
	g_TotalLetters += SIZE(Seq);	// �ܰ�������� +size(Seq)

	const RPResult &Res = Mod.m_Result;		// ��ǰ���ж�ȡ���

// Confidence
	bool ShowConf = false;								// �Ƿ�չʾ���Ŷ�
	if (opt_all)										// opt_all == 1 ʱ , ������չʾ���Ŷ�
		ShowConf = true;
	if (Res.m_HiConf)									// m_HiConf == True ʱ , ������չʾ���Ŷ� 
		ShowConf = true;
	if (opt_loconf && Res.m_Gene != "unclassified")		// opt_loconf == 1 ���������յ÷� > 0 ʱ , ������չʾ���Ŷ� 
		ShowConf = true;

// Gene
	bool ShowGene = false;								// �Ƿ�չʾ����
	if (Res.m_Gene == "RT" && optset_rt)				// �����б��ж�Ϊ RT , �� optset_rt == 1 , չʾ������
		ShowGene = true;
	if (Res.m_Gene == "RdRP" && optset_rdrp)			// �����б��ж�Ϊ RdRP , �� optset_rdrp == 1 , չʾ������
		ShowGene = true;
	if (Res.m_Gene == "nifH" && optset_nifH)			// �����б��ж�Ϊ nifH , �� optset_nifH == 1 , չʾ������
		ShowGene = true;
//	if (!optset_rt && !optset_rdrp)						// optset_rdrp == 0 �� optset_rt == 0 , ����ʾ����ָ��һ������չʾ����
//		Die("Must specify one or both of -rt and -rdrp");
	if (!optset_nifH && !optset_rdrp && !optset_rt)
		Die("Must specify at least one gene.");

	const bool Show = (ShowConf && ShowGene);			// ���ŶȺ��Ƿ�չʾ���о�ͨ�� , �����չʾ

	if (Show)
		{
		uint IntScore = uint(Mod.m_Result.m_FinalScore);	// 1. ������յ÷�
		if (IntScore >= 100)								// 2. �����յ÷ִ��� 99 ����ֵΪ 99
			IntScore = 99;
		++(g_ScoreToCount[IntScore]);						// 3. Score �ֲ�ͳ�� , �� Score ��λ�� +1
		++g_FoundCount;										// 4. ���������ܼ��� +1
		
											// ���ļ�����ǿգ�
		OutputFev(Mod.m_Result);						// �� g_fFev �ļ�����ȶ��ܽ��ļ�, example��score=35.9	query=Y2R2_DROME/498-741	gene=RT	confidence=high	comments=high-PSSM-score.	order=ABC	qlen=244	pp_start=78	pp_end=172	pp_length=95	pssm_total_score=35.9	pssm_min_score=10.1	order=ABC	v1_length=46	v2_length=14	group=RT	motifs=TFVDFKGAFDNVxxxPQGSISGPFIWDILMxxxAYADDLLL
		OutputAln(Mod.m_Result);						// �� g_fRep �ļ�����ȶ����ݣ���� https://github.com/rcedgar/palmscan/blob/main/test/results/PF02123_RdRP_4.txt
		OutputFullAln(Mod.m_Result);					// �� g_fFullAln �ļ�����ȶ����ݣ�ͬ��
		OutputTrimFasta(Mod, Label);					// �� g_fTrimFa �ļ�����޼���� Query ���� (�������80���ַ�)
		OutputTrimFastaNt(Mod, Label);					// �� g_fTrimFaNt �ļ�����޼���� Nt ���� (�� m_TrimmedSeqNt �ǿ�) 
		OutputMotifsFasta(Mod, Label);					// �� g_fMotFa / g_fMotFa2 �ļ���� m_MotifsSeq >Label\n( A(Q) + "xxx" + B(Q) + "xxx" + C(Q) )
		OutputBed(Mod, Label);							// �� g_fBed �ļ���� Bed ��ʽ, example��Label    m_StartPos    m_EndPos    pp-RdRP    finalscore/60*1000    .\n
		}
	else
		SeqToFasta(g_fNotFa, Label.c_str(), Seq.c_str(), SIZE(Seq));

	omp_unset_lock(&OutputLock);	// �ͷŻ�����
	}

static void Search1(RdRpModel &Mod, const string &Label, string &Seq)	// ��װ�������ж��ǰ�������һ��Ǻ��������
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
	const string &QueryFileName = opt_search_pp;	// opt_search_pp Ϊ��ѯ�� Query ���м��ļ����
	const string &ModelFileName = opt_model;		// opt_model Ϊ pssm model �ļ����

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
	omp_init_lock(&OutputLock);							// ��ʼ��������

	const uint ThreadCount = GetRequestedThreadCount();
	vector<ObjMgr *> OMs;
	vector<RdRpModel *> Mods;
	for (int i = 0; i < int(ThreadCount); ++i)			// �����߳������������� 
		{
		ObjMgr *OM = new ObjMgr;						// �������������
		OMs.push_back(OM);
		RdRpModel *Mod = new RdRpModel;					// ���� RdRP ģ��
		Mod->m_Thread = i;								// �����̱߳��
		if (optset_model)
			Mod->FromModelFile(ModelFileName);			// ��ȡģ���ļ������� Mod ����RdRpModel �ࣩ
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
	else if (QueryFileName == ".calibrate.amino.")	// У׼������
		{
		g_Calibrate = true;
		RandomSeqSource *RSS = new RandomSeqSource;
		RSS->m_Nucleo = false;
		RSS->m_SeqCount = UINT_MAX-1;
		SS = RSS;
		MaxSecs = opt_secs;
		}
	else					// ��������ʱ�Ĳ���
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

#pragma omp parallel num_threads(ThreadCount)		// ���б�̣��� ThreadCount ���̲߳��� //
	{
	int ThreadIndex = omp_get_thread_num();			// ��ȡ��ǰ�߳�����
	RdRpModel &Mod = *Mods[ThreadIndex];			// Ϊ��ǰ�̷߳��� Mod ����
	ObjMgr *OM = OMs[ThreadIndex];					// Ϊ��ǰ�̷߳��� OM ���󣨶�����������������Ӽ�����������٣�
	for (;;)										// ѭ��ִ�����м���
		{
		if (Stop)
			break;

		SeqInfo *QSI = OM->GetSeqInfo();			// ����һ�� QSI��SeqInfo�����󣬲��ڵ��ȵ� busy ������

		if (g_QueryCount%100 == 0)
			CurrElapsedSecs = GetElapsedSecs();		// ÿ��ѯ 100 �����У�����һ��������ʱ�� GetElapsedSecs()

		if (g_Calibrate)	// ����У׼ģʽ�£�������ʱ�䳬�� MaxSecs ��ֹͣ
			{
			if (CurrElapsedSecs > MaxSecs)
				{
				Stop = true;
				break;
				}
			}

		if (ThreadIndex == 0 && CurrElapsedSecs > LastElapsedSecs)			// ���Ϊ 0 ���߳��� , ����ǰ����ʱ�䳬����һ�β�ѯ��������ʱ�� , ���ӡ�ļ���ȡ���ͷ�������ռ��
			{
			if (g_Calibrate)
				Progress("Searching %u/%u hits\r", g_FoundCount, g_QueryCount);
			else
				{
				uint Pct10 = SS->GetPctDoneX10();							// ��õ�ǰ�Ѷ�ȡ���ļ��� (��λΪǧ��֮һ)
				double HitPct = GetPct(g_FoundCount, g_QueryCount);			// ��õ�ǰ���ҵ�������ռ�����еı��� (��λΪ %)
				ProgressStep(Pct10, 1000, "Searching %u/%u hits (%.1f%%)",
					g_FoundCount, g_QueryCount, HitPct);
				}
			LastElapsedSecs = CurrElapsedSecs;								// ��¼��ǰ��ѯ��������ʱ��
			}

		bool Ok = SS->GetNext(QSI);			// �� SS��SeqSource�������һ��������Ϣ���� QSI��SeqInfo��
		if (!Ok)							// ����ȡ������һ����ֹͣ
			break;

		const string Label = string(QSI->m_Label);	// ��ȡ���еı�ǩ��Ϣ
		string Seq;
		for (uint i = 0; i < QSI->m_L; ++i)			// ��ȡ����
			Seq += char(QSI->m_Seq[i]);

		Search1(Mod, Label, Seq);			// ��Ҫ���� ���� Mod (RdRpModel ��)���Ա�ǩΪ Label �Ĳ�ѯ���� Seq ���м���
		Output(Mod, Label, Seq);			// ��Ҫ���� ����������
		OM->Down(QSI);
		}
	}
	double HitPct = GetPct(g_FoundCount, g_QueryCount);		// ��ðٷֱ� g_FoundCount / g_QueryCount *100������������ռ�ı�����
	ProgressStep(999, 1000, "Searching %u/%u hits (%.1f%%)",
	  g_FoundCount, g_QueryCount, HitPct);
	
	if (g_Calibrate)				// ��ӡ Score �� 0 �� >99 �ķֲ���Ϣ
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
