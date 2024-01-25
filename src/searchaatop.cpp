#include "myutils.h"
#include "sfasta.h"
#include "pssm.h"
#include "pssmsearch.h"

#define TRACE	0


///
/// �����ѵ÷ֵİ���������Ƭ�Σ����ƣ�����Ϣ������ PSSMHitAA �ṹ��
/// 
/// P������������ PSSM ģ��
/// MinScore��\
/// Label�����������еı�ǩ
/// AASeq������������
/// L�����������еĳ���
/// Hit��PSSMHitAA �ṹ��

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
	if (M > L)		// Query ���еĳ��ȱ������ģ�ͳ���
		return;
	float BestScore = -9e9f;		// ��ʼ����ѵ÷�
	unsigned BestPos = UINT_MAX;	// ��ʼ�������ʼλ��
	for (unsigned Pos2 = 0; Pos2 <= L - M; ++Pos2)
		{
		float Score2 = P.GetScore(AASeq + Pos2);	// ��ַ���ͣ�ָ�벻�Ϻ��ƣ��� 0 �ӵ� L-M
		if (Score2 > BestScore)						// ��õ÷���ߵ�����Ƭ��
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


// ����� Query ���бȶԵ÷���ߵ� Group ���Լ���� group �� PSSMs �ıȶ���Ϣ
void SearchAATop()
	{
	asserta(optset_psm);
	const string InputFileName = string(opt_search_aa_top);		// �����ļ���
	const string OutputFileName = string(opt_output);			// ����ļ��� 
	asserta(InputFileName != "");
	asserta(OutputFileName != "");

	FILE *fOut = CreateStdioFile(OutputFileName);

	PSSM P;
	P.FromFile(opt_psm);						// ���ļ���ȡ PSSM ����
	unsigned M = P.GetColCount();				// ��� PSSM ���г���

	SFasta SF;					
	SF.Open(InputFileName);						// SFasta ��������ȡ fatsa �ļ�

	ProgressStep(0, 1002, "Search aa %s", SF.m_FileName.c_str());
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();		// SF �洢 fasta �ļ��ĸ�����Ϣ��Seq Ϊָ�����п�ͷ��ָ��
		if (Seq == 0)
			break;
		ProgressStep(SF.GetPctDoneX10(), 1002, "Search aa %s", SF.m_FileName.c_str());

		const string Label = SF.GetLabel();
		const unsigned L = SF.GetSeqLength();

		PSSMHitAA Hit;
		GetTopHitAA(P, (float) opt_minscore, Label, (const char *) Seq, L, Hit);	// Hit �������˸� Query ��������ߵ÷� Group �� PSSMs �ıȶ���Ϣ
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
