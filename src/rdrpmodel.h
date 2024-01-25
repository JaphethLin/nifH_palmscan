#pragma once

#include "pssm.h"
#include "pssmsearch.h"
#include "alnparams.h"
#include "xdpmem.h"
#include "pathinfo.h"
#include "objmgr.h"
#include "rphit.h"
#include "rpresult.h"
#include "heuristics.h"

class RdRpModel
	{
public:
	int m_Thread;						// 线程编号
	string m_QueryLabel;
	string m_QuerySeq;

	vector<string> m_GroupNames;		// 'Duplorna', 'Birna', 'Lena', ...
	vector<char> m_MotifLetters;		// 'A', 'B', 'C'

	vector<PSSM> m_PSSMs;				// [1]PSSM, [2]PSSM, [3]PSSM, ...
	vector<string> m_PSSMGroupNames;	// 'Duplorna', 'Duplorna', 'Duplorna', 'Lena', 'Lena', 'Lena', ...
	vector<char> m_PSSMMotifLetters;	// 'A', 'A', 'A', 'B', 'B', 'B', ...
	vector<RPHit> m_HitsU;				// UnGapped Hits
	vector<RPHit> m_HitsG;				// Gapped Hits

	vector<vector<uint> > m_PIV;		// 存储了 motif 和 group 之间的映射关系：m_VIP[motif index][group index] => PSSMIndex

	bool m_Gapped;
	ObjMgr m_OM;
	AlnParams m_AP;		// 比对的参数设置
	XDPMem m_DPMem;

	uint m_MaxX;

	RPResult m_Result;

public:
	RdRpModel()
		{
		m_Gapped = opt_gapped;
		m_MaxX = 10;
		if (optset_maxx)
			m_MaxX = opt_maxx;
		float GapOpen = -10;	// 打开第一个 gap 时候的罚分
		float GapExt = -1;		// 延伸 gap 时候的罚分
		float TermGapOpen = 0;	// 末端打开第一个 gap 时候的罚分
		float TermGapExt = 0;	// 末端延伸 gap 时候的罚分
		m_AP.Init4(0, GapOpen, GapExt, TermGapOpen, TermGapExt);
		m_AP.LOpenA = -999;
		m_AP.ROpenA = -999;
		}

	void Clear()
		{
		m_PSSMs.clear();
		m_PSSMGroupNames.clear();
		m_PSSMMotifLetters.clear();
		}

	const PSSM &GetGetPSSM(uint Index) const	// 获得索引为 Index 的 PSSM 对象 
		{
		asserta(Index < SIZE(m_PSSMs));
		return m_PSSMs[Index];
		}

	uint GetMotifCount() const					// 获得 Motif 的数量
		{
		return SIZE(m_MotifLetters);
		}

	uint GetGroupCount() const					// 获得 Group 的数量
		{
		return SIZE(m_GroupNames);
		}

	uint GetPSSMCount() const					// 获得 PSSM 对象的数量
		{
		return SIZE(m_PSSMs);
		}

	const PSSM &GetPSSM(uint PSSMIndex) const	// 获得索引为 PSSMIndex 的 PSSM 对象 
		{
		asserta(PSSMIndex < SIZE(m_PSSMs));
		return m_PSSMs[PSSMIndex];
		}

	uint GetPSSMIndex(uint MotifIndex, uint GroupIndex) const	// 获得目标 PSSM 的索引
		{
		asserta(MotifIndex < SIZE(m_PIV));
		asserta(GroupIndex < SIZE(m_PIV[0]));
		uint PSSMIndex = m_PIV[MotifIndex][GroupIndex];
		return PSSMIndex;
		}

	const PSSM *GetPSSM(uint MotifIndex, uint GroupIndex) const	// 获得目标 PSSM 对象
		{
		uint PSSMIndex = GetPSSMIndex(MotifIndex, GroupIndex);
		if (PSSMIndex == UINT_MAX)
			return 0;
		asserta(PSSMIndex < SIZE(m_PSSMs));
		return &m_PSSMs[PSSMIndex];
		}

	void GetGroupName(uint GroupIndex, string &GroupName) const	// 获得 motif 分组的名称
		{
		if (GroupIndex >= SIZE(m_GroupNames))
			{
			GroupName = "*";
			return;
			}
		GroupName = m_GroupNames[GroupIndex];
		}

	uint GetPSSMLength(char MotifLetter, uint GroupIndex) const
		{
		return GetPSSMLength(GetMotifIndex(MotifLetter), GroupIndex);
		}

	// 根据 Motif 和 group 索引获得对应的命中对象
	const RPHit *GetHit(char MotifLetter, uint GroupIndex) const
		{
		return GetHit(GetMotifIndex(MotifLetter), GroupIndex);
		}

	uint GetTopGroup() const
		{
		uint TopIx;
		float Score;
		GetTopGroup(TopIx, Score);
		return TopIx;
		}

	const RPHit *GetHit(char MotifLetter) const
		{
		uint GroupIndex = GetTopGroup();
		if (GroupIndex == UINT_MAX)
			return 0;
		return GetHit(MotifLetter, GroupIndex);
		}

	uint GetPSSMLength(char MotifLetter) const
		{
		uint GroupIndex = GetTopGroup();
		if (GroupIndex == UINT_MAX)
			return 0;
		return GetPSSMLength(MotifLetter, GroupIndex);
		}

public:
	void FromSpecFile(const string &FileName);
	void FromModelFile(const string &FileName);
	void FromStrings(const vector<string> &Lines);
	void ToModelFile(const string &FileName) const;

	void SearchAA(const string &QueryLabel, const string &QuerySeq);
	void SearchAA_Seeded(const string &QueryLabel, const string &QuerySeq);

public:
	void SetPIV();
	uint GetMotifIndex(char MotifLetter) const;
	uint GetGroupIndex(const string &GroupName) const;
	uint GetPSSMIndex(char MotifLetter, const string &GroupName) const;
	uint GetPSSMLength(uint MotifIndex, uint GroupIndex) const;
	void SearchAA1(uint PSSMIndex, const string &Seq, RPHit &Hit);
	void SearchAA1_Ungapped(uint PSSMIndex, const string &Seq, RPHit &Hit);
	void SearchAA1_Gapped(uint PSSMIndex, const string &Seq, RPHit &Hit);
	void LogAlnViterbi(uint PSSMIndex, const string &Seq, const RPHit &Hit) const;
	void GetAln(uint MotifIndex, uint GroupIndex,
	  string &Q, string &P, string &Annot,
	  uint &QLo, uint &QHi) const;
	void GetXSeq(char MotifLetter, string &XSeq) const;
	int GetDist(const RPHit *Hit1, const RPHit *Hit2) const;
	int GetDist(char MotifLetter1, char MotifLetter2) const;
	uint GetSpan(const RPHit *Hit1, const RPHit *Hit2) const;

	void LogHitTable() const;
	float GetTotalScore(uint GroupIndex) const;
	float GetMinScore() const;
	uint GetTotalGaps() const;
	void GetTopGroup(uint &TopIx, float &TopScore) const;
	const RPHit *GetHit(uint MotifIndex, uint GroupIndex) const;
	void GetABCOrder(string &XXX) const;
	void GetOrderedHits(string &XXX, const RPHit *&Hit1,
	  const RPHit *&Hit2, const RPHit *&Hit3) const;
	bool GetABCRange(uint &Start, uint &End) const;
	void GetAlnRows(vector<string> &Rows) const;
	void GetFullAln(vector<string> &Rows) const;
	void GetMotifsSeq(string &Seq) const;
	void GetMotifsSeq2(string &Seq) const;
	void GetTrimmedSeq(string &Seq) const;
	float GetCScore() const;
	float GetBScore() const;	// Lim edited
	float GetAScore() const;	// Lim edited
	void GetSuperMotif(const string &MotifsSeq, string &s) const;

public:
	static uint GetNextSeedPos(const string &QuerySeq, uint SeedPos);
	static void SetResult_NoNtHit(const string &QueryLabel, uint QLnt,
	  RPResult &Result);
	static void SetResult_NoAaHit(const string &QueryLabel, uint QLnt,
	  RPResult &Result);

private:
	void SetResult();
	void CheckThread() const { asserta(omp_get_thread_num() == m_Thread); }
	};
