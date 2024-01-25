#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "alpha.h"
#include "sort.h"

void PSSM::SetIsNucleo(const vector<string> &Seqs)
	{
	unsigned CharCount = 0;
	unsigned NucleoCharCount = 0;
	unsigned SeqCount = SIZE(Seqs);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Seq = Seqs[SeqIndex];
		const unsigned L = SIZE(Seq);
		for (unsigned i = 0; i < L; ++i)
			{
			char c = Seq[i];
			++CharCount;
			if (g_IsNucleoChar[c])
				++NucleoCharCount;
			}
		}
	m_IsNucleo = (GetPct(NucleoCharCount, CharCount) > 80.0);
	}


/// 
///		pos1	pos2	pos3	pos4	pos5	pos6	pos7	pos8	pos9	����> x ά�ȱ���
///	A
///	T
///	C
/// G
/// 
///

///
///key algorithm 
///

// ���� SeqDB ���ȡ fasta �ļ������� PSSM ���ȡ SeqDB �е�������Ϣ
void PSSM::FromSeqs(const vector<string> &Labels, const vector<string> &Seqs)
	{
	//m_Labels = Labels;
	//m_Seqs = Seqs;
	unsigned SeqCount = SIZE(Seqs); // ������
	asserta(SeqCount > 0);			// �������� > 0
	m_ColCount = SIZE(Seqs[0]);		// ���г���
	for (unsigned SeqIndex = 1; SeqCount < SeqCount; ++SeqIndex)
		asserta(SIZE(Seqs[SeqIndex]) == m_ColCount);				// �����������г���һ��
	SetIsNucleo(Seqs);				// �ж��Ƿ����������
	SetCounts(Seqs);				// 
	SetFreqs();
	SetPseudoFreqs();
	AddPseudo(0.1f);
	SetScores();
	CalcConsSeq(m_ConsSeq);
	}

// ���� SeqDB ���ȡ fasta �ļ������� PSSM ���ȡ SeqDB �е�������Ϣ
void PSSM::FromFasta(const string &FileName)
	{
	SeqDB DB;
	DB.FromFasta(FileName);
	FromSeqs(DB.m_Labels, DB.m_Seqs);
	}

// PSSM ͳ��Ƶ��
void PSSM::SetCounts(const vector<string> &Seqs)
	{
	unsigned ColCount = GetColCount();		// ������е� m_ColCount ���������г��ȣ�
	m_Counts.clear();						// ��� m_count ��ά vector 
	m_Counts.resize(ColCount);				// ���� m_count ��ά vector ��x ά��Ϊ���г���
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)	// ѭ������ x ά�ȣ����г��ȣ�
		SetCountsCol(Seqs, ColIndex);		// Ϊ x ά�ȵ�ÿ�� vector ��ֵ
	}

// ����ÿ��λ���ϸ�������/������ĳ���Ƶ��
void PSSM::SetCountsCol(const vector<string> &Seqs, unsigned ColIndex)
	{
	unsigned AlphaSize = (m_IsNucleo ? 4 : 20);		// ��ĸ���С
	const byte *CharToLetter = (m_IsNucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);	// �ж���ʲô��ĸ��
	unsigned SeqCount = SIZE(Seqs);		// ��ȡ����������

	vector<unsigned> &Counts = m_Counts[ColIndex];
	Counts.clear();				// ��� y ά��
	Counts.resize(AlphaSize);	// y ά��Ϊ��ĸ���С
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex) // �������ж�ȡ������
		{
		asserta(ColIndex < SIZE(Seqs[SeqIndex]));	// ��ȡ�����г���һ����С�ڣ��������Ѷ�������г���
		char c = Seqs[SeqIndex][ColIndex];			// ��ȡ�� SeqIndex �����У��� ColIndex λ�õ��ַ�
		unsigned Letter = CharToLetter[c];			// ת��Ϊ�������ֵ���
		if (Letter >= AlphaSize)					// ��ֵ���һ��С����ĸ���С
			continue;
		++(Counts[Letter]);		// ��Ӧλ��ͳ��Ƶ�� +1
		}
	}

// PSSM ͳ��Ƶ��
void PSSM::SetFreqs()
	{
	unsigned ColCount = GetColCount();		// ���Ƶ�ξ���
	m_RawFreqs.clear();						// ����Ƶ�ʾ���
	m_RawFreqs.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)	// ѭ������ x ά�ȣ����г��ȣ�
		SetFreqsCol(ColIndex);										// ����ÿ��λ���ϵĸ�������Ƶ��		
	}

// ����ÿ��λ���ϸ��������Ƶ��
void PSSM::SetFreqsCol(unsigned ColIndex)
	{
	unsigned AlphaSize = (m_IsNucleo ? 4 : 20);

	vector<float> &Freqs = m_RawFreqs[ColIndex];
	Freqs.clear();
	Freqs.resize(AlphaSize);

	vector<unsigned> &Counts = m_Counts[ColIndex];
	unsigned N = 0;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)		// �Ե�һλ����ÿ�ְ������Ƶ����ͣ�N��
		{
		unsigned n = Counts[Letter];
		N += n;
		}
	if (N == 0)
		Die("PSSM::SetFreqsCol(%u): no valid letters", ColIndex);	// ����λ��û����ĸ���е��κΰ�������֣��򱨴�

	float Sum = 0.0;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		unsigned n = Counts[Letter];
		float f = N == 0 ? 0.0f : float(n)/N;		// ����Ƶ��
		Sum += f;
		Freqs[Letter] = f;
		}
	asserta(Sum > 0.99f && Sum < 1.01f);	// Ƶ��֮�͵ķ�Χ
	}

/// ���� m_PseudoFreqs ���� ////
void PSSM::SetPseudoFreqs()
	{
	unsigned AlphaSize = GetAlphaSize();	// �ж��ǰ����ỹ�Ǻ���������
	m_PseudoFreqs.clear();
	unsigned ColCount = GetColCount();		// ������г���
	m_PseudoFreqs.resize(ColCount);			// ���ö�ά����
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)	// ѭ������ x ά�ȣ����г��ȣ�
		{
		m_PseudoFreqs[ColIndex].resize(AlphaSize, 0.0f);	// y ά������Ϊ��ĸ���С
		SetPseudoFreqsCol(ColIndex);
		}
	}


/// ����һ��λ���ϸ�������� PseudoFreq { P( x -> y | x ) �� P( x ) }
void PSSM::SetPseudoFreqsCol(unsigned ColIndex)
	{
	unsigned AlphaSize = GetAlphaSize();				// �ж��ǰ����ỹ�Ǻ���������
	const vector<float> &Freqs = m_RawFreqs[ColIndex];	// ����Ƶ��ͳ�ƾ��� ColIndex λ�õ�����
	vector<float> &PseudoFreqs = m_PseudoFreqs[ColIndex];	// ������Ƶ��ͳ�ƾ��� ColIndex λ�õ�����
	asserta(SIZE(Freqs) == AlphaSize);					// ���ݳ���һ��Ϊ��ĸ����
	asserta(SIZE(PseudoFreqs) == AlphaSize);			// ���ݳ���һ��Ϊ��ĸ����
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)	
		{
		float f = Freqs[Letter];						// ��ð����ᣨLetter���ĳ���Ƶ��
		float SumJ = 0.0f;
		for (unsigned Letter2 = 0; Letter2 < AlphaSize; ++Letter2)
			{
			float JP = GetJP(Letter, Letter2);			// ��� Letter ���滻Ϊ Letter2 ����������
			SumJ += JP;
			float Incf = f*JP;
			PseudoFreqs[Letter2] += Incf;
			}
		asserta(SumJ > 0.99f && SumJ < 1.01f);
		}

	float Sum = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Sum += PseudoFreqs[Letter];
	asserta(Sum > 0.99f && Sum < 1.01f);
	}


// joint prob, sums to 1 for each Row, not whole Mx!!
float PSSM::GetJP(unsigned Letter1, unsigned Letter2) const
{
	extern double g_BLOSUM62_JP_Rows[20][20];
	extern double g_NtMx_JP_Rows[4][4];
	if (m_IsNucleo)
	{
		asserta(Letter1 < 4 && Letter2 < 4);
		return (float)g_NtMx_JP_Rows[Letter1][Letter2];
	}
	else
	{
		asserta(Letter1 < 20 && Letter2 < 20);
		return (float)g_BLOSUM62_JP_Rows[Letter1][Letter2];
	}
}


// ���� m_Freqs ���� ���� w ΪȨ�ؼ�Ȩ m_PseudoFreqs ����//
void PSSM::AddPseudo(float w)
{
	m_PseudoWeight = w;
	unsigned ColCount = GetColCount();
	unsigned AlphaSize = GetAlphaSize();
	m_Freqs.clear();
	m_Freqs.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
	{
		m_Freqs[ColIndex].resize(AlphaSize);
		AddPseudoCol(ColIndex, w);
	}
}

// �� w ΪȨ�ؼ�Ȩ m_PseudoFreqs ���󣨵�һλ�ô���
void PSSM::AddPseudoCol(unsigned ColIndex, float w)
{
	unsigned AlphaSize = GetAlphaSize();
	const vector<float>& RawFreqs = m_RawFreqs[ColIndex];
	const vector<float>& PseudoFreqs = m_PseudoFreqs[ColIndex];
	vector<float>& Freqs = m_Freqs[ColIndex];
	asserta(SIZE(RawFreqs) == AlphaSize);
	asserta(SIZE(PseudoFreqs) == AlphaSize);
	asserta(SIZE(Freqs) == AlphaSize);

	// Weighted sum
	float Sum = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
	{
		float RawFreq = RawFreqs[Letter];
		float PseudoFreq = PseudoFreqs[Letter];
		float Freq = RawFreq + w * PseudoFreq;		// �� w ΪȨ�ؼ�Ȩ����Ƶ�ʾ���
		Freqs[Letter] = Freq;
		Sum += Freq;
	}

	// Normalize
	float Sum2 = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
	{
		float Freq = Freqs[Letter] / Sum;
		Freqs[Letter] = Freq;
		Sum2 += Freq;
	}
	asserta(Sum2 > 0.99f && Sum2 < 1.01f);
}


// ��ȡԭʼ����Ƶ��
const vector<float> &PSSM::GetFreqs(unsigned ColIndex) const
	{
	asserta(ColIndex < SIZE(m_RawFreqs));
	return m_RawFreqs[ColIndex];
	}



void PSSM::LogFreqs(const vector<vector<float> > &FreqsVec) const
	{
	const unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();
	const unsigned ColCount = GetColCount();

	Log("Col");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" %5c", LetterToChar[Letter]);
	Log("\n");

	Log("---");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" -----");
	Log("\n");

	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		Log("%3u", ColIndex);
		const vector<float> &Freqs = FreqsVec[ColIndex];
		for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float f = Freqs[Letter];
			if (f == 0)
				Log("     .");
			else
				Log(" %5.3f", f);
			}
		string Logo;
		Log("  %s\n", GetLogo(Freqs, Logo).c_str());
		}
	}

void PSSM::LogMe() const
	{
	unsigned ColCount = GetColCount();
	unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();

	Log("%u cols, alpha %s\n", ColCount, m_IsNucleo ? "nt" : "aa");

	if (!m_RawFreqs.empty())
		{
		Log("\n");
		Log("Raw freqs:\n");
		LogFreqs(m_RawFreqs);

		Log("\n");
		Log("Pseudo freqs, weight %.3f:\n", m_PseudoWeight);
		LogFreqs(m_PseudoFreqs);

		Log("\n");
		Log("Freqs:\n");
		LogFreqs(m_Freqs);
		}

	if (m_Scores.empty())
		return;

	Log("\n");
	Log("Scores:\n");
	Log("Col");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" %5c", LetterToChar[Letter]);
	Log(" %5.5s", "Wild");
	Log("\n");

	Log("---");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" -----");
	Log(" -----");
	Log("\n");
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		Log("%3u", ColIndex);
		const vector<float> &Scores = m_Scores[ColIndex];
		for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float s = m_Scores[ColIndex][Letter];
			Log(" %5.2f", s);
			}
		Log(" %5.2f", Scores[AlphaSize]);
		if (m_Freqs.empty())
			Log(" ");
		else
			{
			string Logo;
			Log("  %s ", GetLogo(m_Freqs[ColIndex], Logo).c_str());
			}

		vector<unsigned> Order;
		SortDescending(Scores, Order);

		for (unsigned i = 0; i < AlphaSize; ++i)
			{
			unsigned Letter = Order[i];
			float Score = Scores[Letter];
			if (Score <= 0.0f)
				break;
			Log(" %c(%.2f)", LetterToChar[Letter], Score);
			}
		Log("\n");
		}

	//Log("\n");
	//Log("Training scores:\n");
	//const unsigned SeqCount = SIZE(m_Seqs);
	//vector<float> Scores;
	//for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
	//	{
	//	const char *Seq = m_Seqs[SeqIndex].c_str();
	//	float Score = GetScore(Seq);
	//	Scores.push_back(Score);
	//	}
	//vector<unsigned> Order;
	//SortDescending(Scores, Order);

	//for (unsigned i = 0; i < SeqCount; ++i)
	//	{
	//	unsigned SeqIndex = Order[i];
	//	float Score = Scores[SeqIndex];
	//	const char *Seq = m_Seqs[SeqIndex].c_str();
	//	Log("%8.4f  %s  >%s\n", Score, Seq, m_Labels[SeqIndex].c_str());
	//	}

	//Log("\n");
	//Log("Random scores:\n");
	//string s;
	//for (unsigned i = 0; i < 16; ++i)
	//	{
	//	const char *Seq = GetRandomSeq(s);
	//	float Score = GetScore(Seq);
	//	Log("%8.4f  %s\n", Score, Seq);
	//	}
	}


// ��ȡÿ��λ��ԭʼ���ָ������İ����ᣬ����ַ��������ָ���>=0.5Ϊ��д����֮Сд��
const char *PSSM::CalcConsSeq(string &Seq) const
	{
	Seq.clear();
	uint ColCount = GetColCount();
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		char c = CalcConsChar(ColIndex);
		Seq += c;
		}
	return Seq.c_str();
	}

// ����ԭʼ����Ƶ�����İ����ᣨĳһλ�ã�
char PSSM::CalcConsChar(uint ColIndex) const
	{
	const vector<float> &Freqs = GetFreqs(ColIndex);	// ���ԭʼ����Ƶ��
	const uint N = SIZE(Freqs);
	asserta(N == 20);
	float TopFreq = 0;
	uint TopLetter = 0;
	for (uint i = 0; i < N; ++i)
		{
		if (Freqs[i] > TopFreq)
			{
			TopFreq = Freqs[i];
			TopLetter = i;
			}
		}
	const char *LetterToChar = GetLetterToChar();
	char c = LetterToChar[TopLetter];
	if (TopFreq < 0.5)
		c = tolower(c);
	return c;
	}

// ��ȡ��һλ����صİ�����
char PSSM::GetConsChar(uint ColIndex) const
	{
	asserta(ColIndex < SIZE(m_ConsSeq));
	char c = m_ConsSeq[ColIndex];
	return c;
	}

const string &PSSM::GetLogo(const vector<float> &Freqs, string &Logo, unsigned n) const
	{
	Logo.clear();
	vector<unsigned> Order;
	SortAscending(Freqs, Order);
	unsigned AlphaSize = GetAlphaSize();
	asserta(SIZE(Order) == AlphaSize);
	const char *LetterToChar = GetLetterToChar();
	char c = '.';
	for (unsigned i = 0; i < AlphaSize; ++i)
		{
		unsigned Letter = Order[i];
		float f = Freqs[Letter];
		c = LetterToChar[Letter];
		unsigned w = unsigned(f*n);
		for (unsigned k = 0; k < w; ++k)
			Logo.push_back(c);
		}
	while (SIZE(Logo) < n)
		Logo.push_back(c);
	return Logo;
	}

///  ��� PSSM ��������
void PSSM::SetScores()
	{
	unsigned AlphaSize = GetAlphaSize();
	m_Scores.clear();
	unsigned ColCount = GetColCount();
	m_Scores.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		m_Scores[ColIndex].resize(AlphaSize+1);
		SetScoresCol(ColIndex);
		}
//	SetBackgroundScore();
	}

///  ��� PSSM �������֣�ĳһλ�ã�
void PSSM::SetScoresCol(unsigned ColIndex)
	{
	extern double g_AAFreqs[20];	// ����������Ȼ�������еĳ���Ƶ��
	unsigned AlphaSize = GetAlphaSize();
	float WildScore = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		asserta(ColIndex < SIZE(m_Freqs));
		asserta(ColIndex < SIZE(m_Scores));
		asserta(Letter < SIZE(m_Freqs[ColIndex]));
		asserta(Letter < SIZE(m_Scores[ColIndex]));
		float f = m_Freqs[ColIndex][Letter];
		if (m_IsNucleo)
			{
			double Score = log(f) - log((float) 0.25);
			WildScore += float(Score*0.25);
			m_Scores[ColIndex][Letter] = float(Score);
			}
		else
			{
			double Score = log(f) - log((float) g_AAFreqs[Letter]);		// log(�������ڸñ��������е�Ԥ��Ƶ�� / ����������Ȼ�������еĳ���Ƶ��)
			WildScore += float(Score*g_AAFreqs[Letter]);				// �� log(�������ڸñ��������е�Ԥ��Ƶ�� / ����������Ȼ�������еĳ���Ƶ��) * ����������Ȼ�������еĳ���Ƶ��
			m_Scores[ColIndex][Letter] = float(Score);
			}
		}
	asserta(AlphaSize < SIZE(m_Scores[ColIndex]));
	if (WildScore > 0.0)
		WildScore = 0.0;
	m_Scores[ColIndex][AlphaSize] = WildScore;
	}

// ��ȡ������������ĳһλ�õĵ÷֣����ڵ�ǰPSSM��
float PSSM::GetScore1(uint ColIndex, char c) const
	{
	const byte *CharToLetter = GetCharToLetter();
	if (c == '*')
		return (float) opt_stop_score;
	unsigned AlphaSize = GetAlphaSize();
	unsigned Letter = CharToLetter[c];
	if (Letter == INVALID_LETTER)
		return m_Scores[ColIndex][AlphaSize];
	asserta(Letter < AlphaSize);
	return m_Scores[ColIndex][Letter];
	}

// ��ȡ���еĵ÷֣����ڵ�ǰPSSM��
float PSSM::GetScore(const char *Seq) const
	{
	unsigned AlphaSize = GetAlphaSize();
	const byte *CharToLetter = GetCharToLetter();
	unsigned ColCount = GetColCount();
	float Score = 0.0f;
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		char c = Seq[ColIndex];
		if (c == '*')
			Score += (float) opt_stop_score;
		else
			{
			unsigned Letter = CharToLetter[c];
			if (Letter == INVALID_LETTER)
				Score += m_Scores[ColIndex][AlphaSize];
			else
				{
				asserta(Letter < AlphaSize);
				Score += m_Scores[ColIndex][Letter];
				}
			}
		}
	return Score;
	}

// ������ͬ���ȵ��������
const char *PSSM::GetRandomSeq(string &s) const
	{
	s.clear();
	unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();
	unsigned ColCount = GetColCount();
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		byte Letter = GetRandomLetter();
		char c = LetterToChar[Letter];
		s.push_back(c);
		}
	return s.c_str();
	}

// ������ĸ��Χ�ڵ��������
byte PSSM::GetRandomLetter() const
	{
	unsigned AlphaSize = GetAlphaSize();
	return rand()%AlphaSize;
	}


// �� PSSM ����д���ļ��������ļ����
void PSSM::ToFile(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	ToFile(f);
	CloseStdioFile(f);
	}

// �� PSSM score���� ��Freq���󣨸�λ�ý��򣩡�score���󣨸�λ�ý���д���ļ�
void PSSM::ToFile(FILE *f) const
	{
	unsigned ColCount = GetColCount();
	unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();

	fprintf(f, "PSSM %s %u\n", m_IsNucleo ? "nt" : "aa", ColCount);
	fprintf(f, "consseq	%s\n", m_ConsSeq.c_str());
	if (m_IsNucleo)
		fprintf(f, "Col     A     C     G     T     N\n");
	else
		fprintf(f, "Col     A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y     X\n");

	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		// ��ӡ score ����
		fprintf(f, "%3u", ColIndex + 1);
		const vector<float> &Scores = m_Scores[ColIndex];
		for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float s = m_Scores[ColIndex][Letter];
			fprintf(f, " %5.2f", s);
			}
		fprintf(f, " %5.2f", Scores[AlphaSize]);

		// ��ӡ Freq ���󣨽���
		if (!m_Freqs.empty())
			{
			string Logo;
			fprintf(f, "  %s ", GetLogo(m_Freqs[ColIndex], Logo).c_str());
			}

		vector<unsigned> Order;
		SortDescending(Scores, Order);

		for (unsigned i = 0; i < AlphaSize; ++i)
			{
			unsigned Letter = Order[i];
			float Score = Scores[Letter];
			if (Score <= 0.0f)
				break;
			fprintf(f, " %c(%.2f)", LetterToChar[Letter], Score);
			}
		fprintf(f, "\n");
		}
	fprintf(f, "//\n");
	}



// ��ȡ PSSM �ļ�
void PSSM::FromFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	FromFile(f);
	CloseStdioFile(f);
	}

void PSSM::FromFile(FILE *f)
	{
	vector<string> Strings;
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		Strings.push_back(Line);	// �� PSSM ������Ϣ�ڶ��п�ʼ��ȡ������ String
		if (Line == "//")			// '//' ֹͣ���ж�Ϊһ�������ȡ���
			break;
		}
	FromStrings("FILE*", Strings);	// ���� String ����
	}


static bool GetNextString(const vector<string>& Strings, const unsigned N,
	unsigned& i, string& s)
{
	if (i >= N)
		return false;
	s = Strings[i++];
	return true;
}

void PSSM::FromCStrings(const string &Name, const char **CStrings)
	{
	vector<string> Strings;
	for (unsigned i = 0; i < 32; ++i)
		{
		const string s = string(CStrings[i]);
		Strings.push_back(s);
		asserta(!s.empty());
		if (s[0] == '/')
			{
			FromStrings(Name, Strings);
			return;
			}
		}
	Die("Missing // in psm %s", Name.c_str());
	}

void PSSM::FromStrings(const string &Name, const vector<string> &Strings)
	{
	Clear();

	const unsigned N = SIZE(Strings);
	string Line;
	unsigned LineNr = 0;

// �����һ��' PSSM aa 9\n '
	bool Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok)
		Die("PSSM::FromFile, empty .psm %s", Name.c_str());

	vector<string> Fields;
	Split(Line, Fields, ' ');		// ���ַ��� Line �� ' ' �ָ�����ַ����� Fields ��
	if (SIZE(Fields) != 3 || Fields[0] != "PSSM")
		Die("Invalid header in .psm %s", Name.c_str());

	unsigned AlphaSize = 0;
	if (Fields[1] == "nt")
		{
		m_IsNucleo = true;
		AlphaSize = 4;
		}
	else if (Fields[1] == "aa")
		{
		m_IsNucleo = false;
		AlphaSize = 20;
		}
	else
		Die("Invalid alphabet in .psm %s", Name.c_str());
// //

	m_ColCount = atou(Fields[2]);	//	��ȡ���г���
	if (m_ColCount == 0)
		Die("Zero cols in .psm file %s", Name.c_str());

// ����ڶ���' consseq\tAbcDEFghIjkLMN\n '
	Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok)
		Die("Premature end-of-file in .psm %s", Name.c_str());

	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2)
		Die("Bad consseq line %s", Line.c_str());
	m_ConsSeq = Fields[1];

// ���������	' Col     A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y     X\n '
	Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok)
		Die("Premature end-of-file in .psm %s", Name.c_str());

	if (!StartsWith(Line, "Col "))
		Die("Invalid 2nd header line in .psm %s", Name.c_str());

// �������������'   1   0.11  0.15  0.09  0.31  ...   0.02\n '
//				    2   0.01  0.35  0.13  0.04  ...   0.16\n '
	m_Scores.resize(m_ColCount);
	for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		{
		m_Scores[ColIndex].resize(AlphaSize+1, 0.0f);
		Ok = GetNextString(Strings, N, LineNr, Line);
		if (!Ok)
			Die("Premature end-of-file in .psm %s", Name.c_str());
		Split(Line, Fields, 0);
		if (SIZE(Fields) < AlphaSize + 2 || atou(Fields[0]) != ColIndex + 1)
			Die("Invalid scores for col %u in .psm  %s", ColIndex, Name.c_str());

		for (unsigned Letter = 0; Letter <= AlphaSize; ++Letter)
			{
			const string &f = Fields[Letter+1];
			float Score = (float) StrToFloat(f.c_str());
			m_Scores[ColIndex][Letter] = Score;
			}
		}
	Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok || Line != "//")
		Die("Missing // at end of .psm %s", Name.c_str());
	}

