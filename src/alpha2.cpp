#include "allocid.h"
static ALLOCID AllocId = ALLOCID_alpha2;
#include "myutils.h"
#include "alpha.h"
#include "timing.h"

/***
IUPAC nucleotide ambiguity codes

------------------------------------------
Symbol       Meaning      Nucleic Acid
------------------------------------------
A            A           Adenine
C            C           Cytosine
G            G           Guanine
T            T           Thymine
U            U           Uracil
M          A or C
R          A or G
W          A or T
S          C or G
Y          C or T
K          G or T
V        A or C or G
H        A or C or T
D        A or G or T
B        C or G or T
X      G or A or T or C
N      G or A or T or C
	

Reference:
IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030. 
***/

struct IUPAC_Code
	{
	char Code;
	const char *Chars;
	char CompCode;
	};

static IUPAC_Code IUPAC_Codes[] =
	{
//	Code	Means				Comp	CompCode
	{ 'M', "AC",   'K' },		// GT		K
	{ 'R', "AG",   'Y' },		// CT		Y
	{ 'W', "AT",   'W' },		// AT		W
	{ 'S', "CG",   'S' },		// CG		S
	{ 'Y', "CT",   'R' },		// AG		R
	{ 'K', "GT",   'M' },		// AC		M
	{ 'V', "ACG",  'B' },		// CGT		B
	{ 'H', "ACT",  'D' },		// AGT		D
	{ 'D', "AGT",  'H' },		// ACT		H
	{ 'B', "CGT",  'V' },		// ACG		V
	{ 'X', "GATC", 'X' },		// ACGT		X
	{ 'N', "GATC", 'N' },		// ACGT		N
	};

// IUPAC 编码数量
static const unsigned IUPAC_Code_Count = sizeof(IUPAC_Codes)/sizeof(IUPAC_Codes[0]);

static byte BIT_A = 0x1;
static byte BIT_C = 0x2;
static byte BIT_G = 0x4;
static byte BIT_T = 0x8;

// 核苷酸字符转2进制查询表
byte g_Nucleo_CharToBit[256];
// IUPAC编码的核苷酸字符转2进制查询表（模棱两可的核苷酸）
byte g_IUPAC_CharToBits[256];
// IUPAC编码的核苷酸字符转反义核苷酸查询表（模棱两可的核苷酸）
byte g_IUPAC_CharToCompChar[256];


/***
Code	    Bits	  CompCode
{ 'A',	"00000001",   'A' },	
{ 'C',	"00000010",   'C' },
{ 'G',	"00000100",   'G' },
{ 'T',	"00001000",   'T' },
{ 'M',	"00000011":"AC",   'K' },	
{ 'R',	"00000101":"AG",   'Y' },
{ 'W',	"00001001":"AT",   'W' },
{ 'S',	"00000110":"CG",   'S' },
{ 'Y',	"00001010":"CT",   'R' },
{ 'K',	"00001100":"GT",   'M' },

... ...

***/

// 初始化 IUPAC 编码的反义核苷酸字符查询表 和 二进制编号查询表
static bool Init_IUPAC()
	{
	for (unsigned k = 0; k < 256; ++k)
		g_IUPAC_CharToCompChar[k] = k;

	g_Nucleo_CharToBit['a'] = BIT_A;
	g_Nucleo_CharToBit['A'] = BIT_A;

	g_Nucleo_CharToBit['c'] = BIT_C;
	g_Nucleo_CharToBit['C'] = BIT_C;

	g_Nucleo_CharToBit['g'] = BIT_G;
	g_Nucleo_CharToBit['G'] = BIT_G;

	g_Nucleo_CharToBit['t'] = BIT_T;
	g_Nucleo_CharToBit['T'] = BIT_T;

	g_Nucleo_CharToBit['u'] = BIT_T;
	g_Nucleo_CharToBit['U'] = BIT_T;

	g_IUPAC_CharToBits['a'] = BIT_A;
	g_IUPAC_CharToBits['A'] = BIT_A;

	g_IUPAC_CharToBits['c'] = BIT_C;
	g_IUPAC_CharToBits['C'] = BIT_C;

	g_IUPAC_CharToBits['g'] = BIT_G;
	g_IUPAC_CharToBits['G'] = BIT_G;

	g_IUPAC_CharToBits['t'] = BIT_T;
	g_IUPAC_CharToBits['T'] = BIT_T;

	g_IUPAC_CharToBits['u'] = BIT_T;
	g_IUPAC_CharToBits['U'] = BIT_T;

	g_IUPAC_CharToCompChar['A'] = 'T';
	g_IUPAC_CharToCompChar['a'] = 't';

	g_IUPAC_CharToCompChar['C'] = 'G';
	g_IUPAC_CharToCompChar['c'] = 'g';

	g_IUPAC_CharToCompChar['G'] = 'C';
	g_IUPAC_CharToCompChar['g'] = 'c';

	g_IUPAC_CharToCompChar['T'] = 'A';
	g_IUPAC_CharToCompChar['t'] = 'a';

	g_IUPAC_CharToCompChar['U'] = 'A';
	g_IUPAC_CharToCompChar['u'] = 'a';

	for (unsigned k = 0; k < IUPAC_Code_Count; ++k)
		{
		IUPAC_Code &IC = IUPAC_Codes[k];
		byte Code = IC.Code;
		byte CompCode = IC.CompCode;
		g_IUPAC_CharToCompChar[Code] = CompCode;

		byte Bits = 0;
		for (const char *p = IC.Chars; *p; ++p)
			{
			byte Char = (byte) *p;
			byte Bit = g_Nucleo_CharToBit[Char];
			Bits |= Bit;
			}
		g_IUPAC_CharToBits[Code] = Bits;
		}

#if	DEBUG
	for (unsigned k = 0; k < 256; ++k)
		assert(g_IUPAC_CharToCompChar[g_IUPAC_CharToCompChar[k]] == k || toupper(k) == 'U');

	assert( IUPAC_Eq('a', 'A'));
	assert( IUPAC_Eq('A', 'A'));

	assert( IUPAC_Eq('A', 'M'));
	assert( IUPAC_Eq('C', 'M'));
	assert(!IUPAC_Eq('G', 'M'));
	assert(!IUPAC_Eq('T', 'M'));

	assert( IUPAC_Eq('A', 'V'));
	assert( IUPAC_Eq('C', 'V'));
	assert( IUPAC_Eq('G', 'V'));
	assert(!IUPAC_Eq('T', 'V'));

	assert( IUPAC_Eq('A', 'N'));
	assert( IUPAC_Eq('C', 'N'));
	assert( IUPAC_Eq('G', 'N'));
	assert( IUPAC_Eq('T', 'N'));
#endif
	return true;
	}

// 初始化 g_IUPAC_CharToCompChar[256]
static bool g_IUPAC_CharToBitsInitDone = Init_IUPAC();


const char *WordToStrAmino(unsigned Word, unsigned WordLength)
	{
	static char Str[32];
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%20;
		Str[WordLength-i-1] = g_LetterToCharAmino[Letter];
		Word /= 20;
		}
	Str[WordLength] = 0;
	return Str;
	}

// (((20 + aa1) * 20 + aa2) * 20 + aa3) * 20 + aa4 ......
// 信息转换，二进制表示的氨基酸序列转氨基酸字符串
const char *WordToStrAmino2(unsigned Word, unsigned WordLength, char *Str)
	{
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%20;
		Str[WordLength-i-1] = g_LetterToCharAmino[Letter];
		Word /= 20;
		}
	Str[WordLength] = 0;
	return Str;
	}

// (((20 + aa1) * 20 + aa2) * 20 + aa3) * 20 + aa4 ......
// 信息转换，二进制表示的核酸序列转核酸字符串
const char *WordToStrNucleo(unsigned Word, unsigned WordLength)
	{
	static char Str[32];
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%4;
		Str[WordLength-i-1] = g_LetterToCharNucleo[Letter];
		Word /= 4;
		}
	Str[WordLength] = 0;
	return Str;
	}

// 信息转换，二进制表示的核酸/氨基酸序列转核酸/氨基酸字符串
const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo)
	{
	return (Nucleo ? WordToStrNucleo : WordToStrAmino)(Word, WordLength);
	}


// 将核酸序列转成反义链（3' 和 5'反转 & 配位相反[A-T;C-G]）
byte *RevCompAlloc(const byte *Seq, unsigned L)
	{
	byte *RCSeq = myalloc(byte, L);

	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_CharToCompChar[Seq[i]];

	return RCSeq;
	}

// 将核酸序列转成反义链（3' 和 5'反转 & 配位相反[A-T;C-G]）
void RevCompSeq(string &Seq)
	{
	unsigned L = SIZE(Seq);
	unsigned L1 = L - 1;
	unsigned L2 = L/2;
	for (unsigned i = 0; i < L2; ++i)
		{
		unsigned j = L1 - i;
		unsigned ci = Seq[i];
		unsigned cj = Seq[j];

		char ri = g_CharToCompChar[ci];
		char rj = g_CharToCompChar[cj];

		Seq[i] = rj;
		Seq[j] = ri;
		}

	if (L%2 == 1)
		Seq[L2] = g_CharToCompChar[Seq[L2]];
	}

// 将核酸序列转成反义链（ 3' 和 5'反转 & 配位相反[A-T;C-G] ）(In place)
void RevCompInPlace(byte *Seq, unsigned L)
	{
	unsigned L1 = L - 1;
	unsigned L2 = L/2;
	for (unsigned i = 0; i < L2; ++i)
		{
		unsigned j = L1 - i;
		unsigned ci = Seq[i];
		unsigned cj = Seq[j];

		unsigned ri = g_CharToCompChar[ci];
		unsigned rj = g_CharToCompChar[cj];

		Seq[i] = rj;
		Seq[j] = ri;
		}

	if (L%2 == 1)
		Seq[L2] = g_CharToCompChar[Seq[L2]];
	}

// 核酸序列转反义链
void RevComp(const byte *Seq, unsigned L, byte *RCSeq)
	{
	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_CharToCompChar[Seq[i]];
	}

// 核酸序列转反义链（IUPAC 编码）
void RevCompIUPAC(const byte *Seq, unsigned L, byte *RCSeq)
	{
	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_IUPAC_CharToCompChar[Seq[i]];
	}

// codon 转氨基酸序列，单字符处理
unsigned char GetAminoCharFrom3NucChars(unsigned char c1, unsigned char c2,
  unsigned char c3)
	{
	unsigned Letter1 = g_CharToLetterNucleo[c1];
	unsigned Letter2 = g_CharToLetterNucleo[c2];
	unsigned Letter3 = g_CharToLetterNucleo[c3];
	unsigned Word = Letter1*(4*4) + Letter2*4 + Letter3;

	unsigned Letter = g_CodonWordToAminoLetter[Word];
	return g_LetterToCharAmino[Letter];
	}
