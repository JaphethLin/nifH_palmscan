#ifndef MY_VERSION
#define MY_VERSION	"24.1.15"
#endif

STR_OPT(OTHER, build_pssm)
STR_OPT(OTHER, search_aa_top)
STR_OPT(OTHER, search_nt_top)
STR_OPT(OTHER, search_nt_all)
STR_OPT(OTHER, search_nt_topx)
STR_OPT(OTHER, search_pair_nt)
STR_OPT(OTHER, search_pair_g)
STR_OPT(OTHER, search_pair_aa)
STR_OPT(OTHER, build_rdrp_model)
STR_OPT(OTHER, otumaps)
STR_OPT(OTHER, alnqc)

STR_OPT(OTHER, report)
STR_OPT(OTHER, fastaout)
STR_OPT(OTHER, alnout)
STR_OPT(OTHER, tsvout)
STR_OPT(OTHER, chain)
STR_OPT(OTHER, db)
STR_OPT(OTHER, freqsout)
STR_OPT(OTHER, relabel)
STR_OPT(OTHER, search_pp)
STR_OPT(OTHER, fevout)
STR_OPT(OTHER, bedout)
STR_OPT(OTHER, model)
STR_OPT(OTHER, ppout)
STR_OPT(OTHER, ppout_nt)
STR_OPT(OTHER, motifs_fastaout)
STR_OPT(OTHER, motifs_fastaout2)
STR_OPT(OTHER, nohit_fastaout)
STR_OPT(OTHER, pssm_alnout)

STR_OPT(OTHER, log)
STR_OPT(OTHER, input)
STR_OPT(OTHER, output)
STR_OPT(OTHER, output2)
STR_OPT(OTHER, output3)
STR_OPT(OTHER, psm)
STR_OPT(OTHER, psm1)
STR_OPT(OTHER, psm2)
STR_OPT(OTHER, spec)

STR_OPT(OTHER, fw_name)
STR_OPT(OTHER, psm_fw)

STR_OPT(OTHER, gapopen)
STR_OPT(OTHER, gapext)
STR_OPT(OTHER, frame)
STR_OPT(OTHER, test)

UNS_OPT(OTHER,	threads,			8,			0,			UINT_MAX)
UNS_OPT(OTHER,	band,				16,			0,			UINT_MAX)
UNS_OPT(OTHER,	hspw,				0,			1,			UINT_MAX)
UNS_OPT(OTHER,	minhsp,				32,			1,			UINT_MAX)
UNS_OPT(OTHER,	iddrop,				8,			1,			UINT_MAX)
UNS_OPT(OTHER,	cluster_maxdiffs,	1,			0,			UINT_MAX)
UNS_OPT(OTHER,	topn,				32,			1,			UINT_MAX)
UNS_OPT(OTHER,	maxx,				10,			1,			UINT_MAX)

UNS_OPT(OTHER,	secs,				60,			1,			UINT_MAX)

UNS_OPT(OTHER,	maxseqlength,		500000000,		1,			UINT_MAX)
UNS_OPT(OTHER,	sfasta_buff_bytes,	512*1024*1024,1024,		UINT_MAX)
UNS_OPT(OTHER,	randseed,			0,			0,			UINT_MAX)
UNS_OPT(OTHER,	usort_w,			6,			1,			UINT_MAX)
UNS_OPT(OTHER,	mingap,				0,			1,			UINT_MAX)
UNS_OPT(OTHER,	maxgap,				999,		1,			UINT_MAX)
UNS_OPT(OTHER,	hiw,				45,			1,			UINT_MAX)
UNS_OPT(OTHER,	low,				16,			1,			UINT_MAX)

UNS_OPT(OTHER,	fastq_truncqual,	UINT_MAX,	0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_minlen,		0,			0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_minovlen,		0,			0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_trunclen,		0,			0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_maxdiffs,		UINT_MAX,	0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_ascii,		33,			0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_qmin,			0,			0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_qmax,			41,			0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_qmaxout,		41,			0,			UINT_MAX)
UNS_OPT(OTHER,	fastq_stripleft,	0,			0,			UINT_MAX)

FLT_OPT(OTHER,	minscore,			0.0,		-9e9,		+9e9)
FLT_OPT(OTHER,	minscore1,			0.0,		-9e9,		+9e9)
FLT_OPT(OTHER,	minscore2,			0.0,		-9e9,		+9e9)
FLT_OPT(OTHER,	minscore_pair,		0.0,		-9e9,		+9e9)
FLT_OPT(OTHER,	stop_score,			-10.0,		-9e9,		+9e9)
FLT_OPT(OTHER,	minscore_fw,		0.0,		-9e9,		+9e9)
FLT_OPT(OTHER,	fastq_maxee,		1.0,		0.0,		9e9)

FLT_OPT(OTHER,	match,				1.0,		0.0,		DBL_MAX)
FLT_OPT(OTHER,	mismatch,			-2.0,		0.0,		DBL_MAX)
FLT_OPT(OTHER,	xdrop_u,			16.0,		0.0,		DBL_MAX)
FLT_OPT(OTHER,	xdrop_g,			32.0,		0.0,		DBL_MAX)
FLT_OPT(OTHER,	xdrop_nw,			16.0,		0.0,		DBL_MAX)
FLT_OPT(APL,	lopen,				10.0,		0.0,		DBL_MAX)
FLT_OPT(APL,	lext,				1.0,		0.0,		DBL_MAX)

FLAG_OPT(OTHER, trunclabels)
FLAG_OPT(OTHER, notrunclabels)
FLAG_OPT(OTHER, compilerinfo)
FLAG_OPT(OTHER, quiet)
FLAG_OPT(OTHER, logmemgrows)
FLAG_OPT(OTHER, fulldp)
FLAG_OPT(OTHER, logquery)
FLAG_OPT(OTHER, verbose)
FLAG_OPT(OTHER, all)
FLAG_OPT(OTHER, rt)
FLAG_OPT(OTHER, rdrp)
FLAG_OPT(OTHER, nifH)
FLAG_OPT(OTHER, hiconf)
FLAG_OPT(OTHER, loconf)
FLAG_OPT(OTHER, gapped)
FLAG_OPT(OTHER, coords)

#undef FLAG_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT
