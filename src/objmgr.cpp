#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "pathinfo.h"
//#include "alignresult.h"
#include "timing.h"

#if	TRACK_OBJ_THREAD
#include <omp.h>
#endif

#undef Up
#undef Down

const char *ObjTypeToStr(ObjType Type)
	{
	switch (Type)
		{
#define T(x)	case OT_##x: return #x;
#include "objtypes.h"
///	case OT_SeqInfo: return "SeqInfo"
///	case OT_PathInfo: return "PathInfo"
		}
	return "OT_??";
	}

const char *ObjTypeToStr2(ObjType Type)
	{
	switch (Type)
		{
	case OT_SeqInfo:
		return "SI";
	case OT_PathInfo:
		return "PI";
	//case OT_AlignResult:
	//	return "AR";
		}
	return "??";
	}

ObjMgr::ObjMgr()
	{
	zero(m_Free, OTCount);
	zero(m_Busy, OTCount);

#if	DEBUG
	m_Validate = false;
	zero(m_BusyCounts, OTCount);
	zero(m_GetCallCounts, OTCount);
	zero(m_AllocCallCounts, OTCount);
	zero(m_FreeCallCounts, OTCount);
#endif
	}

Obj *ObjMgr::AllocNew(ObjType Type)
	{
#if	DEBUG
	++(m_AllocCallCounts[Type]);
#endif
	switch (Type)
		{

#define	T(x)	case OT_##x: return new x(this);
#include "objtypes.h"

//		case OT_SeqInfo: return new SeqInfo(this);
//		case OT_PathInfo: return new SeqInfo(this);

	default:
			assert(false);
		}
	return 0;
	}


// 从 Free 栈中获取一个对象添加到 Busy 栈中，并返回该对象的地址
Obj *ObjMgr::GetObj(ObjType Type)
	{
#if	DEBUG
	++(m_GetCallCounts[Type]);
	assert(Type < OTCount);
#endif
	Obj *NewObj = 0;
	if (m_Free[Type] == 0)
		NewObj = AllocNew(Type);			// 若 SeqInfo 类型 Free 栈空了则创建一个新的 SeqInfo 对象
	else
		{
		NewObj = m_Free[Type];				// 若非空，则获得栈中第一个 SeqInfo 对象

		assert(NewObj->m_RefCount == 0);	// Free 栈中的对象引用数一定为 0 
		m_Free[Type] = m_Free[Type]->m_Fwd;	// 将栈中下一个对象赋值到第一个对象位置，并清除对前一对线的链接
		if (m_Free[Type])
			m_Free[Type]->m_Bwd = 0;		
		}

	if (m_Busy[Type] != 0)					// 若该对象类型在 Busy 栈中非空
		{
		assert(m_Busy[Type]->m_Bwd == 0);	// 栈中第一个对象的前一对象为空
		m_Busy[Type]->m_Bwd = NewObj;	
		}
	NewObj->m_Fwd = m_Busy[Type];			// 将 busy 栈中第一个对象设置为新获得的对象，并处理链接关系
	m_Busy[Type] = NewObj;

	assert(NewObj != 0);
	NewObj->m_RefCount = 1;					// 新对象的引用数赋值为 1
#if	TRACK_OBJ_THREAD
	NewObj->m_OMPThreadIndex = omp_get_thread_num();
#endif

#if	DEBUG
	++(m_BusyCounts[Type]);
	if (m_Validate)
		Validate();
#endif

	return NewObj;
	}

void ObjMgr::FreeObj(Obj *obj)
	{
	assert(obj->m_RefCount == 0);

	ObjType Type = obj->m_Type;
	assert(Type < OTCount);

#if	DEBUG
	++(m_FreeCallCounts[Type]);
#endif

	if (obj == m_Busy[Type])
		m_Busy[Type] = obj->m_Fwd;

	Obj *Prev = obj->m_Bwd;
	Obj *Next = obj->m_Fwd;

	if (Prev != 0)
		Prev->m_Fwd = Next;
	if (Next != 0)
		Next->m_Bwd = Prev;

	if (m_Free[Type] != 0)
		{
		assert(m_Free[Type]->m_Bwd == 0);
		m_Free[Type]->m_Bwd = obj;
		}
	obj->m_Fwd = m_Free[Type];
	obj->m_Bwd = 0;
	m_Free[Type] = obj;

#if	DEBUG
	if (m_Validate)
		Validate();
#endif
	}

// 统计 free 链表的节点数
unsigned ObjMgr::GetFreeCount(ObjType Type) const
	{
	assert(Type < OTCount);
	unsigned Length = 0;
	const Obj *obj = m_Free[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		++Length;
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return Length;
	}

// 统计 Busy 链表的节点数
unsigned ObjMgr::GetBusyCount(ObjType Type) const
	{
	assert(Type < OTCount);
	unsigned Length = 0;
	const Obj *obj = m_Busy[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		++Length;
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return Length;
	}

// 获得 Obj 中最大的引用数量
unsigned ObjMgr::GetMaxRefCount(ObjType Type) const
	{
	assert(Type < OTCount);
	unsigned MaxRefCount = 0;
	const Obj *obj = m_Busy[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		if (obj->m_RefCount > MaxRefCount)
			MaxRefCount = obj->m_RefCount;
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return MaxRefCount;
	}

// 获得 Busy 链表中对象的总内存分配
float ObjMgr::GetTotalMem(ObjType Type) const
	{
	assert(Type < OTCount);
	float Total = 0.0f;
	const Obj *obj = m_Busy[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		Total += obj->GetMemBytes();
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return Total;
	}

#if	DEBUG
void ObjMgr::ValidateType(ObjType Type) const
	{
	unsigned NA = m_AllocCallCounts[Type];
	unsigned NF = m_FreeCallCounts[Type];

	unsigned nb = 0;
	for (const Obj *obj = m_Busy[Type]; obj; obj = obj->m_Fwd)
		{
		++nb;
		assert(nb <= NA);
		assert(obj->m_Type == Type);
		assert(obj->m_RefCount > 0);
		if (obj->m_Bwd)
			assert(obj->m_Bwd->m_Fwd == obj);
		if (obj->m_Fwd)
			assert(obj->m_Fwd->m_Bwd == obj);
		}

	unsigned nf = 0;
	for (const Obj *obj = m_Free[Type]; obj; obj = obj->m_Fwd)
		{
		++nf;
		assert(nf <= NF);
		assert(obj->m_RefCount == 0);
		assert(obj->m_Type == Type);
		if (obj->m_Bwd)
			assert(obj->m_Bwd->m_Fwd == obj);
		if (obj->m_Fwd)
			assert(obj->m_Fwd->m_Bwd == obj);
		}
	assert(nb + nf == NA);
	assert(nb == m_BusyCounts[Type]);
	}

void ObjMgr::Validate() const
	{
	Die("Validate!");
	for (unsigned i = 0; i < OTCount; ++i)
		{
		ObjType Type = (ObjType) i;
		ValidateType(Type);
		}
	}
#endif // DEBUG

#if	TRACE_OBJS

static bool g_HdrDone = false;

static void LogHdr()
	{
	if (g_HdrDone)
		return;

	Log("\n");
	Log("Op        UID  Ty  RefCnt  BsyCnt  Src\n");
	Log("---  --------  --  ------  ------  ---\n");
	g_HdrDone = true;
	}

Obj *ObjMgr::GetObj_(ObjType Type, const char *FileName, unsigned LineNr)
	{
	LogHdr();
	Obj *pObj = GetObj(Type);
	static unsigned UID;
	pObj->m_UID = UID++;
	pObj->m_SourceFileName = FileName;
	pObj->m_SourceLineNr = LineNr;
	Log("Get  %8u  %2.2s  %6u  %6u  %s(%u)\n",
		pObj->m_UID,
		ObjTypeToStr2(pObj->m_Type),
		pObj->GetRefCount(),
		m_BusyCounts[pObj->m_Type],
		BaseName(FileName),
		LineNr);
	return pObj;
	}

void ObjMgr::Up_(Obj *pObj, const char *FileName, unsigned LineNr)
	{
	LogHdr();
	Up(pObj);
	Log("Up   %8u  %2.2s  %6u  %6u  %s(%u)\n",
		pObj->m_UID,
		ObjTypeToStr2(pObj->m_Type),
		pObj->GetRefCount(),
		m_BusyCounts[pObj->m_Type],
		BaseName(FileName),
		LineNr);
	}

void ObjMgr::Down_(Obj *pObj, const char *FileName, unsigned LineNr)
	{
	LogHdr();
	Down(pObj);
	Log("Dwn  %8u  %2.2s  %6u  %6u  %s(%u)\n",
		pObj->m_UID,
		ObjTypeToStr2(pObj->m_Type),
		pObj->GetRefCount(),
		m_BusyCounts[pObj->m_Type],
		BaseName(FileName),
		LineNr);
	}

void ObjMgr::LogBusy() const
	{
	Log("\n");
	Log("Busy:\n");
	Log("     UID  Ty  RefCnt  BsyCnt  Src\n");
	Log("--------  --  ------  ------  ---\n");
	for (unsigned i = 0; i < OTCount; ++i)
		{
		for (const Obj *pObj = m_Busy[i]; pObj; pObj = pObj->m_Fwd)
			Log("%8u  %2.2s  %6u  %s(%u)\n",
				pObj->m_UID,
				ObjTypeToStr2(pObj->m_Type),
				pObj->GetRefCount(),
				pObj->m_SourceFileName,
				pObj->m_SourceLineNr);
		}
	}

#endif // TRACE_OBJS

static unsigned g_FreeCount[OTCount];
static unsigned g_BusyCount[OTCount];
static unsigned g_MaxRefCount[OTCount];
static float g_Mem[OTCount];

#if DEBUG
static unsigned g_GetCallCount[OTCount];
static unsigned g_AllocCallCount[OTCount];
static unsigned g_FreeCallCount[OTCount];
#endif

void ObjMgr::UpdateGlobalStats() const
	{
	for (unsigned i = 0; i < OTCount; ++i)
		{
		ObjType Type = (ObjType) i;
		g_FreeCount[Type] += GetFreeCount(Type);
		g_BusyCount[Type] += GetBusyCount(Type);
		g_MaxRefCount[Type] = max(g_MaxRefCount[Type], GetMaxRefCount(Type));
		g_Mem[Type] += GetTotalMem(Type);
#if	DEBUG
		g_GetCallCount[Type] += m_GetCallCounts[Type];
		g_AllocCallCount[Type] += m_AllocCallCounts[Type];
		g_FreeCallCount[Type] += m_FreeCallCounts[Type];
#endif
		}
	}

void ObjMgr::LogGlobalStats()
	{
	Log("\n");
	Log("            Type        Busy        Free         Mem   MaxRefCnt        Gets      Allocs       Frees\n");
	Log("----------------  ----------  ----------  ----------  ----------  ----------  ----------  ----------\n");
	for (unsigned i = 0; i < OTCount; ++i)
		{
		ObjType Type = (ObjType) i;
		const char *Name = ObjTypeToStr(Type);

		Log("%16.16s", Name);
		Log("  %10.10s", IntToStr(g_BusyCount[Type]));
		Log("  %10.10s", IntToStr(g_FreeCount[Type]));
		Log("  %10.10s", MemBytesToStr(g_Mem[Type]));
		Log("  %10u", g_MaxRefCount[Type]);
#if	DEBUG
		Log("  %10.10s", IntToStr(g_GetCallCount[Type]));
		Log("  %10.10s", IntToStr(g_AllocCallCount[Type]));
		Log("  %10.10s", IntToStr(g_FreeCallCount[Type]));
#endif
		Log("\n");
		}
	}

void ObjMgr::LogMe() const
	{
	Log("\n");
	Log("            Type        Busy        Free         Mem   MaxRefCnt        Gets      Allocs       Frees\n");
	Log("----------------  ----------  ----------  ----------  ----------  ----------  ----------  ----------\n");
	for (unsigned i = 0; i < OTCount; ++i)
		{
		ObjType Type = (ObjType) i;					// Type <- 0, 1
		const char* Name = ObjTypeToStr(Type);		// Name <- 'OT_SeqInfo', 'OT_PathInfo'

		Log("%16.16s", Name);
		Log("  %10.10s", IntToStr(GetBusyCount(Type)));
		Log("  %10.10s", IntToStr(GetFreeCount(Type)));
		Log("  %10.10s", MemBytesToStr(GetTotalMem(Type)));
		Log("  %10u", GetMaxRefCount(Type));
#if	DEBUG
		Log("  %10.10s", IntToStr(m_GetCallCounts[Type]));
		Log("  %10.10s", IntToStr(m_AllocCallCounts[Type]));
		Log("  %10.10s", IntToStr(m_FreeCallCounts[Type]));
#endif
		Log("\n");
		}
#if	TRACE_OBJS
	LogBusy();
#endif
	}
