#ifndef objmgr_h
#define objmgr_h

#include "objtype.h"
/***
enum ObjType
{
	OT_SeqInfo = 0
	OT_PathInfo = 1
	OTCount = 2
};
***/ 

#include "obj.h"
class Obj;

#define T(x)	class x;
#include "objtypes.h"
/// class(SeqInfo)
/// class(PathInfo)

const char *ObjTypeToStr(ObjType Type);

// ��������� ��ÿ���߳�һ�� ������ free �� busy �� Obj ����˫������ ����ʵ�ֶԸ��� Obj ����������������õ����� ��ʹ���ڴ��������ʵʱ���

class ObjMgr
	{
	friend class Obj;

private:
	Obj *m_Free[OTCount];
	Obj *m_Busy[OTCount];

#if	DEBUG
	bool m_Validate;
	unsigned m_BusyCounts[OTCount];
	unsigned m_GetCallCounts[OTCount];
	unsigned m_FreeCallCounts[OTCount];
	unsigned m_AllocCallCounts[OTCount];
#endif

public:
	ObjMgr();
	void LogMe() const;
	void UpdateGlobalStats() const;
	static void LogGlobalStats();
	Obj *GetObj(ObjType Type);
	Obj *GetObj_(ObjType Type, const char *FileName, unsigned LineNr);

#define T(x)	\
		x *Get##x() { return (x *) GetObj(OT_##x); } \
		x *Get##x##_(const char *FileName, unsigned LineNr) { return (x *) GetObj_(OT_##x, FileName, LineNr); }
#include "objtypes.h"
/**

	SeqInfo *GetSeqInfo()
	{ return (SeqInfo *) GetObj(OT_SeqInfo); }

	PathInfo *GetPathInfo_(const char *FileName, unsigned LineNr)
	{ return (PathInfo *) GetObj_(OT_PathInfo, FileName, LineNr); }

**/

	void Up(Obj *pObj)		// Obj ������������ +1
		{
#if	DEBUG
		if (m_Validate)
			Validate();
#endif

		assert(pObj->m_RefCount != 0);
		++pObj->m_RefCount;

#if	DEBUG
		if (m_Validate)
			Validate();
#endif
		}


	void Down(Obj *pObj)	// Obj ������������ -1������������Ϊ 0 ���ͷŶ����ڴ�
		{
#if	DEBUG
		if (m_Validate)
			Validate();
#endif
#if	TRACK_OBJ_THREAD
		asserta(pObj->m_OMPThreadIndex == omp_get_thread_num());
#endif
		assert(pObj->m_RefCount > 0);
		--pObj->m_RefCount;
		if (pObj->m_RefCount == 0)
			{
			ObjType Type = pObj->m_Type;
#if	DEBUG
			assert(m_BusyCounts[Type] > 0);
			--(m_BusyCounts[Type]);
#endif
			FreeObj(pObj);
			pObj->OnZeroRefCount();
			}

#if	DEBUG
		if (m_Validate)
			Validate();
#endif
		}


	/**
		void Up(Obj* pObj)		
		{
			assert(pObj->m_RefCount != 0);
			++pObj->m_RefCount;
		}

		void Down(Obj *pObj)
		{

			assert(pObj->m_RefCount > 0);
			--pObj->m_RefCount;
			if (pObj->m_RefCount == 0)
				{
				ObjType Type = pObj->m_Type;

				FreeObj(pObj);
				pObj->OnZeroRefCount();
				}
		}

	**/


#if	TRACE_OBJS
	void LogBusy() const;
	void Up_(Obj *pObj, const char *FileName, unsigned LineNr);
	void Down_(Obj *pObj, const char *FileName, unsigned LineNr);

#define Up(x)	Up_(x, __FILE__, __LINE__)
#define Down(x)	Down_(x, __FILE__, __LINE__)
#define GetSeqInfo()	GetSeqInfo_(__FILE__, __LINE__)
#define GetPathInfo()	GetPathInfo_(__FILE__, __LINE__)
#define GetAlignResult()	GetAlignResult_(__FILE__, __LINE__)
#endif

#if	DEBUG
	void Validate() const;
#endif

private:
	Obj *AllocNew(ObjType Type);
	void FreeObj(Obj *obj);
	unsigned GetFreeCount(ObjType Type) const;
	unsigned GetBusyCount(ObjType Type) const;
	unsigned GetMaxRefCount(ObjType Type) const;
	float GetTotalMem(ObjType Type) const;
#if	DEBUG
	void ValidateType(ObjType Type) const;
#endif
	};

#endif // objmgr_h
