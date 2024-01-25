#ifndef objtype_h
#define objtype_h

enum ObjType
	{
#define T(x)	OT_##x,
#include "objtypes.h"
	OTCount
	};

#endif // objtype_h


/// enum ObjType
/// {
///		OT_SeqInfo = 0
///		OT_PathInfo = 1
///		OTCount = 2
/// };