#ifndef WINSOCKWRAPPER_H
#define WINSOCKWRAPPER_H

#ifdef WIN32

//do not use the embedded max in vs2015 because of conflicts with  std::numeric_limits<int>::max();
//http://stackoverflow.com/questions/1904635/warning-c4003-and-errors-c2589-and-c2059-on-x-stdnumeric-limitsintmax
#define NOMINMAX


//http://www.zachburlingame.com/2011/05/resolving-redefinition-errors-betwen-ws2def-h-and-winsock-h/
#if _MSC_VER > 1000
    #pragma once
#endif

#ifndef _WINDOWS_
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #undef WIN32_LEAN_AND_MEAN
#endif

#include <winsock2.h>

#pragma comment(lib, "ws2_32.lib")

#endif

#endif // WINSOCKWRAPPER_H
