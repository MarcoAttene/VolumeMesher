#include <stdio.h>

#ifdef _MSC_VER
#include <windows.h>
#include <psapi.h>
// To ensure correct resolution of symbols, add Psapi.lib to TARGETLIBS
// and compile with -DPSAPI_VERSION=1

double getPeakMegabytesUsed()
{
    HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
    if (NULL == hProcess) return 0;

    PROCESS_MEMORY_COUNTERS pmc;
    double mem = 0;
    if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
    {
        mem = pmc.PeakWorkingSetSize / 1048576.0;
    }

    CloseHandle(hProcess);
    return mem;
}

void PrintMemoryInfo()
{
    printf("Peak memory usage (Mb): %f\n", getPeakMegabytesUsed());
}

#else
void PrintMemoryInfo() {
    printf("Peak memory usage unavailable. This feature is implemented on Windows/MSVC only.\n");
}
#endif
