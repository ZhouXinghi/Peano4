// Default variant of CompilerSpecificSettings.h for GCC on Linux.
// For a detailed description of the semantics of the settings, please
// consult LinuxIntel.h which is the most elaborate documentation.

//#define UseTestSpecificCompilerSettings
//#define CompilerCLANG
//#define CompilerICC
#define CompilerHasProcStat
#define CompilerHasUTSName
#define CompilerHasTimespec
#define CompilerHasSysinfo
//#define CompilerDefinesMPIMaxNameString
#if !defined(noMPISupportsSingleSidedCommunication) and !defined(MPISupportsSingleSidedCommunication)
#define MPISupportsSingleSidedCommunication
#endif

#if !defined(UseManualInlining) &&  !defined(noUseManualInlining)
#endif

#define LittleEndian

#if !defined(AlignmentOnHeap)
#define AlignmentOnHeap 16
#endif

#if PeanoDebug == 0
#undef assert
#define assert(...)
#endif

#define KeywordToAvoidDuplicateSymbolsForInlinedFunctions static //inline

#undef OpenMPTaskGroup

#undef OpenMPDependenciesInTaskWait

#define PragmaPush _Pragma("diag_push")
#define SilenceUnknownAttribute _Pragma("diag_suppress unrecognized_attribute")
#define PragmaPop _Pragma("diag_pop")