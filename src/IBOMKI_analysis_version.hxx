#ifndef IBOMKI_analysis_version_hxx_seen
#define IBOMKI_analysis_version_hxx_seen

// Source for IBOMKI_analysis_version.hxx auto-generated using the
// Ipackage_version.hxx.in template file.

#include <IPackageVersion.hxx>

namespace COMET {
    class IBOMKI_analysis_version;
};

/// This has fields for the library version and can add stuff to the TROOT
/// environment.
class COMET::IBOMKI_analysis_version: public COMET::IPackageVersion {
private:
    static IBOMKI_analysis_version* fThis;
    
public:
    IBOMKI_analysis_version();
    ~IBOMKI_analysis_version();

    /// Return a reference to the singleton.
    static IBOMKI_analysis_version* Get(void); 

    /// Return the version of this library.
    const char* GetName(void) const;

    /// Return the version of this library.
    const char* GetVersion(void) const;

    /// Return the date that this library was compiled.
    const char* GetCompilationDate(void) const;
    
    /// Return the host that this library was compiled on.
    const char* GetCompilationHost(void) const;

    /// Return the directory from which this library was compiled.
    const char* GetCompilationDirectory(void) const;

    /// Return the machine information for the machine that compiled this 
    /// library.  On most machines this is generated from "uname -a".
    const char* GetCompilationMachineInfo(void) const; 

    /// Do any initialization needed for the oaUtility library.  This is
    /// called by the IBOMKI_analysis_version constructor.
    void Initialize(void);
    
    ClassDef(IBOMKI_analysis_version,0);
};
#endif
