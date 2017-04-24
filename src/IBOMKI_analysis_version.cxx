#include "IOADatabase.hxx"

// Source for IBOMKI_analysis_version.cxx auto-generated using the
// Ipackage_version.cxx.in template file.

#include "IBOMKI_analysis_version.hxx"
#include "BOMKI_analysis_version.h"

ClassImp(COMET::IBOMKI_analysis_version);

// Trickiness so that the package version is automatically added to the
// list of used packages.
static COMET::IBOMKI_analysis_version BOMKI_analysis_version;

COMET::IBOMKI_analysis_version* COMET::IBOMKI_analysis_version::fThis = NULL;

COMET::IBOMKI_analysis_version::IBOMKI_analysis_version() {
    fThis = COMET::IBOMKI_analysis_version::Get();
}

COMET::IBOMKI_analysis_version::~IBOMKI_analysis_version() {}

void COMET::IBOMKI_analysis_version::Initialize(void) {
    // register this package.
    COMET::IOADatabase::Get().PackageSet().insert(fThis);
}

COMET::IBOMKI_analysis_version* COMET::IBOMKI_analysis_version::Get(void) {
    // Make sure that fThis is initialized;
    if (!fThis) {
        // Make sure that fThis is not null before allocating a real pointer.
        // This cruft is required so that there isn't an infinite recursion
        // while fThis is initialized.
        fThis = (COMET::IBOMKI_analysis_version*) 1;
        // Allocate real space for the fThis pointer.
        fThis = new COMET::IBOMKI_analysis_version;
        // Now initialize
        fThis->Initialize();
    }
    // Return the pointer.
    return fThis;
}

const char* COMET::IBOMKI_analysis_version::GetName(void) const {
    return BOMKI_analysis_NAME;
}

const char* COMET::IBOMKI_analysis_version::GetVersion(void) const {
    return BOMKI_analysis_VERSION;
}

const char* COMET::IBOMKI_analysis_version::GetCompilationDate(void) const {
    return BOMKI_analysis_COMPILE_DATE;
}

const char* COMET::IBOMKI_analysis_version::GetCompilationHost(void) const {
    return BOMKI_analysis_COMPILE_HOST;
}

const char* COMET::IBOMKI_analysis_version::GetCompilationDirectory(void) const {
    return BOMKI_analysis_COMPILE_DIR;
}

const char* COMET::IBOMKI_analysis_version::GetCompilationMachineInfo(void) const {
    return BOMKI_analysis_COMPILE_UNAME;
}
