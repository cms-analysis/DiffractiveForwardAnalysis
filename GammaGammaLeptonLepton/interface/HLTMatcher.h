#ifndef DiffractiveForwardAnalysis_HLTMatcher_h
#define DiffractiveForwardAnalysis_HLTMatcher_h

// system include files
#include <fstream>
#include <memory>
#include <vector>

//
// class declaration
//

class HLTMatcher {
  public:
    explicit HLTMatcher(std::vector<std::string>);
    ~HLTMatcher();
    int TriggerNum(std::string);
  private:
    std::vector<std::string> HLTnames;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

#endif
