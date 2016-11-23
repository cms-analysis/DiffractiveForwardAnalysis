#ifndef DiffractiveForwardAnalysis_HLTMatcher_h
#define DiffractiveForwardAnalysis_HLTMatcher_h

// system include files
#include <fstream>
#include <memory>
#include <vector>

#include <iostream>

//
// class declaration
//

class HLTMatcher {
  public:
    explicit HLTMatcher() {;}
    explicit HLTMatcher(const std::vector<std::string>&);
    ~HLTMatcher();
    int TriggerNum(const std::string&);
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
