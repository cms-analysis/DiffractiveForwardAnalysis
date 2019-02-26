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

namespace ggll
{
  class HLTMatcher {
    public:
      explicit HLTMatcher() {}
      explicit HLTMatcher( const std::vector<std::string>& );
      ~HLTMatcher() {}

      int TriggerNum( const std::string& );
    private:
      std::vector<std::string> HLTnames;
  };
}

#endif

