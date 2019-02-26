// -*- C++ -*-
//
// Package:    HLTMatcher
// Class:      HLTMatcher
// 
/**\class HLTMatcher HLTMatcher.cc DiffractiveForwardAnalysis/GammaGammaLeptonLepton/src/HLTMatcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme,40 4-B20,+41227671567,
//         Created:  Thu Sep 13 15:17:14 CET 2012
// $Id: HLTMatcher.cc,v 1.3 2013/04/28 08:40:45 lforthom Exp $
//
//

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/HLTMatcher.h"

namespace ggll
{
  HLTMatcher::HLTMatcher(const std::vector<std::string>& list)
  {
    for (unsigned int i=0; i<list.size(); i++) {
      HLTnames.push_back(list[i].substr(0, list[i].find_first_of("*")));
    }
#ifdef DEBUG
    for (unsigned int i=0; i<HLTnames.size(); i++) {
      std::cout << i << " ==> " << HLTnames[i] << std::endl;
    }
#endif
  }

  int
  HLTMatcher::TriggerNum(const std::string& name)
  {
    for (unsigned int i=0; i<HLTnames.size(); i++) {
      if (name.find(HLTnames[i])!=std::string::npos) {
        //std::cout << "--> trigger " << name << " matched with " << HLTnames[i] << std::endl;
        return i;
      }
    }
    return -1;
  }
}
