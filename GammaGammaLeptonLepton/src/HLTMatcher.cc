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
//
// constructors and destructor
//
HLTMatcher::HLTMatcher(std::vector<std::string> _HLTlist)
{
  for (unsigned int i=0; i<_HLTlist.size(); i++) {
    HLTnames.push_back(_HLTlist[i].substr(0, _HLTlist[i].find_first_of("*")));
  }
#ifdef DEBUG
  for (unsigned int i=0; i<HLTnames.size(); i++) {
    std::cout << i << " ==> " << HLTnames[i] << std::endl;
  }
#endif
}

HLTMatcher::~HLTMatcher()
{
}

int
HLTMatcher::TriggerNum(std::string _trigName)
{
  for (unsigned int i=0; i<HLTnames.size(); i++) {
    if (_trigName.find(HLTnames[i])!=std::string::npos) return i;
  }
  return -1;
}
