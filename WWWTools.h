#ifndef WWWTools_h
#define WWWTools_h

#include "WWWTree.h"

// Tight
bool isTightSSLepton( int );
bool isTightSSMuon( int );
bool isTightSSElectron( int );

// Loose
bool isLooseSSLepton( int );
bool isLooseSSMuon( int );
bool isLooseSSElectron( int );

// Loose But Not Tight (Lbnt)
bool isLbntSSLepton( int );
bool isLbntSSMuon( int );
bool isLbntSSElectron( int );

// Veto
bool isVetoLepton( int );
bool isVetoMuon( int );
bool isVetoElectron( int );

// Good SS category jet
bool isGoodSSJet( int );

// B-tagging
bool isLooseBJet( int );

#endif
