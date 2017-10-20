#ifndef WWWTools_h
#define WWWTools_h

// WWW_CORE
#include "WWWTree.h"

// ROOT
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

// C++
#include <stdlib.h>
#include <map>
#include <vector>

typedef std::map<TString, std::vector<Int_t>> ObjIdx;
// comparators
bool comparator_abspdgId( int i, int j ) { return abs( wwwbaby.lep_pdgId()[i] ) > abs( wwwbaby.lep_pdgId()[j] ); }
bool comparator_pdgId   ( int i, int j ) { return      wwwbaby.lep_pdgId()[i]   >      wwwbaby.lep_pdgId()[j]  ; }
bool comparator_pt      ( int i, int j ) { return      wwwbaby.lep_p4()[i].pt() >      wwwbaby.lep_p4()[j].pt(); }

extern unsigned int objidx_set_to_eventid;
extern int objidx_set_to_run;
extern int objidx_set_to_lumi;
extern ObjIdx lepidx;
extern ObjIdx jetidx;

extern unsigned int genobjidx_set_to_eventid;
extern int genobjidx_set_to_run;
extern int genobjidx_set_to_lumi;
extern ObjIdx genidx;

extern unsigned int sample_category_set_to_eventid;
extern int sample_category_set_to_run;
extern int sample_category_set_to_lumi;
extern TString sample_category;
extern int sample_priority;

extern unsigned int bkg_category_set_to_eventid;
extern int bkg_category_set_to_run;
extern int bkg_category_set_to_lumi;
extern TString bkg_category;

static int DEBUG;

bool passPresel();

// SS Signal region definitions
bool passSSMM( TString="TightLepton", bool=false, bool=false, bool=false );
bool passSSEM( TString="TightLepton", bool=false, bool=false, bool=false );
bool passSSEE( TString="TightLepton", bool=false, bool=false, bool=false );

// SS Signal region definitions
bool passSSMMPred();
bool passSSEMPred();
bool passSSEEPred();

// Application region definitions
bool passSSARMM();
bool passSSAREM();
bool passSSAREE();

// Application region definitions
bool passSSARMMPred();
bool passSSAREMPred();
bool passSSAREEPred();

// SS B-tag validation region
bool passBTagVRSSMM();
bool passBTagVRSSEM();
bool passBTagVRSSEE();

// SS B-tag validation region
bool passBTagVRSSMMPred();
bool passBTagVRSSEMPred();
bool passBTagVRSSEEPred();

// SS B-tag validation region
bool passBTagARSSMM();
bool passBTagARSSEM();
bool passBTagARSSEE();

// SS B-tag validation region
bool passBTagARSSMMPred();
bool passBTagARSSEMPred();
bool passBTagARSSEEPred();

// SS B-tag validation region
bool passMjjSBVRSSMM();
bool passMjjSBVRSSEM();
bool passMjjSBVRSSEE();

// SS B-tag validation region
bool passMjjSBVRSSMMPred();
bool passMjjSBVRSSEMPred();
bool passMjjSBVRSSEEPred();

// SS B-tag validation region
bool passMjjSBARSSMM();
bool passMjjSBARSSEM();
bool passMjjSBARSSEE();

// SS B-tag validation region
bool passMjjSBARSSMMPred();
bool passMjjSBARSSEMPred();
bool passMjjSBARSSEEPred();

// SS B-tag validation region
bool passMjjSBPRVRSSMM();
bool passMjjSBPRVRSSEM();
bool passMjjSBPRVRSSEE();

// SS B-tag validation region
bool passMjjSBPRVRSSMMPred();
bool passMjjSBPRVRSSEMPred();
bool passMjjSBPRVRSSEEPred();

// SS B-tag validation region
bool passMjjSBPRARSSMM();
bool passMjjSBPRARSSEM();
bool passMjjSBPRARSSEE();

// SS B-tag validation region
bool passMjjSBPRARSSMMPred();
bool passMjjSBPRARSSEMPred();
bool passMjjSBPRARSSEEPred();

// 3L Signal region definitions
bool pass3L0SFOS( TString="Tight3lLepton" );
bool pass3L1SFOS( TString="Tight3lLepton", bool=false, bool=false );
bool pass3L2SFOS( TString="Tight3lLepton", bool=false, bool=false );

// Application region definitions
bool pass3LAR0SFOS();
bool pass3LAR1SFOS();
bool pass3LAR2SFOS();

// Trigger selections
bool passTrigMM();
bool passTrigEM();
bool passTrigEE();

// SS Signal region definitions
bool passPRVRSSMM();
bool passPRVRSSEM();
bool passPRVRSSEE();

// SS Signal region definitions
bool passPRVRSSMMPred();
bool passPRVRSSEMPred();
bool passPRVRSSEEPred();

// Application region definitions
bool passPRARSSMM();
bool passPRARSSEM();
bool passPRARSSEE();

// Application region definitions
bool passPRARSSMMPred();
bool passPRARSSEMPred();
bool passPRARSSEEPred();

// WZ 3 lepton CR
bool passWZCR();

// Event weights
float weight( bool=false, int=0 );

// Print event information
void printEvent();

// Check whether this event's object selections have been performed
bool isObjectSelected();

// Set object indices for this event
void setObjectIndices( bool preload=false );

// Returns a collection of indices for selected leptons
ObjIdx makeEmptyLepidx();
ObjIdx makeEmptyJetidx();

// Returns a collection of indices for selected leptons
ObjIdx getLeptonsIndices();
ObjIdx getJetsIndices();

// Loads collection of indices for selected objects
void loadLeptonIndices();
void loadJetsIndices();

// Loose
bool isVetoLepton( int );
bool isVetoMuon( int );
bool isVetoElec( int );

// Loose
bool isLooseLepton( int );
bool isLooseMuon( int );
bool isLooseElec( int );

// Loose But Not Tight (Lbnt)
bool isLbntLepton( int );
bool isLbntMuon( int );
bool isLbntElec( int );

// Tight
bool isTightLepton( int );
bool isTightMuon( int );
bool isTightElec( int );

// Tight 3L electrons
bool isTight3lLepton( int );
bool isTight3lElec( int );

// Loose
bool isLooseLepton( int );
bool isLooseMuon( int );
bool isLooseElec( int );

// Loose But Not Tight (Lbnt)
bool isLbnt3lLepton( int );
bool isLbnt3lMuon( int );
bool isLbnt3lElec( int );

// Good SS category jet
bool isGoodSSJet( int );

// Good 3L category jet
bool isGood3LJet( int );

// B-tagging
bool isLooseBJet( int );

// Check whether this event's object selections have been performed
bool isGenObjectSelected();

// Set object indices for this event
void setGenObjectIndices();

// Parse the baby version based on the input file path
int getBabyVersion( TString );

// Parse the path and flag if it's www sample
bool isSignal( TString );

// Parse the path and flag if it's data sample
bool isData( TString );

// Mjj calculation with two indices for jets_p4 in WWW baby
float Mjj( int, int );

// Mjj for the tagged W-boson
float MjjW();

// Mjj for the leading two jets
float MjjLead();

// DEtajj for the leading two jets
float DEtajjLead();

// MET
float MET();

// Mll for the leading two leptons in SS channel (return -999 if not two leptons)
float Mll( TString="SignalLepton" );

// Delta Phi for the leading two leptons in SS channel (return -999 if not two leptons)
float DPhill();

// Delta Eta for the leading two leptons in SS channel (return -999 if not two leptons)
float DEtall();

// Ptlll (Pt3l)
float Pt3l();

// Mlll
float M3l();

// DPhilll,MET
float DPhi3lMET();

// MT max of the two SS leptons
float MTmax( TString="SignalLepton" );

// MT with the leading lepton
float MT0( TString="SignalLepton" );

// MT with the sub-leading lepton
float MT1( TString="SignalLepton" );

// MT with the sub-sub-leading lepton
float MT2( TString="SignalLepton" );

// MT with the lepton in given index
float MTidx( int, TString="SignalLepton" );

// Get Lep Flavor
int LepFlavProduct();

// Lepton sign
bool isLeadPlus();
bool isSubleadPlus();
bool isLeadMinus();
bool isSubleadMinus();
bool isPlusPlus();
bool isMinusMinus();

// SS events flavor categorization
bool isOS(); // is two lepton opposite sign event?
bool isSS(); // is two lepton opposite sign event?
bool isSSEE(); // is ee event
bool isSSEM(); // is em event
bool isSSMM(); // is mm event

// Total charge
int totalCharge();

// Number of SFOS pairs
int getNumSFOS();

// Get the Mll of same-flavor pair
float get0SFOSMll();

// Get the Mll of same-electron pair if it doesn't exist, it returns -999
float get0SFOSMee();

// Get the Mll of SFOS pair in 1SFOS
float get1SFOSMll();

// Get the Mll of first SFOS pair in 2SFOS
// (Check the pair in order of (0,1), (0,2) (1,2)
//  and returns the first pair that works)
float get2SFOSMll0();

// Get the Mll of second SFOS pair in 2SFOS
// (Check the pair in order of (1,2), (0,2) (0,1)
//  and returns the first pair that works)
float get2SFOSMll1();

// Four fermion mass
float M4();

// If the sample is VHNonbb, check whether this is WHWWW or non-WHWWW
bool passWHWWW();

// Check the cache
bool isSampleCategorySet();

// Current events background categorization
TString sampleCategory( int& priority=DEBUG );

// Check the cache
bool isBkgCategorySet();

// Current events background categorization (based on what kind of bkg, e.g. trueSS, lost-lepton (LL), or fakes etc.)
TString bkgCategory();

// Get the category of sample based on input filename
TString categ( TString );

// is a a MC sample to be replaced by fake estimation
// If so, we only accept two real from W/Z leptons
bool isFakeSample();

// is a a MC sample to be replaced by fake estimation
// If so, we only accept two real from W/Z leptons
bool isPhotonConversionSample();

// True if the leptons are from all W/Z leptons
bool isRealLeptonEvent();

// True if a lepton comes from photon
bool isPhotonConversionEvent();

// Print event ID in order of evt, run, lumi
void printEventID();

// Fake "Factor"
double fakerate( int, int=0 );
double fakefactor( int, int=0 );


#endif
