#ifndef WWWTools_h
#define WWWTools_h

// WWW_CORE
#include "WWWTree.h"

// ROOT
#include "Math/VectorUtil.h"

// C++
#include <stdlib.h>
#include <map>
#include <vector>

typedef std::map<TString, std::vector<unsigned int>> ObjIdx;
// comparators
bool comparator_abspdgId( int i, int j ) { return abs( wwwbaby.lep_pdgId()[i] ) > abs( wwwbaby.lep_pdgId()[j] ); }
bool comparator_pdgId   ( int i, int j ) { return      wwwbaby.lep_pdgId()[i]   >      wwwbaby.lep_pdgId()[j]  ; }
bool comparator_pt      ( int i, int j ) { return      wwwbaby.lep_p4()[i].pt() >      wwwbaby.lep_p4()[j].pt(); }

extern unsigned int objidx_set_to_eventid;
extern ObjIdx lepidx;
extern ObjIdx jetidx;

// SS Signal region definitions (Nminus1_cut_idx == -1 means do not drop any cuts)
bool passSSpresel( int Nminus1_cut_idx = -1 ); // Just object countings
bool passSScommon( int Nminus1_cut_idx = -1 ); // Presel + Common selections (W-tagging, b-veto, high Mjj rejection)
bool passSSMM( int Nminus1_cut_idx = -1 );
bool passSSEM( int Nminus1_cut_idx = -1 );
bool passSSEE( int Nminus1_cut_idx = -1 );

// SS Signal region definitions (Nminus1_cut_idx == -1 means do not drop any cuts)
bool pass3Lpresel( int Nminus1_cut_idx = -1 ); // Just object countings
bool pass3Lcommon( int Nminus1_cut_idx = -1 ); // Presel + Common selections (W-tagging, b-veto, high Mjj rejection)
bool pass3L0SFOS( int Nminus1_cut_idx = -1 );
bool pass3L1SFOS( int Nminus1_cut_idx = -1 );
bool pass3L2SFOS( int Nminus1_cut_idx = -1 );

// Check whether this event's object selections have been performed
bool isObjectSelected();

// Set object indices for this event
void setObjectIndices();

// Returns a collection of indices for selected leptons
ObjIdx getLeptonsIndices();
ObjIdx getJetsIndices();

// Tight
bool isTightSSLepton( int );
bool isTightSSMuon( int );
bool isTightSSElec( int );

// Loose
bool isLooseSSLepton( int );
bool isLooseSSMuon( int );
bool isLooseSSElec( int );

// Loose But Not Tight (Lbnt)
bool isLbntSSLepton( int );
bool isLbntSSMuon( int );
bool isLbntSSElec( int );

// Tight
bool isTightSSLepton( int );
bool isTightSSMuon( int );
bool isTightSSElec( int );

// Loose
bool isLooseSSLepton( int );
bool isLooseSSMuon( int );
bool isLooseSSElec( int );

// Loose But Not Tight (Lbnt)
bool isLbntSSLepton( int );
bool isLbntSSMuon( int );
bool isLbntSSElec( int );

// Good SS category jet
bool isGoodSSJet( int );

// Good 3L category jet
bool isGood3LJet( int );

// B-tagging
bool isLooseBJet( int );

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

// DEtajj for the leading two jets
float MllSS();

// Ptlll (Pt3l)
float Pt3l();

// Get Lep Flavor
int LepFlavProduct();

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

#endif
