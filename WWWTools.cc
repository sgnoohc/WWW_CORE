#include "WWWTools.h"

unsigned int objidx_set_to_eventid;
int objidx_set_to_run;
int objidx_set_to_lumi;
ObjIdx lepidx;
ObjIdx jetidx;

#define MZ 91.1876
#define LUMI 35.87

//======================================================================================
//
//
//
// Same sign analysis signal region selection functions
//
//
//
//======================================================================================

//______________________________________________________________________________________
// The SS dimuon signal region definition
bool passSSMM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                            ==   2     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()         >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()         >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]         >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                          >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                        ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto()       ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                              ==   0     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                                  >   60. &&
                                 MjjW()                                                  <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                               <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                            <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]         == 169     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                                   + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass() >   40.    ) )) { cutidx = 13; return false; }
    cutidx = 14;
    return true;
}

//______________________________________________________________________________________
// The SS emu signal region definition
bool passSSEM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                            ==   2     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()         >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()         >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]         >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                          >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                        ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto()       ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                              ==   0     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                                  >   60. &&
                                 MjjW()                                                  <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                               <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                            <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]         == 143     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                                   + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass() >   30.    ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 15 ) || ( wwwbaby.met_pt()                                        >   40.    ) )) { cutidx = 14; return false; }
    if (!( ( Nm1idx == 16 ) || ( MTmax( "TightLepton" )                                  >   90.    ) )) { cutidx = 15; return false; }
    cutidx = 16;
    return true;
}

//______________________________________________________________________________________
// The SS dielectron signal region definition
bool passSSEE( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                            ==   2     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()         >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()         >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]         >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                          >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                        ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto()       ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                              ==   0     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                                  >   60. &&
                                 MjjW()                                                  <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                               <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                            <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]         == 121     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                                   + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass() >   40.    ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 15 ) || ( wwwbaby.met_pt()                                        >   40.    ) )) { cutidx = 14; return false; }
    if (!( ( Nm1idx == 16 ) || ( ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                                   + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass() <   80. ||
                                 ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                                   + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass() >  100.    ) )) { cutidx = 15; return false; }
    cutidx = 16;
    return true;
}

//======================================================================================
//
//
//
// Same sign analysis fake rate b-tagging validation region
//
//
//
//======================================================================================

//______________________________________________________________________________________
bool passBTagVRSSMM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                      ==   2     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()   >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()   >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]   >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                    >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                        >=   1     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                            >   60. &&
                                 MjjW()                                            <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                         <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                      <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]   == 169     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( Mll( "TightLepton" )                              >   40.    ) )) { cutidx = 13; return false; }
    cutidx = 14;
    return true;
}

//______________________________________________________________________________________
bool passBTagVRSSEM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                      ==   2     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()   >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()   >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]   >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                    >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                        >=   1     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                            >   60. &&
                                 MjjW()                                            <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                         <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                      <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]   == 143     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( Mll( "TightLepton" )                              >   30.    ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 15 ) || ( wwwbaby.met_pt()                                  >   40.    ) )) { cutidx = 14; return false; }
    if (!( ( Nm1idx == 16 ) || ( MTmax( "TightLepton" )                            >   90.    ) )) { cutidx = 15; return false; }
    cutidx = 16;
    return true;
}

//______________________________________________________________________________________
bool passBTagVRSSEE( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                      ==   2     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()   >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()   >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]   >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                    >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                        >=   1     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                            >   60. &&
                                 MjjW()                                            <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                         <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                      <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]]   == 121     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( Mll( "TightLepton" )                              >   40.    ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 15 ) || ( wwwbaby.met_pt()                                  >   40.    ) )) { cutidx = 14; return false; }
    if (!( ( Nm1idx == 16 ) || ( Mll( "TightLepton" )                              <   80. ||
                                 Mll( "TightLepton" )                              >  100.    ) )) { cutidx = 15; return false; }
    cutidx = 16;
    return true;
}

//======================================================================================
//
//
//
// Same sign analysis fake rate b-tagging application region
//
//
//
//======================================================================================

//______________________________________________________________________________________
bool passBTagARSSMM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  0 ) || ( lepidx["LooseLepton"].size()                      ==   2     ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                      ==   1     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]].pt()   >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]].pt()   >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]]   >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                    >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                        >=   1     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                            >   60. &&
                                 MjjW()                                            <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                         <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                      <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]]   == 169     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( Mll( "LooseLepton" )                              >   40.    ) )) { cutidx = 13; return false; }
    cutidx = 14;
    return true;
}

//______________________________________________________________________________________
bool passBTagARSSEM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  0 ) || ( lepidx["LooseLepton"].size()                      ==   2     ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                      ==   1     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]].pt()   >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]].pt()   >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]]   >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                    >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                        >=   1     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                            >   60. &&
                                 MjjW()                                            <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                         <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                      <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]]   == 143     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( Mll( "LooseLepton" )                              >   30.    ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 15 ) || ( wwwbaby.met_pt()                                  >   40.    ) )) { cutidx = 14; return false; }
    if (!( ( Nm1idx == 16 ) || ( MTmax( "LooseLepton" )                            >   90.    ) )) { cutidx = 15; return false; }
    cutidx = 16;
    return true;
}

//______________________________________________________________________________________
bool passBTagARSSEE( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx ==  0 ) || ( lepidx["LooseLepton"].size()                      ==   2     ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx ==  1 ) || ( lepidx["TightLepton"].size()                      ==   1     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx ==  2 ) || ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]].pt()   >   30.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx ==  3 ) || ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]].pt()   >   30.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx ==  4 ) || ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]]   >    0     ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx ==  5 ) || ( jetidx["GoodSSJet"]    .size()                    >=   2     ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx ==  6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   2     ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx ==  7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0     ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx ==  8 ) || ( jetidx["LooseBJet"].size()                        >=   1     ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx ==  9 ) || ( MjjW()                                            >   60. &&
                                 MjjW()                                            <  100.    ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 11 ) || ( MjjLead()                                         <  400.    ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 12 ) || ( DEtajjLead()                                      <    1.5   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 13 ) || ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]]
                                 * wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]]   == 121     ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 14 ) || ( Mll( "LooseLepton" )                              >   40.    ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 15 ) || ( wwwbaby.met_pt()                                  >   40.    ) )) { cutidx = 14; return false; }
    if (!( ( Nm1idx == 16 ) || ( Mll( "LooseLepton" )                              <   80. ||
                                 Mll( "LooseLepton" )                              >  100.    ) )) { cutidx = 15; return false; }
    cutidx = 16;
    return true;
}

//======================================================================================
//
//
//
// Three lepton analysis signal region selection functions
//
//
//
//======================================================================================

//______________________________________________________________________________________
// The 3L0SFOS signal region definition
bool pass3L0SFOS( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1  ) || ( lepidx["3LTightLepton"].size()                    ==   3   ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx == 2  ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   3   ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx == 3  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][0]].pt() >   25.  ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx == 4  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][1]].pt() >   20.  ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx == 5  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][2]].pt() >   20.  ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx == 6  ) || ( abs( totalCharge() )                              ==   1   ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx == 7  ) || ( jetidx["Good3LJet"].size()                        <=   1   ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx == 8  ) || ( pass3Lpresel()                                             ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx == 9  ) || ( Pt3l()                                            >   60.  ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 10 ) || ( jetidx["LooseBJet"].size()                        ==   0   ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 11 ) || ( getNumSFOS()                                      ==   0   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 12 ) || ( get0SFOSMll()                                     >   20.  ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 13 ) || ( fabs( get0SFOSMee() - MZ )                        >   15.  ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 14 ) || ( DPhi3lMET()                                       >    2.7 ) )) { cutidx = 14; return false; }
    cutidx = 14;
    return true;
}

//______________________________________________________________________________________
// The 3L1SFOS signal region definition
bool pass3L1SFOS( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1  ) || ( lepidx["3LTightLepton"].size()                    ==   3   ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx == 2  ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   3   ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx == 3  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][0]].pt() >   25.  ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx == 4  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][1]].pt() >   20.  ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx == 5  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][2]].pt() >   20.  ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx == 6  ) || ( abs( totalCharge() )                              ==   1   ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx == 7  ) || ( jetidx["Good3LJet"].size()                        <=   1   ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx == 8  ) || ( pass3Lpresel()                                             ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx == 9  ) || ( Pt3l()                                            >   60.  ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 10 ) || ( jetidx["LooseBJet"].size()                        ==   0   ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 11 ) || ( getNumSFOS()                                      ==   1   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 12 ) || ( wwwbaby.met_pt()                                  >   45.  ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 13 ) || ( get1SFOSMll()                                     <   55.||
                                 get1SFOSMll()                                     >  100.  ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 14 ) || ( DPhi3lMET()                                       >    2.5 ) )) { cutidx = 14; return false; }
    cutidx = 14;
    return true;
}

//______________________________________________________________________________________
// The 3L2SFOS signal region definition
bool pass3L2SFOS( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1  ) || ( lepidx["3LTightLepton"].size()                    ==   3   ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx == 2  ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==   3   ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx == 3  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][0]].pt() >   25.  ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx == 4  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][1]].pt() >   20.  ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx == 5  ) || ( wwwbaby.lep_p4()[lepidx["3LTightLepton"][2]].pt() >   20.  ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx == 6  ) || ( abs( totalCharge() )                              ==   1   ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx == 7  ) || ( jetidx["Good3LJet"].size()                        <=   1   ) )) { cutidx = 7 ; return false; }
    if (!( ( Nm1idx == 8  ) || ( pass3Lpresel()                                             ) )) { cutidx = 8 ; return false; }
    if (!( ( Nm1idx == 9  ) || ( Pt3l()                                            >   60.  ) )) { cutidx = 9 ; return false; }
    if (!( ( Nm1idx == 10 ) || ( jetidx["LooseBJet"].size()                        ==   0   ) )) { cutidx = 10; return false; }
    if (!( ( Nm1idx == 11 ) || ( getNumSFOS()                                      ==   2   ) )) { cutidx = 11; return false; }
    if (!( ( Nm1idx == 12 ) || ( wwwbaby.met_pt()                                  >   55.  ) )) { cutidx = 12; return false; }
    if (!( ( Nm1idx == 13 ) || ( fabs( get2SFOSMll0() - MZ )                        >   20. ) )) { cutidx = 13; return false; }
    if (!( ( Nm1idx == 14 ) || ( fabs( get2SFOSMll1() - MZ )                        >   20. ) )) { cutidx = 14; return false; }
    if (!( ( Nm1idx == 15 ) || ( DPhi3lMET()                                       >    2.5 ) )) { cutidx = 15; return false; }
    cutidx = 15;
    return true;
}

//======================================================================================
// Some common selections that are not mainly used
//======================================================================================

//______________________________________________________________________________________
bool passSSpresel( int Nm1idx, int& cutidx )
{
    setObjectIndices();
//    if (!( ( Nm1idx == 0 ) || ( wwwbaby.HLT_DoubleMu_noiso() || wwwbaby.HLT_DoubleMu()   ) )) { cutidx = 0; return false; }
    if (!( ( Nm1idx == 1 ) || ( lepidx["SignalLepton"].size()                     ==  2  ) )) { cutidx = 1; return false; }
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]].pt()  >  30. ) )) { cutidx = 2; return false; }
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.lep_p4()[lepidx["SignalLepton"][1]].pt()  >  30. ) )) { cutidx = 3; return false; }
    if (!( ( Nm1idx == 4 ) || ( isSS()                                                   ) )) { cutidx = 4; return false; }
    if (!( ( Nm1idx == 5 ) || ( jetidx["GoodSSJet"]    .size()                    >=  2  ) )) { cutidx = 5; return false; }
    if (!( ( Nm1idx == 6 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==  2  ) )) { cutidx = 6; return false; }
    if (!( ( Nm1idx == 7 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==  0  ) )) { cutidx = 7; return false; }
    cutidx = 8;
    return true;
}

//______________________________________________________________________________________
bool passSScommon( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSSpresel()                       ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( jetidx["LooseBJet"].size() == 0      ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( MjjW()                     > 60.  &&
                                MjjW()                     < 100.    ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( MjjLead()                  < 400.    ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( DEtajjLead()               < 1.5     ) )) { cutidx = 4 ; return false; }
    cutidx = 5;
    return true;
}

//______________________________________________________________________________________
bool pass3Lpresel( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1 ) || ( lepidx["SignalLepton"].size()                     ==  3  ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==  3  ) )) { cutidx = 2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]].pt()  >  25. ) )) { cutidx = 3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( wwwbaby.lep_p4()[lepidx["SignalLepton"][1]].pt()  >  20. ) )) { cutidx = 4 ; return false; }
    if (!( ( Nm1idx == 5 ) || ( wwwbaby.lep_p4()[lepidx["SignalLepton"][2]].pt()  >  20. ) )) { cutidx = 5 ; return false; }
    if (!( ( Nm1idx == 6 ) || ( abs( totalCharge() )                              ==  1  ) )) { cutidx = 6 ; return false; }
    if (!( ( Nm1idx == 7 ) || ( jetidx["Good3LJet"].size()                        <=  1  ) )) { cutidx = 7 ; return false; }
    cutidx = 8;
    return true;
}

//______________________________________________________________________________________
bool pass3Lcommon( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( pass3Lpresel()                   ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( Pt3l() > 60.                     ) )) { cutidx = 1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( jetidx["LooseBJet"].size() ==  0 ) )) { cutidx = 2 ; return false; }
    cutidx = 3;
    return true;
}

//______________________________________________________________________________________
bool passSSARMM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    setSignalLeptonToOneTightAndOneLbnt();
    bool pass = passSSMM( Nm1idx, cutidx );
    setSignalLeptonToTightLepton();
    return pass;
}

//______________________________________________________________________________________
bool passSSAREM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    setSignalLeptonToOneTightAndOneLbnt();
    bool pass = passSSEM( Nm1idx, cutidx );
    setSignalLeptonToTightLepton();
    return pass;
}

//______________________________________________________________________________________
bool passSSAREE( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    setSignalLeptonToOneTightAndOneLbnt();
    bool pass = passSSEE( Nm1idx, cutidx );
    setSignalLeptonToTightLepton();
    return pass;
}

//______________________________________________________________________________________
float weight()
{
    float wgt = wwwbaby.isData() ? 1 : wwwbaby.evt_scale1fb() * LUMI;

    if ( wwwbaby.evt_dataset()[0].Contains( "/TEST-www/wwwext-Private80X-v1/USER" ) )
        wgt *= 0.066805*164800./(91900.+164800.);
    if ( wwwbaby.evt_dataset()[0].Contains( "/TEST-www/www-Private80X-v1/USER" ) )
        wgt *= 0.066805*91900./(91900.+164800.);

    return wgt;
}

//______________________________________________________________________________________
void printEvent()
{
    std::cout << "=======================================================" << std::endl;
    std::cout << " run  = " << wwwbaby.run() << std::endl;
    std::cout << " lumi = " << wwwbaby.lumi() << std::endl;
    std::cout << " evt  = " << wwwbaby.evt() << std::endl;
    for ( auto& lidx : lepidx )
    {
        std::cout << "Printing ... " << lidx.first << std::endl;
        for ( auto& ilep : lidx.second )
        {
            std::cout << wwwbaby.lep_p4()[ilep].pt()     << " ";
            std::cout << wwwbaby.lep_p4()[ilep].eta()    << " ";
            std::cout << wwwbaby.lep_p4()[ilep].phi()    << " ";
            std::cout << wwwbaby.lep_p4()[ilep].energy() << " ";
            std::cout << wwwbaby.lep_p4()[ilep].mass()   << " ";
            std::cout << wwwbaby.lep_pdgId()[ilep]       << " ";
            std::cout << std::endl;
        }
    }
    for ( auto& jidx : jetidx )
    {
        std::cout << "Printing ... " << jidx.first << std::endl;
        for ( auto& ijet : jidx.second )
        {
            std::cout << wwwbaby.jets_p4()[ijet].pt()     << " ";
            std::cout << wwwbaby.jets_p4()[ijet].eta()    << " ";
            std::cout << wwwbaby.jets_p4()[ijet].phi()    << " ";
            std::cout << wwwbaby.jets_p4()[ijet].energy() << " ";
            std::cout << wwwbaby.jets_p4()[ijet].mass()   << " ";
            std::cout << wwwbaby.jets_csv()[ijet]         << " ";
            std::cout << std::endl;
        }
    }
    int wwwdebug = 0;
    std::cout  << "passSSEM      " << " : " << passSSEM       (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "passSSpresel  " << " : " << passSSpresel   (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "passSScommon  " << " : " << passSScommon   (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "passSSMM      " << " : " << passSSMM       (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "passSSEM      " << " : " << passSSEM       (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "passSSEE      " << " : " << passSSEE       (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "pass3Lpresel  " << " : " << pass3Lpresel   (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "pass3Lcommon  " << " : " << pass3Lcommon   (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "pass3L0SFOS   " << " : " << pass3L0SFOS    (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "pass3L1SFOS   " << " : " << pass3L1SFOS    (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "pass3L2SFOS   " << " : " << pass3L2SFOS    (-1, wwwdebug) << " ";
    std::cout  << wwwdebug << std::endl;
    std::cout  << "MjjW          " << " : " << MjjW           () << std::endl;
    std::cout  << "MjjLead       " << " : " << MjjLead        () << std::endl;
    std::cout  << "DEtajjLead    " << " : " << DEtajjLead     () << std::endl;
    std::cout  << "Mll           " << " : " << Mll            () << std::endl;
    std::cout  << "Pt3l          " << " : " << Pt3l           () << std::endl;
    std::cout  << "DPhi3lMET     " << " : " << DPhi3lMET      () << std::endl;
    std::cout  << "LepFlavProduct" << " : " << LepFlavProduct () << std::endl;
    std::cout  << "isOS          " << " : " << isOS           () << std::endl;
    std::cout  << "isSS          " << " : " << isSS           () << std::endl;
    std::cout  << "isSSEE        " << " : " << isSSEE         () << std::endl;
    std::cout  << "isSSEM        " << " : " << isSSEM         () << std::endl;
    std::cout  << "isSSMM        " << " : " << isSSMM         () << std::endl;
    std::cout  << "totalCharge   " << " : " << totalCharge    () << std::endl;
    std::cout  << "getNumSFOS    " << " : " << getNumSFOS     () << std::endl;
    std::cout  << "get0SFOSMll   " << " : " << get0SFOSMll    () << std::endl;
    std::cout  << "get0SFOSMee   " << " : " << get0SFOSMee    () << std::endl;
    std::cout  << "get1SFOSMll   " << " : " << get1SFOSMll    () << std::endl;
    std::cout  << "get2SFOSMll0  " << " : " << get2SFOSMll0   () << std::endl;
    std::cout  << "get2SFOSMll1  " << " : " << get2SFOSMll1   () << std::endl;
    std::cout << "=======================================================" << std::endl;
}

//______________________________________________________________________________________
void setObjectIndices()
{
    // If the object already selected skip
    if ( isObjectSelected() ) return;
    //otherwise perform object selections
    lepidx = getLeptonsIndices();
    jetidx = getJetsIndices();
    // Now set the event id of the event the objects are selected from
    objidx_set_to_eventid = wwwbaby.evt();
    objidx_set_to_run     = wwwbaby.run();
    objidx_set_to_lumi    = wwwbaby.lumi();
}

//______________________________________________________________________________________
void setSignalLeptonToOneTightAndOneLbnt()
{
    lepidx["SignalLepton"].clear();
    if ( lepidx["TightLepton"].size() == 1 &&
         lepidx["LbntLepton"] .size() == 1 )
    {
        lepidx["SignalLepton"].push_back( lepidx["TightLepton"][0] );
        lepidx["SignalLepton"].push_back( lepidx["LbntLepton"] [0] );
    }
}

//______________________________________________________________________________________
void setSignalLeptonToTightLepton()
{
    lepidx["SignalLepton"] = lepidx["TightLepton"];
    lepidx["SignalMuon"] = lepidx["TightMuon"];
    lepidx["SignalElec"] = lepidx["TightElec"];
}

//______________________________________________________________________________________
bool isObjectSelected()
{
    if ( objidx_set_to_eventid == wwwbaby.evt() &&
         objidx_set_to_run     == wwwbaby.run() &&
         objidx_set_to_lumi    == wwwbaby.lumi() )
        return true;
    else
        return false;
}

//______________________________________________________________________________________
ObjIdx getLeptonsIndices()
{
    ObjIdx idx;
    for ( unsigned int ilep = 0; ilep < wwwbaby.lep_p4().size(); ++ilep )
    {
        if ( is3LTightLepton( ilep ) ) idx["3LTightLepton"].push_back( ilep );
        if ( is3LTightElec  ( ilep ) ) idx["3LTightElec"]  .push_back( ilep );
        if ( isTightLepton  ( ilep ) ) idx["TightLepton"]  .push_back( ilep );
        if ( isTightMuon    ( ilep ) ) idx["TightMuon"]    .push_back( ilep );
        if ( isTightElec    ( ilep ) ) idx["TightElec"]    .push_back( ilep );
        if ( isLooseLepton  ( ilep ) ) idx["LooseLepton"]  .push_back( ilep );
        if ( isLooseMuon    ( ilep ) ) idx["LooseMuon"]    .push_back( ilep );
        if ( isLooseElec    ( ilep ) ) idx["LooseElec"]    .push_back( ilep );
        if ( isLbntLepton   ( ilep ) ) idx["LbntLepton"]   .push_back( ilep );
        if ( isLbntMuon     ( ilep ) ) idx["LbntMuon"]     .push_back( ilep );
        if ( isLbntElec     ( ilep ) ) idx["LbntElec"]     .push_back( ilep );
    }
    // Default is to set signal leptons to Tight leptons
    idx["SignalLepton"] = idx["TightLepton"];
    idx["SignalMuon"] = idx["TightMuon"];
    idx["SignalElec"] = idx["TightElec"];
    return idx;
}

//______________________________________________________________________________________
ObjIdx getJetsIndices()
{
    ObjIdx idx;
    for ( unsigned int ijet = 0; ijet < wwwbaby.jets_p4().size(); ++ijet )
    {
        if ( isGoodSSJet( ijet ) ) idx["GoodSSJet"].push_back( ijet );
        if ( isGood3LJet( ijet ) ) idx["Good3LJet"].push_back( ijet );
        if ( isLooseBJet( ijet ) ) idx["LooseBJet"].push_back( ijet );
    }

    // If there are less than two good jets we skip the W-boson jet selection
    if ( idx["GoodSSJet"].size() < 2 )
        return idx;

    // Special selection select W-boson jets
    float mindr = 999;
    std::pair<int, int> wjetpair;
    for ( unsigned int ijet = 0; ijet < idx["GoodSSJet"].size(); ++ijet )
    {
        for ( unsigned int jjet = ijet + 1; jjet < idx["GoodSSJet"].size(); ++jjet )
        {
            if ( ijet == jjet )
                continue;
            LorentzVector ijet_p4 = wwwbaby.jets_p4()[idx["GoodSSJet"][ijet]];
            LorentzVector jjet_p4 = wwwbaby.jets_p4()[idx["GoodSSJet"][jjet]];
            float dr = ROOT::Math::VectorUtil::DeltaR( ijet_p4, jjet_p4 );
            if ( dr < mindr )
            {
                mindr = dr;
                wjetpair = make_pair( idx["GoodSSJet"][ijet], idx["GoodSSJet"][jjet] );
            }
        }
    }

    // Now set the value
    if ( mindr < 999 )
    {
        idx["GoodSSWJet"].push_back( wjetpair.first );
        idx["GoodSSWJet"].push_back( wjetpair.second );
    }

    return idx;
}

//______________________________________________________________________________________
bool isTightLepton( int ilep )
{
    return isTightMuon( ilep ) || isTightElec( ilep );
}

//______________________________________________________________________________________
bool isTightMuon( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 13     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isTightElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 11     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0     )) return false;
    return true;
}

//______________________________________________________________________________________
bool is3LTightLepton( int ilep )
{
    return isTightMuon( ilep ) || is3LTightElec( ilep );
}

//______________________________________________________________________________________
bool is3LTightElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 11     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isSUSYSSTightMuon( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 13     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseLepton( int ilep )
{
    return isLooseMuon( ilep ) || isLooseElec( ilep );
}

//______________________________________________________________________________________
bool isLooseMuon( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_fo()[ilep]    )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 13     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_fo_noiso()[ilep] )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 11        )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5      )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6      )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0        )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLbntLepton( int ilep )
{
    return isLbntMuon( ilep ) || isLbntElec( ilep );
}

//______________________________________________________________________________________
bool isLbntMuon( int ilep )
{
    return isLooseMuon( ilep ) && !isTightMuon( ilep );
}

//______________________________________________________________________________________
bool isLbntElec( int ilep )
{
    return isLooseElec( ilep ) && !isTightElec( ilep );
}

//______________________________________________________________________________________
bool isGoodSSJet( int ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].pt()      > 30.  )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].eta() )   <  2.5 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isGood3LJet( int ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].pt()      > 30.  )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].eta() )   <  5.0 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseBJet( int ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].pt()      > 20.     )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].eta() )   <  2.4    )) return false;
    if (!(       wwwbaby.jets_csv()[ijet]          >  0.5426 )) return false;
    return true;
}

//______________________________________________________________________________________
int getBabyVersion( TString path )
{
    if ( path.Contains( "v0.1.5" ) ) return 5;
    if ( path.Contains( "v0.1.9" ) ) return 9;
    if ( path.Contains( "v0.1.11" ) ) return 11;
    if ( path.Contains( "v0_1_5" ) ) return 5;
    if ( path.Contains( "v0_1_9" ) ) return 9;
    if ( path.Contains( "v0_1_11" ) ) return 11;
    std::cout << "Error:: Unrecognized version number!" << std::endl;
    exit(-1);
    return -1;
}

//______________________________________________________________________________________
bool isSignal( TString path )
{
    if ( path.Contains( "www_2l" ) ) return true;
    if ( path.Contains( "www_incl" ) ) return true;
    return false;
}

//______________________________________________________________________________________
bool isData( TString path )
{
    if ( path.Contains( "data_Run" ) ) return true;
    return false;
}

//______________________________________________________________________________________
float Mjj( int ijet, int jjet )
{
    return ( wwwbaby.jets_p4()[ijet] + wwwbaby.jets_p4()[jjet] ).mass();
}

//______________________________________________________________________________________
float MjjW()
{
    if ( jetidx["GoodSSWJet"].size() < 2 )
        return -999;
    return Mjj( jetidx["GoodSSWJet"][0], jetidx["GoodSSWJet"][1] );
}

//______________________________________________________________________________________
float MjjLead()
{
    if ( jetidx["GoodSSJet"].size() < 2 )
        return -999;
    return Mjj( jetidx["GoodSSJet"][0], jetidx["GoodSSJet"][1] );
}

//______________________________________________________________________________________
float DEtajjLead()
{
    if ( jetidx["GoodSSJet"].size() < 2 )
        return -999;
    return fabs(
            wwwbaby.jets_p4()[jetidx["GoodSSJet"][0]].eta()
            - wwwbaby.jets_p4()[jetidx["GoodSSJet"][1]].eta() );
}

//______________________________________________________________________________________
float Mll( TString leptype )
{
    if ( lepidx[leptype].size() < 2 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx[leptype][0]]
             + wwwbaby.lep_p4()[lepidx[leptype][1]] ).mass();
}

//______________________________________________________________________________________
float DPhill()
{
    if ( lepidx["SignalLepton"].size() < 2 )
        return -999;
    return fabs( ROOT::Math::VectorUtil::DeltaPhi(
            wwwbaby.lep_p4()[lepidx["SignalLepton"][0]],
            wwwbaby.lep_p4()[lepidx["SignalLepton"][1]] ));
}

//______________________________________________________________________________________
float DEtall()
{
    if ( lepidx["SignalLepton"].size() < 2 )
        return -999;
    return fabs( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]].eta()
             - wwwbaby.lep_p4()[lepidx["SignalLepton"][1]].eta() );
}

//______________________________________________________________________________________
float Pt3l()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
             + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).pt();
}

//______________________________________________________________________________________
float DPhi3lMET()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    TLorentzVector met;
    met.SetPtEtaPhiM( wwwbaby.met_pt(), 0, wwwbaby.met_phi(), 0 );
    LorentzVector metlv;
    metlv.SetPxPyPzE( met.Px(), met.Py(), met.Pz(), met.E() );
    return fabs( ROOT::Math::VectorUtil::DeltaPhi(
                wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
                + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]],
                met ));
}

//______________________________________________________________________________________
float MTmax( TString leptype )
{
    if ( lepidx[leptype].size() < 2 )
        return -999;
    float mt0 = MT0();
    float mt1 = MT1();
    return mt0 > mt1 ? mt0 : mt1;
}

//______________________________________________________________________________________
float MT0( TString leptype )
{
    return MTidx( 0, leptype );
}

//______________________________________________________________________________________
float MT1( TString leptype )
{
    return MTidx( 1, leptype );
}

//______________________________________________________________________________________
float MT2( TString leptype )
{
    return MTidx( 2, leptype );
}

//______________________________________________________________________________________
float MTidx( int idx, TString leptype )
{
    if ( lepidx[leptype].size() < (unsigned int) idx + 1 )
        return -999;
    // Set MET
    TLorentzVector met;
    met.SetPtEtaPhiM( wwwbaby.met_pt(), 0, wwwbaby.met_phi(), 0 );
    LorentzVector metlv;
    metlv.SetPxPyPzE( met.Px(), met.Py(), met.Pz(), met.E() );
    // Get Lepton p4
    LorentzVector p4;
    p4 = wwwbaby.lep_p4()[lepidx[leptype][idx]];
    return sqrt( 2 * p4.pt() * metlv.pt() * ( 1 - cos( ROOT::Math::VectorUtil::DeltaPhi( metlv, p4 ) ) ) );
}

//______________________________________________________________________________________
int LepFlavProduct()
{
    if ( lepidx["SignalLepton"].size() < 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]]
           * wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]];
}

//______________________________________________________________________________________
bool isLeadPlus()
{
    if ( lepidx["SignalLepton"].size() < 1 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]] < 0;
}

//______________________________________________________________________________________
bool isSubleadPlus()
{
    if ( lepidx["SignalLepton"].size() < 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]] < 0;
}

//______________________________________________________________________________________
bool isLeadMinus()
{
    if ( lepidx["SignalLepton"].size() < 1 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]] > 0;
}

//______________________________________________________________________________________
bool isSubleadMinus()
{
    if ( lepidx["SignalLepton"].size() < 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]] > 0;
}

//______________________________________________________________________________________
bool isPlusPlus()
{
    if ( isLeadPlus() && isSubleadPlus() )
        return true;
    else
        return false;
}

//______________________________________________________________________________________
bool isMinusMinus()
{
    if ( isLeadMinus() && isSubleadMinus() )
        return true;
    else
        return false;
}

//______________________________________________________________________________________
bool isOS() { return LepFlavProduct() < 0; }
bool isSS() { return LepFlavProduct() > 0; }
bool isSSEE() { return LepFlavProduct() == 121; }
bool isSSEM() { return LepFlavProduct() == 143; }
bool isSSMM() { return LepFlavProduct() == 169; }

//______________________________________________________________________________________
int totalCharge()
{
    int charge = 0;
    for ( auto& ilep : lepidx["SignalLepton"] )
        charge += wwwbaby.lep_charge()[ilep];
    return charge;
}

//______________________________________________________________________________________
int getNumSFOS()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    /* Loops through pairs of entries in the lep_pdgId vector
       and counts how many have opposite value*/
    int num_SFOS = 0;
    for ( int i = 0; i < ( int ) lepidx["SignalLepton"].size(); i++ )
    {
        for ( int j = i + 1; j < ( int ) lepidx["SignalLepton"].size(); j++ )
        {
            if ( wwwbaby.lep_pdgId()[lepidx["SignalLepton"][i]] ==
                    -( wwwbaby.lep_pdgId()[lepidx["SignalLepton"][j]] ) )
                num_SFOS++;
        }
    }
    return num_SFOS;
}

//______________________________________________________________________________________
float get0SFOSMll()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][2]];
    if ( pdgid0 == pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get0SFOSMee()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][2]];
    if ( pdgid0 == pdgid1 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 && abs( pdgid1 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get1SFOSMll()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get2SFOSMll0()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get2SFOSMll1()
{
    if ( lepidx["SignalLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["SignalLepton"][2]];
    if ( pdgid2 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid0 )
        return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["SignalLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float M4()
{
    if ( lepidx["SignalLepton"].size() != 2 )
        return -999;
    if ( jetidx["GoodSSWJet"].size() != 2 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx["SignalLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["SignalLepton"][1]]
             + wwwbaby.jets_p4()[jetidx["GoodSSWJet"][0]]
             + wwwbaby.jets_p4()[jetidx["GoodSSWJet"][1]]
             ).mass();
}

//______________________________________________________________________________________
bool passWHWWW()
{
    // VH Non-bb sample counting
    int nW = 0;
    int nZ = 0;
    int nWfromH = 0;
    int nZfromH = 0;

    if ( wwwbaby.evt_dataset()[0].Contains( "/VHToNonbb" ) )
    {
        for ( unsigned int igen = 0; igen < wwwbaby.genPart_pdgId().size(); ++igen )
        {
            if ( abs( wwwbaby.genPart_pdgId().at( igen ) )  == 24
                   && wwwbaby.genPart_status().at( igen )   == 22
                   && wwwbaby.genPart_motherId().at( igen ) != 25 )
                nW++;

            if ( abs( wwwbaby.genPart_pdgId().at( igen ) )  == 23
                   && wwwbaby.genPart_status().at( igen )   == 22
                   && wwwbaby.genPart_motherId().at( igen ) != 25 )
                nZ++;

            if ( abs( wwwbaby.genPart_pdgId().at( igen ) )  == 24
                   && wwwbaby.genPart_status().at( igen )   == 22
                   && wwwbaby.genPart_motherId().at( igen ) == 25 )
                nWfromH++;

            if ( abs( wwwbaby.genPart_pdgId().at( igen ) )  == 23
                   && wwwbaby.genPart_status().at( igen )   == 22
                   && wwwbaby.genPart_motherId().at( igen ) == 25 )
                nZfromH++;
        }

        // WH
        if ( nW == 1 )
        {
            if ( nWfromH != 2 )
                return false;
            else
                return true;
        }
        // ZH
        else if ( nZ == 1 )
            return false;
        // unknown
        else
            return false;
    }

    return true;
}

//______________________________________________________________________________________
TString sampleCategory( int& priority )
{
    const TString& dsname = wwwbaby.evt_dataset()[0];
    if ( dsname.Contains( "/QCD_HT" ) )                                            { priority = 2; return "qcd"; }
    if ( dsname.Contains( "/GJets" ) )                                             { priority = 2; return "gj"; }
    if ( dsname.Contains( "/WJetsToLNu" ) )                                        { priority = 1; return "wj"; }
    if ( dsname.Contains( "/DYJetsToLL_M" ) && dsname.Contains( "madgraphMLM" ) )  { priority = 1; return "dy"; }
    if ( dsname.Contains( "/DYJetsToLL_M" ) && dsname.Contains( "amcatnloFXFX" ) ) { priority = 2; return "dynlo"; }
    if ( dsname.Contains( "/TTJets_Single" ) && dsname.Contains( "ext1" ) )        { priority = 1; return "tt1l"; }
    if ( dsname.Contains( "/TTJets_Single" ) && !dsname.Contains( "ext1" ) )       { priority = 2; return "tt1lnonext"; }
    if ( dsname.Contains( "/TTJets_Di" ) && dsname.Contains( "ext1" ) )            { priority = 1; return "tt2l"; }
    if ( dsname.Contains( "/TTJets_Di" ) && !dsname.Contains( "ext1") )            { priority = 2; return "tt2lnonext"; }
    if ( dsname.Contains( "/TTTo2L2Nu_" ) )                                        { priority = 2; return "tt2lpowheg"; }
    if ( dsname.Contains( "TTTo2L2Nu_" ) )                                         { priority = 2; return "tt2lpowheg"; } // some weird ones have weird names
    if ( dsname.Contains( "/ST_" ) )                                               { priority = 1; return "singletop"; }
    if ( dsname.Contains( "/tZq" ) )                                               { priority = 1; return "singletop"; }
    if ( dsname.Contains( "/WWTo" ) && dsname.Contains( "powheg" ) )               { priority = 1; return "ww"; }
    if ( dsname.Contains( "/WWTo" ) && dsname.Contains( "DoubleScattering" ) )     { priority = 1; return "wwdpi"; }
    if ( dsname.Contains( "/WpWpJJ_EWK" ) )                                        { priority = 1; return "vbsww"; }
    if ( dsname.Contains( "/WmWmJJ_13" ) )                                         { priority = 2; return "wwjj"; }
    if ( dsname.Contains( "/WpWpJJ_13" ) )                                         { priority = 2; return "wwjj"; }
    if ( dsname.Contains( "/WZTo" ) )                                              { priority = 1; return "wz"; }
    if ( dsname.Contains( "/GluGluHToZZTo4L" ) )                                   { priority = 2; return "gghzz"; }
    if ( dsname.Contains( "/EWKWPlus2Jets_WToLNu_M" ) )                            { priority = 2; return "wjj"; }
    if ( dsname.Contains( "/EWKWMinus2Jets_WToLNu_M" ) )                           { priority = 2; return "wjj"; }
    if ( dsname.Contains( "/EWKZ2Jets_ZToLL_M" ) )                                 { priority = 2; return "zjj"; }
    if ( dsname.Contains( "/ZZTo" ) )                                              { priority = 1; return "zz"; }
    if ( dsname.Contains( "/ttZJets_13TeV_madgraphMLM" ) )                         { priority = 1; return "ttz"; }
    if ( dsname.Contains( "/TTZToLLNuNu_M" ) )                                     { priority = 2; return "ttznlo"; }
    if ( dsname.Contains( "/ttWJets_13TeV_madgraphMLM" ) )                         { priority = 1; return "ttw"; }
    if ( dsname.Contains( "/TTWJetsTo" ) )                                         { priority = 2; return "ttwnlo"; }
    if ( dsname.Contains( "/ttH_HToZZ_4LFilter" ) )                                { priority = 2; return "tthzz"; }
    if ( dsname.Contains( "/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV" ) )    { priority = 1; return "tth"; }
    if ( dsname.Contains( "/ttHTobb_M125_13TeV_powheg_pythia8" ) )                 { priority = 1; return "tth"; }
    if ( dsname.Contains( "/TTGJets_TuneCUETP8M1_13TeV" ) )                        { priority = 1; return "ttg"; }
    if ( dsname.Contains( "/TTTT" ) )                                              { priority = 2; return "tttt"; }
    if ( dsname.Contains( "/WGToLNuG" ) && dsname.Contains( "madgraphMLM" ) )      { priority = 2; return "wg"; }
    if ( dsname.Contains( "/WGJets_"  ) )                                          { priority = 2; return "wgpt40"; }
    if ( dsname.Contains( "/WGToLNuG" ) && dsname.Contains( "amcatnloFXFX" ) )     { priority = 2; return "wgnlo"; }
    if ( dsname.Contains( "/WGstarTo" ) )                                          { priority = 2; return "wgstar"; }
    if ( dsname.Contains( "/ZGTo2LG"  ) )                                          { priority = 2; return "zg"; }
    if ( dsname.Contains( "/WminusH_HToBB" ) )                                     { priority = 2; return "whbb"; }
    if ( dsname.Contains( "/WplusH_HToBB" ) )                                      { priority = 2; return "whbb"; }
    if ( dsname.Contains( "/VHToNonbb" ) && !passWHWWW() )                         { priority = 1; return "vh"; }
    if ( dsname.Contains( "/VHToNonbb" ) && passWHWWW() )                          { priority = 1; return "whwww"; }
    if ( dsname.Contains( "/WWW" ) )                                               { priority = 2; return "www_incl"; }
    if ( dsname.Contains( "/WWZ" ) )                                               { priority = 1; return "wwz_incl"; }
    if ( dsname.Contains( "/WZZ" ) )                                               { priority = 1; return "wzz_incl"; }
    if ( dsname.Contains( "/ZZZ" ) )                                               { priority = 1; return "zzz_incl"; }
    if ( dsname.Contains( "/WZG" ) )                                               { priority = 2; return "wzg_incl"; }
    if ( dsname.Contains( "/WWG" ) )                                               { priority = 2; return "wwg_incl"; }
    if ( dsname.Contains( "/TEST-www" ) )                                          { priority = 1; return "www"; }
    if ( dsname.Contains( "/TEST-tth" ) )                                          { priority = 2; return "tthmia"; }
    if ( dsname.Contains( "Run2016" ) && dsname.Contains( "DoubleMuon" ) )         { priority = 1; return "data_mm"; }
    if ( dsname.Contains( "Run2016" ) && dsname.Contains( "DoubleEG" ) )           { priority = 1; return "data_ee"; }
    if ( dsname.Contains( "Run2016" ) && dsname.Contains( "MuonEG" ) )             { priority = 1; return "data_em"; }
    if ( dsname.Contains( "Run2016" ) )                                            { priority = 1; return "data_other"; }
    if ( wwwbaby.evt() == 1532413419 && wwwbaby.lumi() == 962 && wwwbaby.run() == 281797 ) { priority = 1; return "data_ee"; }

    if ( dsname.Length() == 0 )
    {
        std::cout << std::endl;
        std::cout << "evt_dataset is empty! ";
        std::cout << wwwbaby.evt() << " " << wwwbaby.lumi() << " " << wwwbaby.run() << std::endl;
    }

    TString othername = dsname + "_UNCATEGORIZED";
    priority = 2;
    return othername.ReplaceAll( "/", "_" );
}

//______________________________________________________________________________________
TString bkgCategory()
{

    if ( wwwbaby.isData() ) return sampleCategory();

    const TString& dsname = wwwbaby.evt_dataset()[0];
    if ( dsname.Contains( "/VHToNonbb" ) && passWHWWW() ) return "whwww";
    if ( dsname.Contains( "/TEST-www" ) )                 return "www";

    std::vector<unsigned int> iSS = lepidx["TightLepton"];
    std::vector<unsigned int> iaSS = lepidx["LbntLepton"];
    std::vector<unsigned int> ilSS = lepidx["LooseLepton"];

    TString sn = "";
    TString sn2 = "";

    if ( ( ( iSS.size() + iaSS.size() ) >= 2 ) || ( ( iSS.size() + ilSS.size() ) >= 2 ) )
    {
        int l1( -1 ), l2( -1 );

        if ( iSS.size() >= 2 )
        {
            l1 = iSS[0];
            l2 = iSS[1];
        }
        else if ( iSS.size() == 1 && iaSS.size() >= 1 )
        {
            l1 = iSS[0];
            l2 = iaSS[0];
        }
        else if ( iSS.size() == 1 && ilSS.size() >= 1 )
        {
            l1 = iSS[0];
            l2 = ilSS[0];
        }
        else if ( iaSS.size() >= 2 )
        {
            l1 = iaSS[0];
            l2 = iaSS[1];
        }
        else if ( ilSS.size() >= 2 )
        {
            l1 = ilSS[0];
            l2 = ilSS[1];
        }

        int nW( 0 ), nZ( 0 ), nG( 0 ), nF( 0 );

        if ( wwwbaby.lep_isFromW()[l1] )
            ++nW;
        else if ( wwwbaby.lep_isFromZ()[l1] )
            ++nZ;
        else if ( wwwbaby.lep_isFromB()[l1] || wwwbaby.lep_isFromC()[l1] || wwwbaby.lep_isFromL()[l1] || wwwbaby.lep_isFromLF()[l1] )
            ++nF;
        else if ( wwwbaby.lep_motherIdSS()[l1] == ( -3 ) )
            ++nG;

        if ( wwwbaby.lep_isFromW()[l2] )
            ++nW;
        else if ( wwwbaby.lep_isFromZ()[l2] )
            ++nZ;
        else if ( wwwbaby.lep_isFromB()[l2] || wwwbaby.lep_isFromC()[l2] || wwwbaby.lep_isFromL()[l2] || wwwbaby.lep_isFromLF()[l2] )
            ++nF;
        else if ( wwwbaby.lep_motherIdSS()[l2] == ( -3 ) )
            ++nG;

        if ( nW == 2 && wwwbaby.lep_mc_Id()[l1]*wwwbaby.lep_mc_Id()[l2] > 0 )
            sn = "trueSS";//W+W+
        else if ( nW == 2 )
            sn = "chargeflips";//W+W-
        else if ( nZ == 2 && wwwbaby.lep_mc_Id()[l1]*wwwbaby.lep_mc_Id()[l2] <= 0 )
            sn = "chargeflips";//Z
        else if ( nZ == 2 )
            sn = "SSLL";//ZZ both with a lost lepton
        else if ( nW == 1 && nZ == 1 )
            sn = "SSLL";//WZ
        else if ( ( nW + nZ ) == 1 && nG == 1 )
            sn = "photonfakes";
        else if ( ( nW + nZ ) == 1 )
            sn = "fakes";
        else if ( ( nW + nZ ) == 0 && nG == 2 )
            sn = "photondoublefakes";
        else if ( ( nW + nZ ) == 0 && nG == 1 )
            sn = "fakesphotonfakes";
        else if ( ( nW + nZ ) == 0 )
            sn = "doublefakes";
        else
        {
            if ( nG >= 1 )
                sn = "otherphotonfakes";
            else
                sn = "others";
        }
    }
    else
    {
        std::cout << std::endl;
        std::cout << "bkgCategoryFull() Why am I here?";
        std::cout << wwwbaby.evt() << " " << wwwbaby.lumi() << " " << wwwbaby.run() << std::endl;
        return "others";
    }

    return sn;

//    if ( ( i3l.size() >= 3 ) || ( i3l.size() == 2 && looseEle >= 0 ) )
//    {
//        int l1( -1 ), l2( -1 ), l3( -1 );
//
//        if ( i3l.size() >= 3 )
//        {
//            l1 = i3l[0];
//            l2 = i3l[1];
//            l3 = i3l[2];
//        }
//
//        if ( i3l.size() == 2 && looseEle >= 0 )
//        {
//            l1 = i3l[0];
//            l2 = i3l[1];
//            l3 = looseEle;
//        }
//
//        int nW( 0 ), nZ( 0 ), nG( 0 ), nF( 0 );
//
//        if ( wwwbaby.lep_isFromW()[l1] )
//            ++nW;
//        else if ( wwwbaby.lep_isFromZ()[l1] )
//            ++nZ;
//        else if ( wwwbaby.lep_isFromB()[l1] || wwwbaby.lep_isFromC()[l1] || wwwbaby.lep_isFromL()[l1] || wwwbaby.lep_isFromLF()[l1] )
//            ++nF;
//        else if ( wwwbaby.lep_motherIdSS()[l1] == ( -3 ) )
//            ++nG;
//
//        if ( wwwbaby.lep_isFromW()[l2] )
//            ++nW;
//        else if ( wwwbaby.lep_isFromZ()[l2] )
//            ++nZ;
//        else if ( wwwbaby.lep_isFromB()[l2] || wwwbaby.lep_isFromC()[l2] || wwwbaby.lep_isFromL()[l2] || wwwbaby.lep_isFromLF()[l2] )
//            ++nF;
//        else if ( wwwbaby.lep_motherIdSS()[l2] == ( -3 ) )
//            ++nG;
//
//        if ( wwwbaby.lep_isFromW()[l3] )
//            ++nW;
//        else if ( wwwbaby.lep_isFromZ()[l3] )
//            ++nZ;
//        else if ( wwwbaby.lep_isFromB()[l3] || wwwbaby.lep_isFromC()[l3] || wwwbaby.lep_isFromL()[l3] || wwwbaby.lep_isFromLF()[l3] )
//            ++nF;
//        else if ( wwwbaby.lep_motherIdSS()[l3] == ( -3 ) )
//            ++nG;
//
//        if ( nW == 3 && ( wwwbaby.lep_mc_Id()[l1] > 0 && wwwbaby.lep_mc_Id()[l2] > 0 && wwwbaby.lep_mc_Id()[l3] > 0 ) )
//            sn2 = "chargeflips";//W+W+W+ - it could be +++ final state, but at the end this final state will be vetoed, so if reco is ++- (e.g.), then this is a chargeflip
//        else if ( nW == 3 && ( wwwbaby.lep_mc_Id()[l1] < 0 && wwwbaby.lep_mc_Id()[l2] < 0 && wwwbaby.lep_mc_Id()[l3] < 0 ) )
//            sn2 = "chargeflips";//W+W+W+ - it could be +++ final state, but at the end this final state will be vetoed, so if reco is ++- (e.g.), then this is a chargeflip
//        else if ( nW == 3 )
//            sn2 = "trueWWW";
//        else if ( nW == 2 && nZ == 1 )
//            sn2 = "3lLL";//ttZ w/ LL
//        else if ( nW == 1 && nZ == 2 )
//            sn2 = "true3L";//WZ, neglect WZZ as LL
//        else if ( nZ == 3 )
//            sn2 = "3lLL";//ZZ
//        else if ( ( nW + nZ ) == 2 )
//        {
//            if ( nG == 1 )
//                sn2 = "photonfakes";
//            else
//                sn2 = "fakes";
//        }
//        else if ( ( nW + nZ ) == 1 )
//        {
//            if ( nG == 2 )
//                sn2 = "photondoublefakes";
//            else if ( nG == 1 )
//                sn2 = "fakesphotonfakes";
//            else
//                sn2 = "doublefakes";
//        }
//        else
//        {
//            if ( nG == 3 )
//                sn2 = "photontriplefakes";
//            else if ( nG >= 1 )
//                sn2 = "otherphotonfakes";
//            else
//                sn2 = "others";//could be triple fakes
//        }
//    }
}

//______________________________________________________________________________________
void printEventID()
{
    std::cout << std::endl;
    std::cout << "(evt, run, lumi) = (" << wwwbaby.evt() << ", " << wwwbaby.run() << ", " << wwwbaby.lumi() << ")" << std::endl;
    std::cout << wwwbaby.run() << ":" << wwwbaby.lumi() << ":" << wwwbaby.evt() << std::endl;
}

