#include "WWWTools.h"

unsigned int objidx_set_to_eventid = 0;
int objidx_set_to_run = 0;
int objidx_set_to_lumi = 0;
ObjIdx lepidx;
ObjIdx jetidx;

unsigned int genobjidx_set_to_eventid = 0;
int genobjidx_set_to_run = 0;
int genobjidx_set_to_lumi = 0;
ObjIdx genidx;

unsigned int sample_category_set_to_eventid = 0;
int sample_category_set_to_run = 0;
int sample_category_set_to_lumi = 0;
TString sample_category;
int sample_priority;

unsigned int bkg_category_set_to_eventid = 0;
int bkg_category_set_to_run = 0;
int bkg_category_set_to_lumi = 0;
TString bkg_category;

#define MZ 91.1876
#define LUMI 35.87

//_________________________________________________________________________________________________
bool passPresel()
{
    if (
        // Same sign two lepton preselection
        (
         ( lepidx["TightLepton"].size() == 2 ) &&
         ( lepidx["GoodSSJet"].size() >= 2 ) &&
         ( lepidx["VetoLepton"].size() == 2 ) &&
         ( ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]] * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]] ) > 0 ) &&
         ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt() >= 30.)  &&
         ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt() >= 30.) )
        ||
        // Three lepton preselection
        (
         ( lepidx["TightLepton"].size() == 3                      )  &&
         ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt() >= 25. )  &&
         ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt() >= 20. )  &&
         ( wwwbaby.lep_p4()[lepidx["TightLepton"][2]].pt() >= 20. ) )
        ||
        // One loose but not tight + one tight preselection
        (
         ( lepidx["TightLepton"].size() == 1 ) &&
         ( lepidx["LooseLepton"].size() == 2 ) &&
         ( lepidx["GoodSSJet"].size() >= 2 ) &&
         ( lepidx["VetoLepton"].size() == 2 ) &&
         ( ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]] * wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]] ) > 0 ) &&
         ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]].pt() >= 30.)  &&
         ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]].pt() >= 30.) )
        ||
        // One loose but not tight + two tight preselection
        (
         ( lepidx["TightLepton"].size() == 2 ) &&
         ( lepidx["LooseLepton"].size() == 3 ) &&
         ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]].pt() >= 25.)  &&
         ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]].pt() >= 20.)  &&
         ( wwwbaby.lep_p4()[lepidx["LooseLepton"][2]].pt() >= 20.) )
       )
        return true;
    else
        return false;
}

//______________________________________________________________________________________
bool passSSMM( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passTrigMM()                                                             )) return false;
    if (!( lepidx[lepid].size()                              ==   2                 )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           == 169                 )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                                            )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           >    0                 )) return false;
    if (!( jetidx["GoodSSJet"].size()                        >=   2                 )) return false;
    if (!( lepidx["VetoLepton"].size()                       ==   2                 )) return false;
    if (!( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0                 )) return false;
    if (!( jetidx["LooseBJet"].size()                        ==   0     || dropbtag )) return false;
    if (!( ( MjjW()                                          >   60. &&
             MjjW()                                          <  100. )  || dropmjj  )) return false;
    if (!( MjjLead()                                         <  400.    || dropjet  )) return false;
    if (!( DEtajjLead()                                      <    1.5   || dropjet  )) return false;
    if (!( ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   >   40.                )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSEM( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passTrigEM()                                                             )) return false;
    if (!( lepidx[lepid].size()                              ==   2                 )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           == 143                 )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                                            )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           >    0                 )) return false;
    if (!( jetidx["GoodSSJet"].size()                        >=   2                 )) return false;
    if (!( lepidx["VetoLepton"].size()                       ==   2                 )) return false;
    if (!( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0                 )) return false;
    if (!( jetidx["LooseBJet"].size()                        ==   0     || dropbtag )) return false;
    if (!( ( MjjW()                                          >   60. &&
             MjjW()                                          <  100. )  || dropmjj  )) return false;
    if (!( MjjLead()                                         <  400.    || dropjet  )) return false;
    if (!( DEtajjLead()                                      <    1.5   || dropjet  )) return false;
    if (!( ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   >   30.                )) return false;
    if (!( wwwbaby.met_pt()                                  >   40.                )) return false;
    if (!( MTmax( lepid )                                    >   90.                )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSEE( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passTrigEE()                                                              )) return false;
    if (!( lepidx[lepid].size()                              ==   2                  )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           == 121                  )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                                             )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt()           >   30.                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt()           >   30.                 )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           >    0                  )) return false;
    if (!( jetidx["GoodSSJet"].size()                        >=   2                  )) return false;
    if (!( lepidx["VetoLepton"].size()                       ==   2                  )) return false;
    if (!( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0                  )) return false;
    if (!( jetidx["LooseBJet"].size()                        ==   0      || dropbtag )) return false;
    if (!( ( MjjW()                                          >   60. &&
             MjjW()                                          <  100. )   || dropmjj  )) return false;
    if (!( MjjLead()                                         <  400.     || dropjet  )) return false;
    if (!( DEtajjLead()                                      <    1.5    || dropjet  )) return false;
    if (!( ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   >   40.                 )) return false;
    if (!( wwwbaby.met_pt()                                  >   40.                 )) return false;
    if (!( ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   <   80. ||
           ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   >  100.    )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSMMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passSSMM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passSSMM();
        else
            return false;
    }
    else
    {
        return passSSMM();
    }
}

//______________________________________________________________________________________
bool passSSEMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passSSEM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passSSEM();
        else
            return false;
    }
    else
    {
        return passSSEM();
    }
}

//______________________________________________________________________________________
bool passSSEEPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passSSEE();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passSSEE();
        else
            return false;
    }
    else
    {
        return passSSEE();
    }
}

//______________________________________________________________________________________
bool passSSARMMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passSSARMM();
    else
        return false;
}

//______________________________________________________________________________________
bool passSSAREMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passSSAREM();
    else
        return false;
}

//______________________________________________________________________________________
bool passSSAREEPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passSSAREE();
    else
        return false;
}

//______________________________________________________________________________________
bool passSSARMM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1 )) return false;
    if (!( passSSMM( "LooseLepton" )         )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSAREM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1 )) return false;
    if (!( passSSEM( "LooseLepton" )         )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSAREE()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1 )) return false;
    if (!( passSSEE( "LooseLepton" )         )) return false;
    return true;
}

//______________________________________________________________________________________
bool passBTagVRSSMM()
{
    setObjectIndices();
    if (!( passSSMM( "TightLepton", true ) )) return false;
    if (!( jetidx["LooseBJet"].size() >= 1 )) return false;
    return true;
}

//______________________________________________________________________________________
bool passBTagVRSSEM()
{
    setObjectIndices();
    if (!( passSSEM( "TightLepton", true ) )) return false;
    if (!( jetidx["LooseBJet"].size() >= 1 )) return false;
    return true;
}

//______________________________________________________________________________________
bool passBTagVRSSEE()
{
    setObjectIndices();
    if (!( passSSEE( "TightLepton", true ) )) return false;
    if (!( jetidx["LooseBJet"].size() >= 1 )) return false;
    return true;
}

//______________________________________________________________________________________
bool passBTagARSSMM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1 )) return false;
    if (!( passSSMM( "LooseLepton", true )   )) return false;
    if (!( jetidx["LooseBJet"].size()   >= 1 )) return false;
    return true;
}

//______________________________________________________________________________________
bool passBTagARSSEM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1 )) return false;
    if (!( passSSEM( "LooseLepton", true )   )) return false;
    if (!( jetidx["LooseBJet"].size()   >= 1 )) return false;
    return true;
}

//______________________________________________________________________________________
bool passBTagARSSEE()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1 )) return false;
    if (!( passSSEE( "LooseLepton", true )   )) return false;
    if (!( jetidx["LooseBJet"].size()   >= 1 )) return false;
    return true;
}


//______________________________________________________________________________________
bool passBTagVRSSMMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passBTagVRSSMM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passBTagVRSSMM();
        else
            return false;
    }
    else
    {
        return passBTagVRSSMM();
    }
}

//______________________________________________________________________________________
bool passBTagVRSSEMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passBTagVRSSEM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passBTagVRSSEM();
        else
            return false;
    }
    else
    {
        return passBTagVRSSEM();
    }
}

//______________________________________________________________________________________
bool passBTagVRSSEEPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passBTagVRSSEE();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passBTagVRSSEE();
        else
            return false;
    }
    else
    {
        return passBTagVRSSEE();
    }
}

//______________________________________________________________________________________
bool passBTagARSSMMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passBTagARSSMM();
    else
        return false;
}

//______________________________________________________________________________________
bool passBTagARSSEMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passBTagARSSEM();
    else
        return false;
}

//______________________________________________________________________________________
bool passBTagARSSEEPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passBTagARSSEE();
    else
        return false;
}


//______________________________________________________________________________________
bool passMjjSBVRSSMM()
{
    setObjectIndices();
    if (!( passSSMM( "TightLepton", false, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBVRSSEM()
{
    setObjectIndices();
    if (!( passSSEM( "TightLepton", false, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBVRSSEE()
{
    setObjectIndices();
    if (!( passSSEE( "TightLepton", false, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBARSSMM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSMM( "LooseLepton", false, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBARSSEM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSEM( "LooseLepton", false, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBARSSEE()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSEE( "LooseLepton", false, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBVRSSMMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passMjjSBVRSSMM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passMjjSBVRSSMM();
        else
            return false;
    }
    else
    {
        return passMjjSBVRSSMM();
    }
}

//______________________________________________________________________________________
bool passMjjSBVRSSEMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passMjjSBVRSSEM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passMjjSBVRSSEM();
        else
            return false;
    }
    else
    {
        return passMjjSBVRSSEM();
    }
}

//______________________________________________________________________________________
bool passMjjSBVRSSEEPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passMjjSBVRSSEE();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passMjjSBVRSSEE();
        else
            return false;
    }
    else
    {
        return passMjjSBVRSSEE();
    }
}

//______________________________________________________________________________________
bool passMjjSBARSSMMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passMjjSBARSSMM();
    else
        return false;
}

//______________________________________________________________________________________
bool passMjjSBARSSEMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passMjjSBARSSEM();
    else
        return false;
}

//______________________________________________________________________________________
bool passMjjSBARSSEEPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passMjjSBARSSEE();
    else
        return false;
}

//______________________________________________________________________________________
bool passMjjSBPRVRSSMM()
{
    setObjectIndices();
    if (!( passSSMM( "TightLepton", false, true, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBPRVRSSEM()
{
    setObjectIndices();
    if (!( passSSEM( "TightLepton", false, true, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBPRVRSSEE()
{
    setObjectIndices();
    if (!( passSSEE( "TightLepton", false, true, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBPRARSSMM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSMM( "LooseLepton", false, true, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBPRARSSEM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSEM( "LooseLepton", false, true, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBPRARSSEE()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSEE( "LooseLepton", false, true, true ) )) return false;
    if (!( MjjW() <= 60. || MjjW() >= 100.        )) return false;
    return true;
}

//______________________________________________________________________________________
bool passMjjSBPRVRSSMMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passMjjSBPRVRSSMM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passMjjSBPRVRSSMM();
        else
            return false;
    }
    else
    {
        return passMjjSBPRVRSSMM();
    }
}

//______________________________________________________________________________________
bool passMjjSBPRVRSSEMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passMjjSBPRVRSSEM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passMjjSBPRVRSSEM();
        else
            return false;
    }
    else
    {
        return passMjjSBPRVRSSEM();
    }
}

//______________________________________________________________________________________
bool passMjjSBPRVRSSEEPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passMjjSBPRVRSSEE();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passMjjSBPRVRSSEE();
        else
            return false;
    }
    else
    {
        return passMjjSBPRVRSSEE();
    }
}

//______________________________________________________________________________________
bool passMjjSBPRARSSMMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passMjjSBPRARSSMM();
    else
        return false;
}

//______________________________________________________________________________________
bool passMjjSBPRARSSEMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passMjjSBPRARSSEM();
    else
        return false;
}

//______________________________________________________________________________________
bool passMjjSBPRARSSEEPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passMjjSBPRARSSEE();
    else
        return false;
}

//______________________________________________________________________________________
bool passPRVRSSMM()
{
    setObjectIndices();
    if (!( passSSMM( "TightLepton", false, true, true ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passPRVRSSEM()
{
    setObjectIndices();
    if (!( passSSEM( "TightLepton", false, true, true ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passPRVRSSEE()
{
    setObjectIndices();
    if (!( passSSEE( "TightLepton", false, true, true ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passPRARSSMM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSMM( "LooseLepton", false, true, true ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passPRARSSEM()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSEM( "LooseLepton", false, true, true ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passPRARSSEE()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 1      )) return false;
    if (!( passSSEE( "LooseLepton", false, true, true ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passPRVRSSMMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passPRVRSSMM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passPRVRSSMM();
        else
            return false;
    }
    else
    {
        return passPRVRSSMM();
    }
}

//______________________________________________________________________________________
bool passPRVRSSEMPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passPRVRSSEM();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passPRVRSSEM();
        else
            return false;
    }
    else
    {
        return passPRVRSSEM();
    }
}

//______________________________________________________________________________________
bool passPRVRSSEEPred()
{
    setObjectIndices();
    if ( isFakeSample() )
    {
        if ( isRealLeptonEvent() )
            return passPRVRSSEE();
        else
            return false;
    }
    else if ( isPhotonConversionSample() )
    {
        if ( isPhotonConversionEvent() )
            return passPRVRSSEE();
        else
            return false;
    }
    else
    {
        return passPRVRSSEE();
    }
}

//______________________________________________________________________________________
bool passPRARSSMMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passPRARSSMM();
    else
        return false;
}

//______________________________________________________________________________________
bool passPRARSSEMPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passPRARSSEM();
    else
        return false;
}

//______________________________________________________________________________________
bool passPRARSSEEPred()
{
    setObjectIndices();
    if ( sampleCategory().Contains( "data" ) )
        return passPRARSSEE();
    else
        return false;
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
bool pass3L0SFOS( TString lepid )
{
    setObjectIndices();
    if (!( wwwbaby.evt_passgoodrunlist()                    )) return false;
    if (!( wwwbaby.firstgoodvertex()               ==   0   )) return false;
    if (!( lepidx[lepid].size()                    ==   3   )) return false;
    if (!( lepidx["VetoLepton"].size()             ==   3   )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt() >   25.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt() >   20.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt() >   20.  )) return false;
    if (!( abs( totalCharge() )                    ==   1   )) return false;
    if (!( jetidx["Good3LJet"].size()              <=   1   )) return false;
    if (!( Pt3l()                                  >   60.  )) return false;
    if (!( jetidx["LooseBJet"].size()              ==   0   )) return false;
    if (!( getNumSFOS()                            ==   0   )) return false;
    if (!( get0SFOSMll()                           >   20.  )) return false;
    if (!( fabs( get0SFOSMee() - MZ )              >   15.  )) return false;
    if (!( DPhi3lMET()                             >    2.7 )) return false;
    return true;
}

//______________________________________________________________________________________
// The 3L1SFOS signal region definition
bool pass3L1SFOS( TString lepid )
{
    setObjectIndices();
    if (!( wwwbaby.evt_passgoodrunlist()                    )) return false;
    if (!( wwwbaby.firstgoodvertex()               ==   0   )) return false;
    if (!( lepidx[lepid].size()                    ==   3   )) return false;
    if (!( lepidx["VetoLepton"].size()             ==   3   )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt() >   25.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt() >   20.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt() >   20.  )) return false;
    if (!( abs( totalCharge() )                    ==   1   )) return false;
    if (!( jetidx["Good3LJet"].size()              <=   1   )) return false;
    if (!( Pt3l()                                  >   60.  )) return false;
    if (!( jetidx["LooseBJet"].size()              ==   0   )) return false;
    if (!( getNumSFOS()                            ==   1   )) return false;
    if (!( wwwbaby.met_pt()                        >   45.  )) return false;
    if (!( get1SFOSMll()                           <   55. ||
           get1SFOSMll()                           >  100.  )) return false;
    if (!( DPhi3lMET()                             >    2.5 )) return false;
    return true;
}

//______________________________________________________________________________________
// The 3L2SFOS signal region definition
bool pass3L2SFOS( TString lepid )
{
    setObjectIndices();
    if (!( wwwbaby.evt_passgoodrunlist()                    )) return false;
    if (!( wwwbaby.firstgoodvertex()               ==   0   )) return false;
    if (!( lepidx[lepid].size()                    ==   3   )) return false;
    if (!( lepidx["VetoLepton"].size()             ==   3   )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt() >   25.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt() >   20.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt() >   20.  )) return false;
    if (!( abs( totalCharge() )                    ==   1   )) return false;
    if (!( jetidx["Good3LJet"].size()              <=   1   )) return false;
    if (!( Pt3l()                                  >   60.  )) return false;
    if (!( jetidx["LooseBJet"].size()              ==   0   )) return false;
    if (!( getNumSFOS()                            ==   2   )) return false;
    if (!( wwwbaby.met_pt()                        >   55.  )) return false;
    if (!( fabs( get2SFOSMll0() - MZ )             >   20.  )) return false;
    if (!( fabs( get2SFOSMll1() - MZ )             >   20.  )) return false;
    if (!( DPhi3lMET()                             >    2.5 )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3LAR0SFOS()
{
    if (!( lepidx["Tight3lLepton"].size() == 2 )) return false;
    if (!( pass3L0SFOS( "LooseLepton" )        )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3LAR1SFOS()
{
    if (!( lepidx["Tight3lLepton"].size() == 2 )) return false;
    if (!( pass3L1SFOS( "LooseLepton" )        )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3LAR2SFOS()
{
    if (!( lepidx["Tight3lLepton"].size() == 2 )) return false;
    if (!( pass3L2SFOS( "LooseLepton" )        )) return false;
    return true;
}

//______________________________________________________________________________________
//bool passTrigMM() { return wwwbaby.HLT_DoubleMu() || wwwbaby.HLT_DoubleMu_noiso(); }
//bool passTrigEM() { return wwwbaby.HLT_MuEG() || wwwbaby.HLT_MuEG_noiso(); }
//bool passTrigEE() { return wwwbaby.HLT_DoubleEl() || wwwbaby.HLT_DoubleEl_DZ() || wwwbaby.HLT_DoubleEl_DZ_2() || wwwbaby.HLT_DoubleEl_noiso(); }
bool passTrigMM() { return wwwbaby.HLT_DoubleMu(); }
bool passTrigEM() { return wwwbaby.HLT_MuEG(); }
bool passTrigEE() { return wwwbaby.HLT_DoubleEl() || wwwbaby.HLT_DoubleEl_DZ() || wwwbaby.HLT_DoubleEl_DZ_2(); }

//______________________________________________________________________________________
float weight( bool applyfakefactor, int isyst )
{
    float wgt = wwwbaby.isData() ? 1 : wwwbaby.evt_scale1fb() * LUMI;

    if ( wwwbaby.evt_dataset()[0].Contains( "/TEST-www/wwwext-Private80X-v1/USER" ) )
        wgt *= 0.066805*164800./(91900.+164800.);
    if ( wwwbaby.evt_dataset()[0].Contains( "/TEST-www/www-Private80X-v1/USER" ) )
        wgt *= 0.066805*91900./(91900.+164800.);

    int ffsyst = 0;
    // Handle isyst
    if      ( isyst == 0 ) ffsyst =  0;
    else if ( isyst == 1 ) ffsyst = -1;
    else if ( isyst == 2 ) ffsyst =  1;
    else
        std::cerr << "Error: Unrecognized systematics: " << isyst << std::endl;


    if ( applyfakefactor )
    {
        if ( lepidx["TightLepton"].size() == 1 && lepidx["LooseLepton"].size() == 2 )
            wgt *= fakefactor( lepidx["LbntLepton"][0], ffsyst );
        else if ( lepidx["Tight3lLepton"].size() == 2 && lepidx["LooseLepton"].size() == 3 )
            wgt *= fakefactor( lepidx["Lbn3tLepton"][0], ffsyst );
    }

    return wgt;
}

//______________________________________________________________________________________
void printEvent()
{
    std::cout << "=======================================================" << std::endl;
    std::cout << " run  = " << wwwbaby.run() << std::endl;
    std::cout << " lumi = " << wwwbaby.lumi() << std::endl;
    std::cout << " evt  = " << wwwbaby.evt() << std::endl;
    std::cout << " dsname  = " << wwwbaby.evt_dataset()[0] << std::endl;
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
            std::cout << wwwbaby.lep_relIso03EA()[ilep]  << " ";
            std::cout << wwwbaby.lep_relIso03EAv2()[ilep]<< " ";
            std::cout << wwwbaby.lep_ip3d()[ilep]  << " ";
            std::cout << wwwbaby.lep_dxy()[ilep]<< " ";
            std::cout << wwwbaby.lep_dz()[ilep]<< " ";
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
    std::cout  << "passSSMM      " << " : " << passSSMM       () << std::endl;
    std::cout  << "passSSEM      " << " : " << passSSEM       () << std::endl;
    std::cout  << "passSSEE      " << " : " << passSSEE       () << std::endl;
    std::cout  << "passSSARMM    " << " : " << passSSARMM     () << std::endl;
    std::cout  << "passSSAREM    " << " : " << passSSAREM     () << std::endl;
    std::cout  << "passSSAREE    " << " : " << passSSAREE     () << std::endl;
    std::cout  << "passBTagVRSSMM" << " : " << passBTagVRSSMM () << std::endl;
    std::cout  << "passBTagVRSSEM" << " : " << passBTagVRSSEM () << std::endl;
    std::cout  << "passBTagVRSSEE" << " : " << passBTagVRSSEE () << std::endl;
    std::cout  << "passBTagARSSMM" << " : " << passBTagARSSMM () << std::endl;
    std::cout  << "passBTagARSSEM" << " : " << passBTagARSSEM () << std::endl;
    std::cout  << "passBTagARSSEE" << " : " << passBTagARSSEE () << std::endl;
    std::cout  << "pass3L0SFOS   " << " : " << pass3L0SFOS    () << std::endl;
    std::cout  << "pass3L1SFOS   " << " : " << pass3L1SFOS    () << std::endl;
    std::cout  << "pass3L2SFOS   " << " : " << pass3L2SFOS    () << std::endl;
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
void setObjectIndices( bool preload )
{
    // If the object already selected skip
    if ( isObjectSelected() ) return;
    if ( preload )
    {
        loadLeptonIndices();
        loadJetsIndices();
    }
    else
    {
        //otherwise perform object selections
        lepidx = getLeptonsIndices();
        jetidx = getJetsIndices();
    }
    // Now set the event id of the event the objects are selected from
    objidx_set_to_eventid = wwwbaby.evt();
    objidx_set_to_run     = wwwbaby.run();
    objidx_set_to_lumi    = wwwbaby.lumi();
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
ObjIdx makeEmptyLepidx()
{
    ObjIdx idx;
    idx["Tight3lLepton"].clear();
    idx["Tight3lElec"].clear();
    idx["TightLepton"].clear();
    idx["TightMuon"].clear();
    idx["TightElec"].clear();
    idx["LooseLepton"].clear();
    idx["LooseMuon"].clear();
    idx["LooseElec"].clear();
    idx["VetoLepton"].clear();
    idx["VetoMuon"].clear();
    idx["VetoElec"].clear();
    idx["LbntLepton"].clear();
    idx["LbntMuon"].clear();
    idx["LbntElec"].clear();
    idx["Lbnt3lLepton"].clear();
    idx["Lbnt3lMuon"].clear();
    idx["Lbnt3lElec"].clear();
    idx["SignalLepton"].clear();
    idx["SignalMuon"].clear();
    idx["SignalElec"].clear();
    return idx;
}

//______________________________________________________________________________________
ObjIdx makeEmptyJetidx()
{
    ObjIdx idx;
    idx["GoodSSJet"].clear();
    idx["Good3LJet"].clear();
    idx["LooseBJet"].clear();
    idx["GoodSSWJet"].clear();
    return idx;
}

//______________________________________________________________________________________
ObjIdx getLeptonsIndices()
{
    ObjIdx idx = makeEmptyLepidx();
    for ( unsigned int ilep = 0; ilep < wwwbaby.lep_p4().size(); ++ilep )
    {
        if ( isTight3lLepton( ilep ) ) idx["Tight3lLepton"].push_back( ilep );
        if ( isTight3lElec  ( ilep ) ) idx["Tight3lElec"]  .push_back( ilep );
        if ( isTightLepton  ( ilep ) ) idx["TightLepton"]  .push_back( ilep );
        if ( isTightMuon    ( ilep ) ) idx["TightMuon"]    .push_back( ilep );
        if ( isTightElec    ( ilep ) ) idx["TightElec"]    .push_back( ilep );
        if ( isLooseLepton  ( ilep ) ) idx["LooseLepton"]  .push_back( ilep );
        if ( isLooseMuon    ( ilep ) ) idx["LooseMuon"]    .push_back( ilep );
        if ( isLooseElec    ( ilep ) ) idx["LooseElec"]    .push_back( ilep );
        if ( isVetoLepton   ( ilep ) ) idx["VetoLepton"]   .push_back( ilep );
        if ( isVetoMuon     ( ilep ) ) idx["VetoMuon"]     .push_back( ilep );
        if ( isVetoElec     ( ilep ) ) idx["VetoElec"]     .push_back( ilep );
        if ( isLbntLepton   ( ilep ) ) idx["LbntLepton"]   .push_back( ilep );
        if ( isLbntMuon     ( ilep ) ) idx["LbntMuon"]     .push_back( ilep );
        if ( isLbntElec     ( ilep ) ) idx["LbntElec"]     .push_back( ilep );
        if ( isLbnt3lLepton ( ilep ) ) idx["Lbnt3lLepton"] .push_back( ilep );
        if ( isLbnt3lMuon   ( ilep ) ) idx["Lbnt3lMuon"]   .push_back( ilep );
        if ( isLbnt3lElec   ( ilep ) ) idx["Lbnt3lElec"]   .push_back( ilep );
    }
    // Default is to set signal leptons to Tight leptons
    idx["SignalLepton"] = idx["LooseLepton"];
    idx["SignalMuon"] = idx["LooseMuon"];
    idx["SignalElec"] = idx["LooseElec"];
    return idx;
}

//______________________________________________________________________________________
void loadLeptonIndices()
{
    lepidx["Tight3lLepton"] = wwwbaby.Tight3lLepton();
    lepidx["Tight3lElec"]   = wwwbaby.Tight3lElec();
    lepidx["TightLepton"]   = wwwbaby.TightLepton();
    lepidx["TightMuon"]     = wwwbaby.TightMuon();
    lepidx["TightElec"]     = wwwbaby.TightElec();
    lepidx["LooseLepton"]   = wwwbaby.LooseLepton();
    lepidx["LooseMuon"]     = wwwbaby.LooseMuon();
    lepidx["LooseElec"]     = wwwbaby.LooseElec();
    lepidx["VetoLepton"]    = wwwbaby.VetoLepton();
    lepidx["VetoMuon"]      = wwwbaby.VetoMuon();
    lepidx["VetoElec"]      = wwwbaby.VetoElec();
    lepidx["LbntLepton"]    = wwwbaby.LbntLepton();
    lepidx["LbntMuon"]      = wwwbaby.LbntMuon();
    lepidx["LbntElec"]      = wwwbaby.LbntElec();
    lepidx["Lbnt3lLepton"]  = wwwbaby.Lbnt3lLepton();
    lepidx["Lbnt3lMuon"]    = wwwbaby.Lbnt3lMuon();
    lepidx["Lbnt3lElec"]    = wwwbaby.Lbnt3lElec();
    lepidx["SignalLepton"]  = wwwbaby.SignalLepton();
    lepidx["SignalMuon"]    = wwwbaby.SignalMuon();
    lepidx["SignalElec"]    = wwwbaby.SignalElec();
}

//______________________________________________________________________________________
ObjIdx getJetsIndices()
{
    ObjIdx idx = makeEmptyJetidx();
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
void loadJetsIndices()
{
    jetidx["GoodSSJet"]  = wwwbaby.GoodSSJet();
    jetidx["Good3LJet"]  = wwwbaby.Good3LJet();
    jetidx["LooseBJet"]  = wwwbaby.LooseBJet();
    jetidx["GoodSSWJet"] = wwwbaby.GoodSSWJet();
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
    if (!( fabs( wwwbaby.lep_p4()[ilep].pt()    )  >  20.0   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
//    if (!(       wwwbaby.lep_relIso03EAv2()[ilep]        <   0.06  )) return false;
//    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )        <   0.015 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isTightElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 11     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].pt()    )  >  20.0   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6   )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0     )) return false;
    if (!(       wwwbaby.lep_isTriggerSafe_v1()[ilep]        )) return false;
//    if ( fabs( wwwbaby.lep_etaSC()[ilep] ) <= 1.479 )
//    {
//        if (!(   wwwbaby.lep_relIso03EAv2()[ilep]  <   0.0588  )) return false;
//    }
//    else
//    {
//        if (!(   wwwbaby.lep_relIso03EAv2()[ilep]  <   0.0571  )) return false;
//    }
//    if (!(       wwwbaby.lep_relIso03EAv2()[ilep]  <   0.06  )) return false;
//    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
//    if (!(       wwwbaby.lep_lostHits()[ilep]      ==  0     )) return false;
    return true;
}

//______________________________________________________________________________________
bool isTight3lLepton( int ilep )
{
    return isTightMuon( ilep ) || isTight3lElec( ilep );
}

//______________________________________________________________________________________
bool isTight3lElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 11     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].pt()    )  >  20.0   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6   )) return false;
    if (!(       wwwbaby.lep_isTriggerSafe_v1()[ilep]        )) return false;
//    if ( fabs( wwwbaby.lep_etaSC()[ilep] ) <= 1.479 )
//    {
//        if (!(   wwwbaby.lep_relIso03EAv2()[ilep]  <   0.0588  )) return false;
//    }
//    else
//    {
//        if (!(   wwwbaby.lep_relIso03EAv2()[ilep]  <   0.0571  )) return false;
//    }
//    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
//    if (!(       wwwbaby.lep_lostHits()[ilep]      ==  0     )) return false;
    return true;
}

//______________________________________________________________________________________
bool isVetoLepton( int ilep )
{
    return isVetoMuon( ilep ) || isVetoElec( ilep );
}

//______________________________________________________________________________________
bool isVetoMuon( int ilep )
{
    // https://github.com/cmstas/CORE/blob/1a7909bdb178a1bea19ed812e695d17154e231b9/MuonSelections.cc#L846
    if (!(       wwwbaby.lep_pass_VVV_cutbased_veto()[ilep]  )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 13     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].pt()    )  >  10.0   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    return true;
}

//______________________________________________________________________________________
bool isVetoElec( int ilep )
{
    // https://github.com/cmstas/CORE/blob/1a7909bdb178a1bea19ed812e695d17154e231b9/ElectronSelections.cc#L1527
    if (!(       wwwbaby.lep_pass_VVV_cutbased_veto()[ilep]  )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 11     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].pt()    )  >  10.0   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6   )) return false;
//    if (!(       wwwbaby.lep_3ch_agree()[ilep]       !=  0        )) return false;
//    if (!(       wwwbaby.lep_relIso03EAv2()[ilep]    <   0.4      )) return false;
//    if (!(       wwwbaby.lep_isTriggerSafe_v1()[ilep]             )) return false;
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
    if (!(       wwwbaby.lep_pass_VVV_cutbased_fo()[ilep]        )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 13         )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].pt()    )  >  20.0       )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4       )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015     )) return false;
//    if (!(       wwwbaby.lep_relIso03EAv2()[ilep]  <   0.4   )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_fo()[ilep]       )) return false;
    if (!(  abs( wwwbaby.lep_pdgId()[ilep]      )  == 11        )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].pt()    )  >  20.0      )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5      )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6      )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0        )) return false;
    if (!(       wwwbaby.lep_lostHits()[ilep]      ==  0        )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015    )) return false;
    if (!(       wwwbaby.lep_isTriggerSafe_v1()[ilep]           )) return false;
//    if (!(       wwwbaby.lep_relIso03EAv2()[ilep]  <   0.2      )) return false;
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
bool isLbnt3lLepton( int ilep )
{
    return isLbnt3lMuon( ilep ) || isLbnt3lElec( ilep );
}

//______________________________________________________________________________________
bool isLbnt3lMuon( int ilep )
{
    return isLooseMuon( ilep ) && !isTightMuon( ilep );
}

//______________________________________________________________________________________
bool isLbnt3lElec( int ilep )
{
    return isLooseElec( ilep ) && !isTight3lElec( ilep );
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
void setGenObjectIndices()
{
    // If the object already selected skip
    if ( isGenObjectSelected() ) return;
    //otherwise perform object selections
    std::cout << std::endl;
    std::cout << wwwbaby.ngenLep() << " " << bkgCategory() << " " << sampleCategory() << std::endl;
    for ( unsigned int igen = 0; igen < wwwbaby.genPart_pdgId().size(); ++igen )
    {
        if ( abs( wwwbaby.genPart_motherId()[igen] ) != 23 &&
             abs( wwwbaby.genPart_motherId()[igen] ) != 24 ) 
            continue;
        std::cout << wwwbaby.genPart_p4()[igen].pt()  << " ";
        std::cout << wwwbaby.genPart_p4()[igen].eta() << " ";
        std::cout << wwwbaby.genPart_p4()[igen].phi() << " ";
        std::cout << wwwbaby.genPart_pdgId()[igen] << " ";
        std::cout << wwwbaby.genPart_motherId()[igen] << " ";
        std::cout << std::endl;
    }
    // Now set the event id of the event the objects are selected from
    genobjidx_set_to_eventid = wwwbaby.evt();
    genobjidx_set_to_run     = wwwbaby.run();
    genobjidx_set_to_lumi    = wwwbaby.lumi();
}

//______________________________________________________________________________________
bool isGenObjectSelected()
{
    if ( genobjidx_set_to_eventid == wwwbaby.evt() &&
         genobjidx_set_to_run     == wwwbaby.run() &&
         genobjidx_set_to_lumi    == wwwbaby.lumi() )
        return true;
    else
        return false;
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
    if ( path.Contains( "v0.1.16" ) ) return 16;
    if ( path.Contains( "v0_1_16" ) ) return 16;
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
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]]
             + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).pt();
}

//______________________________________________________________________________________
float DPhi3lMET()
{
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    TLorentzVector met;
    met.SetPtEtaPhiM( wwwbaby.met_pt(), 0, wwwbaby.met_phi(), 0 );
    LorentzVector metlv;
    metlv.SetPxPyPzE( met.Px(), met.Py(), met.Pz(), met.E() );
    return fabs( ROOT::Math::VectorUtil::DeltaPhi(
                wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]]
                + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]],
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
    for ( auto& ilep : lepidx["Tight3lLepton"] )
        charge += wwwbaby.lep_charge()[ilep];
    return charge;
}

//______________________________________________________________________________________
int getNumSFOS()
{
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    /* Loops through pairs of entries in the lep_pdgId vector
       and counts how many have opposite value*/
    int num_SFOS = 0;
    for ( int i = 0; i < ( int ) lepidx["Tight3lLepton"].size(); i++ )
    {
        for ( int j = i + 1; j < ( int ) lepidx["Tight3lLepton"].size(); j++ )
        {
            if ( wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][i]] ==
                    -( wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][j]] ) )
                num_SFOS++;
        }
    }
    return num_SFOS;
}

//______________________________________________________________________________________
float get0SFOSMll()
{
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][2]];
    if ( pdgid0 == pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    TFile* file = 0;
    file->GetName();
    return -999;
}

//______________________________________________________________________________________
float get0SFOSMee()
{
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][2]];
    if ( pdgid0 == pdgid1 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 && abs( pdgid1 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    return -999;
}

//______________________________________________________________________________________
float get1SFOSMll()
{
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    TFile* file = 0;
    file->GetName();
    return -999;
}

//______________________________________________________________________________________
float get2SFOSMll0()
{
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    TFile* file = 0;
    file->GetName();
    return -999;
}

//______________________________________________________________________________________
float get2SFOSMll1()
{
    if ( lepidx["Tight3lLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["Tight3lLepton"][2]];
    if ( pdgid2 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid0 )
        return ( wwwbaby.lep_p4()[lepidx["Tight3lLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["Tight3lLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    TFile* file = 0;
    file->GetName();
    return -999;
}

//______________________________________________________________________________________
float M4()
{
    if ( lepidx["SignalLepton"].size() < 2 )
        return -999;
    if ( jetidx["GoodSSJet"].size() < 2 )
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
bool isSampleCategorySet()
{
    if ( sample_category_set_to_eventid == wwwbaby.evt() &&
         sample_category_set_to_run     == wwwbaby.run() &&
         sample_category_set_to_lumi    == wwwbaby.lumi() )
        return true;
    else
        return false;
}

//______________________________________________________________________________________
TString sampleCategory( int& priority )
{

    if ( isSampleCategorySet() )
    {
        priority = sample_priority;
        return sample_category;
    }

    // Now set the event id of the event the objects are selected from
    sample_category_set_to_eventid = wwwbaby.evt();
    sample_category_set_to_run     = wwwbaby.run();
    sample_category_set_to_lumi    = wwwbaby.lumi();

    const TString& dsname = wwwbaby.evt_dataset()[0];
    if      ( dsname.Contains( "/QCD_HT" ) )                                            { priority = 2; sample_category =  "qcd"; }
    else if ( dsname.Contains( "/GJets" ) )                                             { priority = 2; sample_category =  "gj"; }
    else if ( dsname.Contains( "/WJetsToLNu" ) )                                        { priority = 1; sample_category =  "wj"; }
    else if ( dsname.Contains( "/DYJetsToLL_M" ) && dsname.Contains( "madgraphMLM" ) )  { priority = 1; sample_category =  "dy"; }
    else if ( dsname.Contains( "/DYJetsToLL_M" ) && dsname.Contains( "amcatnloFXFX" ) ) { priority = 2; sample_category =  "dynlo"; }
    else if ( dsname.Contains( "/TTJets_Single" ) && dsname.Contains( "ext1" ) )        { priority = 1; sample_category =  "tt1l"; }
    else if ( dsname.Contains( "/TTJets_Single" ) && !dsname.Contains( "ext1" ) )       { priority = 2; sample_category =  "tt1lnonext"; }
    else if ( dsname.Contains( "/TTJets_Di" ) && dsname.Contains( "ext1" ) )            { priority = 1; sample_category =  "tt2l"; }
    else if ( dsname.Contains( "/TTJets_Di" ) && !dsname.Contains( "ext1") )            { priority = 2; sample_category =  "tt2lnonext"; }
    else if ( dsname.Contains( "/TTTo2L2Nu_" ) )                                        { priority = 2; sample_category =  "tt2lpowheg"; }
    else if ( dsname.Contains( "TTTo2L2Nu_" ) )                                         { priority = 2; sample_category =  "tt2lpowheg"; } // some weird ones have weird names
    else if ( dsname.Contains( "/ST_" ) )                                               { priority = 1; sample_category =  "singletop"; }
    else if ( dsname.Contains( "/tZq" ) )                                               { priority = 1; sample_category =  "singletop"; }
    else if ( dsname.Contains( "/WWTo" ) && dsname.Contains( "powheg" ) )               { priority = 1; sample_category =  "ww"; }
    else if ( dsname.Contains( "/WWTo" ) && dsname.Contains( "DoubleScattering" ) )     { priority = 1; sample_category =  "wwdpi"; }
    else if ( dsname.Contains( "/WpWpJJ_EWK" ) )                                        { priority = 1; sample_category =  "vbsww"; }
    else if ( dsname.Contains( "/WmWmJJ_13" ) )                                         { priority = 2; sample_category =  "wwjj"; }
    else if ( dsname.Contains( "/WpWpJJ_13" ) )                                         { priority = 2; sample_category =  "wwjj"; }
    else if ( dsname.Contains( "/WZTo" ) )                                              { priority = 1; sample_category =  "wz"; }
    else if ( dsname.Contains( "/GluGluHToZZTo4L" ) )                                   { priority = 2; sample_category =  "gghzz"; }
    else if ( dsname.Contains( "/EWKWPlus2Jets_WToLNu_M" ) )                            { priority = 2; sample_category =  "wjj"; }
    else if ( dsname.Contains( "/EWKWMinus2Jets_WToLNu_M" ) )                           { priority = 2; sample_category =  "wjj"; }
    else if ( dsname.Contains( "/EWKZ2Jets_ZToLL_M" ) )                                 { priority = 2; sample_category =  "zjj"; }
    else if ( dsname.Contains( "/ZZTo" ) )                                              { priority = 1; sample_category =  "zz"; }
    else if ( dsname.Contains( "/ttZJets_13TeV_madgraphMLM" ) )                         { priority = 1; sample_category =  "ttz"; }
    else if ( dsname.Contains( "/TTZToLLNuNu_M" ) )                                     { priority = 2; sample_category =  "ttznlo"; }
    else if ( dsname.Contains( "/ttWJets_13TeV_madgraphMLM" ) )                         { priority = 1; sample_category =  "ttw"; }
    else if ( dsname.Contains( "/TTWJetsTo" ) )                                         { priority = 2; sample_category =  "ttwnlo"; }
    else if ( dsname.Contains( "/ttH_HToZZ_4LFilter" ) )                                { priority = 2; sample_category =  "tthzz"; }
    else if ( dsname.Contains( "/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV" ) )    { priority = 1; sample_category =  "tth"; }
    else if ( dsname.Contains( "/ttHTobb_M125_13TeV_powheg_pythia8" ) )                 { priority = 1; sample_category =  "tth"; }
    else if ( dsname.Contains( "/TTGJets_TuneCUETP8M1_13TeV" ) )                        { priority = 1; sample_category =  "ttg"; }
    else if ( dsname.Contains( "/TTTT" ) )                                              { priority = 2; sample_category =  "tttt"; }
    else if ( dsname.Contains( "/WGToLNuG" ) && dsname.Contains( "madgraphMLM" ) )      { priority = 1; sample_category =  "wg"; }
    else if ( dsname.Contains( "/WGJets_"  ) )                                          { priority = 2; sample_category =  "wgpt40"; }
    else if ( dsname.Contains( "/WGToLNuG" ) && dsname.Contains( "amcatnloFXFX" ) )     { priority = 2; sample_category =  "wgnlo"; }
    else if ( dsname.Contains( "/WGstarTo" ) )                                          { priority = 2; sample_category =  "wgstar"; }
    else if ( dsname.Contains( "/ZGTo2LG"  ) )                                          { priority = 1; sample_category =  "zg"; }
    else if ( dsname.Contains( "/WminusH_HToBB" ) )                                     { priority = 2; sample_category =  "whbb"; }
    else if ( dsname.Contains( "/WplusH_HToBB" ) )                                      { priority = 2; sample_category =  "whbb"; }
    else if ( dsname.Contains( "/VHToNonbb" ) && !passWHWWW() )                         { priority = 1; sample_category =  "vh"; }
    else if ( dsname.Contains( "/VHToNonbb" ) && passWHWWW() )                          { priority = 1; sample_category =  "whwww"; }
    else if ( dsname.Contains( "/WWW" ) )                                               { priority = 2; sample_category =  "www_incl"; }
    else if ( dsname.Contains( "/WWZ" ) )                                               { priority = 1; sample_category =  "wwz_incl"; }
    else if ( dsname.Contains( "/WZZ" ) )                                               { priority = 1; sample_category =  "wzz_incl"; }
    else if ( dsname.Contains( "/ZZZ" ) )                                               { priority = 1; sample_category =  "zzz_incl"; }
    else if ( dsname.Contains( "/WZG" ) )                                               { priority = 2; sample_category =  "wzg_incl"; }
    else if ( dsname.Contains( "/WWG" ) )                                               { priority = 2; sample_category =  "wwg_incl"; }
    else if ( dsname.Contains( "/TEST-www" ) )                                          { priority = 1; sample_category =  "www"; }
    else if ( dsname.Contains( "/TEST-tth" ) )                                          { priority = 2; sample_category =  "tthmia"; }
    else if ( dsname.Contains( "Run2016" ) && dsname.Contains( "DoubleMuon" ) )         { priority = 1; sample_category =  "data_mm"; }
    else if ( dsname.Contains( "Run2016" ) && dsname.Contains( "DoubleEG" ) )           { priority = 1; sample_category =  "data_ee"; }
    else if ( dsname.Contains( "Run2016" ) && dsname.Contains( "MuonEG" ) )             { priority = 1; sample_category =  "data_em"; }
    else if ( dsname.Contains( "Run2016" ) )                                            { priority = 1; sample_category =  "data_other"; }
    else if ( wwwbaby.evt() == 1532413419 && wwwbaby.lumi() == 962 && wwwbaby.run() == 281797 ) { priority = 1; sample_category =  "data_ee"; }
    sample_priority = priority;
    return sample_category;

    if ( dsname.Length() == 0 )
    {
        std::cout << std::endl;
        std::cout << "evt_dataset is empty! ";
        std::cout << wwwbaby.evt() << " " << wwwbaby.lumi() << " " << wwwbaby.run() << std::endl;
    }

    TString othername = dsname + "_UNCATEGORIZED";
    sample_priority = 2;
    sample_category = othername.ReplaceAll( "/", "_" );
    return sample_category;
}

//______________________________________________________________________________________
bool isBkgCategorySet()
{
    if ( bkg_category_set_to_eventid == wwwbaby.evt() &&
         bkg_category_set_to_run     == wwwbaby.run() &&
         bkg_category_set_to_lumi    == wwwbaby.lumi() )
        return true;
    else
        return false;
}

//______________________________________________________________________________________
TString bkgCategory()
{

    if ( isBkgCategorySet() ) return bkg_category;

    // Now set the event id of the event the objects are selected from
    bkg_category_set_to_eventid = wwwbaby.evt();
    bkg_category_set_to_run     = wwwbaby.run();
    bkg_category_set_to_lumi    = wwwbaby.lumi();

    if ( wwwbaby.isData() ) 
    {
        bkg_category = sampleCategory();
        return sampleCategory();
    }

    const TString& dsname = wwwbaby.evt_dataset()[0];
    if ( dsname.Contains( "/VHToNonbb" ) && passWHWWW() )
    {
        bkg_category = "whwww";
        return bkg_category;
    }
    else if ( dsname.Contains( "/TEST-www" ) )
    {
        bkg_category = "www";
        return bkg_category;
    }

    std::vector<Int_t> iSS = lepidx["TightLepton"];
    std::vector<Int_t> iaSS = lepidx["LbntLepton"];
    std::vector<Int_t> ilSS = lepidx["LooseLepton"];

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
//        std::cout << iSS.size() << " " << iaSS.size() << " " << ilSS.size() << std::endl;

//        std::cout << std::endl;
//        std::cout << "bkgCategory() Why am I here?";
//        std::cout << wwwbaby.evt() << " " << wwwbaby.lumi() << " " << wwwbaby.run() << std::endl;
        bkg_category = "others";
        return bkg_category;
    }

    bkg_category = sn;

    return bkg_category;

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
bool isFakeSample()
{
    TString sample = sampleCategory();
    if ( sample.Contains( "tt1l" ) ) return true;
    if ( sample.Contains( "singletop" ) ) return true;
    if ( sample.Contains( "wj" ) ) return true;
    if ( sample.Contains( "dy" ) ) return true;
    return false;
}

//______________________________________________________________________________________
bool isPhotonConversionSample()
{
    TString sample = sampleCategory();
    if ( sample.EqualTo( "ttg" ) ) return true;
    if ( sample.EqualTo( "wg" ) ) return true;
    if ( sample.EqualTo( "zg" ) ) return true;
    return false;
}

//______________________________________________________________________________________
bool isRealLeptonEvent()
{
    if ( bkgCategory().Contains( "SSLL" ) ) return true;
    if ( bkgCategory().Contains( "3lLL" ) ) return true;
    if ( bkgCategory().Contains( "trueSS" ) ) return true;
    if ( bkgCategory().Contains( "true3L" ) ) return true;
    return false;
}

//______________________________________________________________________________________
bool isPhotonConversionEvent()
{
    if ( bkgCategory().Contains( "photonfakes" ) ) return true;
    if ( bkgCategory().Contains( "photondoublefakes" ) ) return true;
    return false;
}

//______________________________________________________________________________________
void printEventID()
{
    std::cout << std::endl;
    std::cout << "(evt, run, lumi) = (" << wwwbaby.evt() << ", " << wwwbaby.run() << ", " << wwwbaby.lumi() << ")" << std::endl;
    std::cout << wwwbaby.run() << ":" << wwwbaby.lumi() << ":" << wwwbaby.evt() << std::endl;
}

//______________________________________________________________________________________
double fakerate( int idx, int syst )
{
    if ( abs( wwwbaby.lep_pdgId()[idx] ) == 11 )
    {
        float relIso = wwwbaby.lep_relIso03EAv2()[idx];
        float coneptcorr = 0;
        if ( abs( wwwbaby.lep_etaSC()[idx] ) <= 1.479 )
            coneptcorr = std::max( 0., relIso - 0.0588 );
        else
            coneptcorr = std::max( 0., relIso - 0.0571 );
        float corrpt = wwwbaby.lep_p4()[idx].pt() * ( 1 + coneptcorr );
        float abseta = fabs( wwwbaby.lep_p4()[idx].eta() );

        // fSumw[4][3]=0.170493 , x=85 , y=1.9895 , error=0.0370097
        // fSumw[4][2]=0.101855 , x=85 , y=1.1395 , error=0.0601733
        // fSumw[4][1]=0.199399 , x=85 , y=0.4    , error=0.0770095

        // fSumw[3][3]=0.1736   , x=40 , y=1.9895 , error=0.014733
        // fSumw[3][2]=0.164039 , x=40 , y=1.1395 , error=0.0213855
        // fSumw[3][1]=0.169051 , x=40 , y=0.4    , error=0.0207381

        // fSumw[2][3]=0.217486 , x=25 , y=1.9895 , error=0.0248149
        // fSumw[2][2]=0.191887 , x=25 , y=1.1395 , error=0.0232803
        // fSumw[2][1]=0.195776 , x=25 , y=0.4    , error=0.0282233

        if ( corrpt >= 50. )
        {
            if      ( abseta >= 1.479 ) return 0.170493 + syst * 0.0370097;
            else if ( abseta >= 0.8 )   return 0.101855 + syst * 0.0601733;
            else                        return 0.199399 + syst * 0.0770095;
        }
        else if ( corrpt >= 30. )
        {
            if      ( abseta >= 1.479 ) return 0.1736   + syst * 0.014733 ;
            else if ( abseta >= 0.8 )   return 0.164039 + syst * 0.0213855;
            else                        return 0.169051 + syst * 0.0207381;
        }
        else if ( corrpt >= 20. )
        {
            if      ( abseta >= 1.479 ) return 0.217486 + syst * 0.0248149;
            else if ( abseta >= 0.8 )   return 0.191887 + syst * 0.0232803;
            else                        return 0.195776 + syst * 0.0282233;
        }
        else
            return 0;
    }
    else if ( abs( wwwbaby.lep_pdgId()[idx] ) == 13 )
    {
        float relIso = wwwbaby.lep_relIso03EAv2()[idx];
        float coneptcorr = std::max( 0., relIso - 0.06 );
        float corrpt = wwwbaby.lep_p4()[idx].pt() * ( 1 + coneptcorr );
        float abseta = fabs( wwwbaby.lep_p4()[idx].eta() );

        // Data with fabs(eta) not applied by accident
        // fSumw[4][3]=0.0229585 , x=85 , y=2.25 , error=0.0417386
        // fSumw[4][2]=0.0229585 , x=85 , y=1.65 , error=0.0417386
        // fSumw[4][1]=0.0760604 , x=85 , y=0.6  , error=0.0334087

        // fSumw[3][3]=0.105424  , x=40 , y=2.25 , error=0.0156435
        // fSumw[3][2]=0.0727306 , x=40 , y=1.65 , error=0.00683751
        // fSumw[3][1]=0.0459721 , x=40 , y=0.6  , error=0.00482489

        // fSumw[2][3]=0.12044   , x=25 , y=2.25 , error=0.0282658
        // fSumw[2][2]=0.0866916 , x=25 , y=1.65 , error=0.0093254
        // fSumw[2][1]=0.054216  , x=25 , y=0.6  , error=0.00640781

        // QCD MC with correct abseta
        // fSumw[4][3]=0.0583592, x=85, y=2.25, error=0.0227116
        // fSumw[4][2]=0.0554221, x=85, y=1.65, error=0.0171941
        // fSumw[4][1]=0.0317323, x=85, y=0.6, error=0.00831419

        // fSumw[3][3]=0.0674571, x=40, y=2.25, error=0.0151419
        // fSumw[3][2]=0.0673591, x=40, y=1.65, error=0.00792283
        // fSumw[3][1]=0.045541, x=40, y=0.6, error=0.00546544

        // fSumw[2][3]=0.124618, x=25, y=2.25, error=0.0126553
        // fSumw[2][2]=0.109371, x=25, y=1.65, error=0.00733365
        // fSumw[2][1]=0.0592033, x=25, y=0.6, error=0.00397246

//        if ( corrpt >= 50. )
//        {
//            if      ( abseta >= 2.1 ) return 0.0229585 + syst * 0.0417386;
//            else if ( abseta >= 1.2 ) return 0.0229585 + syst * 0.0417386;
//            else                      return 0.0760604 + syst * 0.0334087;
//        }
//        else if ( corrpt >= 30. )
//        {
//            if      ( abseta >= 2.1 ) return 0.105424  + syst * 0.0156435 ;
//            else if ( abseta >= 1.2 ) return 0.0727306 + syst * 0.00683751;
//            else                      return 0.0459721 + syst * 0.00482489;
//        }
//        else if ( corrpt >= 20. )
//        {
//            if      ( abseta >= 2.1 ) return 0.12044   + syst * 0.0282658 ;
//            else if ( abseta >= 1.2 ) return 0.0866916 + syst * 0.0093254 ;
//            else                      return 0.054216  + syst * 0.00640781;
//        }
//        else
//            return 0;


        if ( corrpt >= 50. )
        {
            if      ( abseta >= 2.1 ) return 0.0583592 + syst * 0.0227116 ;
            else if ( abseta >= 1.2 ) return 0.0554221 + syst * 0.0171941 ;
            else                      return 0.0317323 + syst * 0.00831419;
        }
        else if ( corrpt >= 30. )
        {
            if      ( abseta >= 2.1 ) return 0.0674571 + syst * 0.0151419 ;
            else if ( abseta >= 1.2 ) return 0.0673591 + syst * 0.00792283;
            else                      return 0.045541  + syst * 0.00546544;
        }
        else if ( corrpt >= 20. )
        {
            if      ( abseta >= 2.1 ) return 0.124618  + syst * 0.0126553 ;
            else if ( abseta >= 1.2 ) return 0.109371  + syst * 0.00733365;
            else                      return 0.0592033 + syst * 0.00397246;
        }
        else
            return 0;

        // QCD
        // fSumw[4][3]=0.0583592, x=85, y=2.25, error=0.0227116
        // fSumw[4][2]=0.0554221, x=85, y=1.65, error=0.0171941
        // fSumw[4][1]=0.0317323, x=85, y=0.6 , error=0.00831419

        // fSumw[3][3]=0.0674571, x=40, y=2.25, error=0.0151419
        // fSumw[3][2]=0.0673591, x=40, y=1.65, error=0.00792283
        // fSumw[3][1]=0.045541 , x=40, y=0.6 , error=0.00546544

        // fSumw[2][3]=0.124618 , x=25, y=2.25, error=0.0126553
        // fSumw[2][2]=0.109371 , x=25, y=1.65, error=0.00733365
        // fSumw[2][1]=0.0592033, x=25, y=0.6 , error=0.00397246

// fSumw[4][3]=0.0583592, x=85, y=2.25, error=0.0227116
// fSumw[4][2]=0.0554221, x=85, y=1.65, error=0.0171941
// fSumw[4][1]=0.0317323, x=85, y=0.6, error=0.00831419

// fSumw[3][3]=0.0674571, x=40, y=2.25, error=0.0151419
// fSumw[3][2]=0.0673591, x=40, y=1.65, error=0.00792283
// fSumw[3][1]=0.045541, x=40, y=0.6, error=0.00546544

// fSumw[2][3]=0.124618, x=25, y=2.25, error=0.0126553
// fSumw[2][2]=0.109371, x=25, y=1.65, error=0.00733365
// fSumw[2][1]=0.0592033, x=25, y=0.6, error=0.00397246

    }
    else
        return 0;
}

//______________________________________________________________________________________
double fakefactor( int idx, int syst )
{
    return fakerate( idx, syst ) / ( 1. - fakerate( idx, syst ) );
}

