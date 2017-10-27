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

//_________________________________________________________________________________________________
bool passSSPresel(TString lepid)
{
    setObjectIndices();
    if (!( passTrigMM() || passTrigEM() || passTrigEE()              )) return false;
    if (!( lepidx[lepid].size()                              ==   2  )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           == 169  )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                             )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt()           >   30. )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt()           >   30. )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           >    0  )) return false;
    if (!( jetidx["GoodSSJet"].size()                        >=   2  )) return false;
    if (!( lepidx["VetoLepton"].size()                       ==   2  )) return false;
    if (!( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0  )) return false;
    if (!( jetidx["LooseBJet"].size()                        ==   0  )) return false;
    return true;
}

//_________________________________________________________________________________________________
bool passSSMM( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
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

//_________________________________________________________________________________________________
bool passSSEM( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
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
    if (!( passFilters()                                                             )) return false;
    if (!( passTrigEE()                                                              )) return false;
    if (!( lepidx[lepid].size()                              ==   2                  )) return false;
    if (!( wwwbaby.lep_pdgId()[lepidx[lepid][0]]
           * wwwbaby.lep_pdgId()[lepidx[lepid][1]]           == 121                  )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                                             )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0                  )) return false;
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
    if (!( fabs(
             ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass() - MZ ) >   10.            )) return false;
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

//______________________________________________________________________________________
bool passWZCR1SFOS()
{
    int numsfos = getNumSFOS();

    if (numsfos == 1)
    {
        if (!( pass3L1SFOS( "TightLepton", true, true )  )) return false;
        if (!( get1SFOSMll() > 55. && get1SFOSMll() < 110. )) return false;
        return true;
    }
    else
        return false;
}

//______________________________________________________________________________________
bool passWZCR2SFOS()
{
    int numsfos = getNumSFOS();

    if (numsfos == 2)
    {
        if (!( pass3L2SFOS( "TightLepton", true, true ) )) return false;
        if (!( ( fabs( get2SFOSMll0() - MZ ) < 20. ) ||
               ( fabs( get2SFOSMll1() - MZ ) < 20. )      )) return false;
        return true;
    }
    else
        return false;
}

//_________________________________________________________________________________________________
bool passSSWZCRMM( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
    if (!( passTrigMM()                                                             )) return false;
    if (!( lepidx[lepid].size()                              ==   3                 )) return false;
    if (!( LepFlavProductOfSSfor3L()                         == 169                 )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                                            )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt()           >   20.                )) return false;
    if (!( jetidx["GoodSSJet"].size()                        >=   2                 )) return false;
    if (!( lepidx["VetoLepton"].size()                       ==   3                 )) return false;
    if (!( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0                 )) return false;
    if (!( jetidx["LooseBJet"].size()                        ==   0     || dropbtag )) return false;
    if (!( MjjLead()                                         <  400.    || dropjet  )) return false;
    if (!( DEtajjLead()                                      <    1.5   || dropjet  )) return false;
    if (!( ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   >   40.                )) return false;
    if (getNumSFOS() == 1)
    {
        if (!( fabs( get1SFOSMll() - MZ ) < 10.0 )) return false;
    }
    else if (getNumSFOS() == 2)
    {
        if (!( ( fabs( get2SFOSMll0() - MZ ) < 10. ) ||
               ( fabs( get2SFOSMll1() - MZ ) < 10. )      )) return false;
    }
    else
        return false;
    return true;
}

//_________________________________________________________________________________________________
bool passSSWZCREM( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
    if (!( passTrigEM()                                                             )) return false;
    if (!( lepidx[lepid].size()                              ==   3                 )) return false;
    if (!( LepFlavProductOfSSfor3L()                         == 143                 )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                                            )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt()           >   30.                )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt()           >   20.                )) return false;
    if (!( jetidx["GoodSSJet"].size()                        >=   2                 )) return false;
    if (!( lepidx["VetoLepton"].size()                       ==   3                 )) return false;
    if (!( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0                 )) return false;
    if (!( jetidx["LooseBJet"].size()                        ==   0     || dropbtag )) return false;
    if (!( MjjLead()                                         <  400.    || dropjet  )) return false;
    if (!( DEtajjLead()                                      <    1.5   || dropjet  )) return false;
    if (!( ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   >   30.                )) return false;
    if (!( wwwbaby.met_pt()                                  >   40.                )) return false;
    if (!( MTmaxSS()                                         >   90.                )) return false;
    if (getNumSFOS() == 1)
    {
        if (!( fabs( get1SFOSMll() - MZ ) < 10.0 )) return false;
    }
    else if (getNumSFOS() == 2)
    {
        if (!( ( fabs( get2SFOSMll0() - MZ ) < 10. ) ||
               ( fabs( get2SFOSMll1() - MZ ) < 10. )      )) return false;
    }
    else
        return false;
    return true;
}

//______________________________________________________________________________________
bool passSSWZCREE( TString lepid, bool dropbtag, bool dropmjj, bool dropjet )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
    if (!( passTrigEE()                                                              )) return false;
    if (!( lepidx[lepid].size()                              ==   3                  )) return false;
    if (!( LepFlavProductOfSSfor3L()                         == 121                  )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                                             )) return false;
    if (!( wwwbaby.firstgoodvertex()                         ==   0                  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt()           >   30.                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt()           >   30.                 )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt()           >   20.                 )) return false;
    if (!( jetidx["GoodSSJet"].size()                        >=   2                  )) return false;
    if (!( lepidx["VetoLepton"].size()                       ==   3                  )) return false;
    if (!( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==   0                  )) return false;
    if (!( jetidx["LooseBJet"].size()                        ==   0      || dropbtag )) return false;
    if (!( MjjLead()                                         <  400.     || dropjet  )) return false;
    if (!( DEtajjLead()                                      <    1.5    || dropjet  )) return false;
    if (!( ( wwwbaby.lep_p4()[lepidx[lepid][0]]
             + wwwbaby.lep_p4()[lepidx[lepid][1]] ).mass()   >   40.                 )) return false;
    if (!( wwwbaby.met_pt()                                  >   40.                 )) return false;
    if (!( fabs( MllSS() - MZ )                              >   10.                 )) return false;
    if (getNumSFOS() == 1)
    {
        if (!( fabs( get1SFOSMll() - MZ ) < 10.0 )) return false;
    }
    else if (getNumSFOS() == 2)
    {
        if (!( ( fabs( get2SFOSMll0() - MZ ) < 10. ) ||
               ( fabs( get2SFOSMll1() - MZ ) < 10. )      )) return false;
    }
    else
        return false;
    return true;
}

//______________________________________________________________________________________
bool passFilters()
{
    if (!( wwwbaby.isData()                                   )) return true;//don't apply filters on simulation
    if (!( wwwbaby.Flag_EcalDeadCellTriggerPrimitiveFilter()  )) return false;
    if (!( wwwbaby.Flag_badChargedCandidateFilter         ()  )) return false;
    if (!( wwwbaby.Flag_HBHENoiseFilter                   ()  )) return false;
    if (!( wwwbaby.Flag_HBHEIsoNoiseFilter                ()  )) return false;
    if (!( wwwbaby.Flag_goodVertices                      ()  )) return false;
    if (!( wwwbaby.Flag_eeBadScFilter                     ()  )) return false;
    if (!( wwwbaby.Flag_globalTightHalo2016               ()  )) return false;
    if (!( wwwbaby.Flag_badMuonFilter                     ()  )) return false;
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
bool pass3Lpresel( TString lepid )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                    )) return false;
    if (!( wwwbaby.firstgoodvertex()               ==   0   )) return false;
    if (!( lepidx[lepid].size()                    ==   3   )) return false;
    if (!( lepidx["VetoLepton"].size()             ==   3   )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt() >   25.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt() >   20.  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt() >   20.  )) return false;
    if (!( abs( totalCharge() )                    ==   1   )) return false;
    if (!( jetidx["Good3LJet"].size()              <=   1   )) return false;
    if (!( jetidx["LooseBJet"].size()              ==   0   )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3LARpresel()
{
    setObjectIndices();
    if (!( lepidx["TightLepton"].size() == 2 )) return false;
    if (!( pass3Lpresel( "LooseLepton" )        )) return false;
    return true;
}

//______________________________________________________________________________________
// The 3L0SFOS signal region definition
bool pass3L0SFOS( TString lepid )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
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
    if (!( fabs( M3l() - MZ )                      >   10.  )) return false;
    return true;
}

//______________________________________________________________________________________
// The 3L1SFOS signal region definition
bool pass3L1SFOS( TString lepid, bool dropZ, bool dropV )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                            )) return false;
    if (!( wwwbaby.firstgoodvertex()               ==   0           )) return false;
    if (!( lepidx[lepid].size()                    ==   3           )) return false;
    if (!( lepidx["VetoLepton"].size()             ==   3  || dropV )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt() >   25.          )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt() >   20.          )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt() >   20.          )) return false;
    if (!( abs( totalCharge() )                    ==   1           )) return false;
    if (!( jetidx["Good3LJet"].size()              <=   1           )) return false;
    if (!( Pt3l()                                  >   60.          )) return false;
    if (!( jetidx["LooseBJet"].size()              ==   0           )) return false;
    if (!( getNumSFOS()                            ==   1           )) return false;
    if (!( wwwbaby.met_pt()                        >   45.          )) return false;
    if (!( get1SFOSMll()                           >   20.          )) return false;
    if (!( get1SFOSMll()                           <   55. ||
           get1SFOSMll()                           >  100. || dropZ )) return false;
    if (!( DPhi3lMET()                             >    2.5         )) return false;
    if (!( fabs( M3l() - MZ )                      >   10.          )) return false;
    return true;
}

//______________________________________________________________________________________
// The 3L2SFOS signal region definition
bool pass3L2SFOS( TString lepid, bool dropZ, bool dropV )
{
    setObjectIndices();
    if (!( passFilters()                                                            )) return false;
    if (!( wwwbaby.evt_passgoodrunlist()                            )) return false;
    if (!( wwwbaby.firstgoodvertex()               ==   0           )) return false;
    if (!( lepidx[lepid].size()                    ==   3           )) return false;
    if (!( lepidx["VetoLepton"].size()             ==   3 || dropV  )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][0]].pt() >   25.          )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][1]].pt() >   20.          )) return false;
    if (!( wwwbaby.lep_p4()[lepidx[lepid][2]].pt() >   20.          )) return false;
    if (!( abs( totalCharge() )                    ==   1           )) return false;
    if (!( jetidx["Good3LJet"].size()              <=   1           )) return false;
    if (!( Pt3l()                                  >   60.          )) return false;
    if (!( jetidx["LooseBJet"].size()              ==   0           )) return false;
    if (!( getNumSFOS()                            ==   2           )) return false;
    if (!( wwwbaby.met_pt()                        >   55.          )) return false;
    if (!(       get2SFOSMll0()                    >   20.          )) return false;
    if (!(       get2SFOSMll1()                    >   20.          )) return false;
    if (!( fabs( get2SFOSMll0() - MZ )             >   20. || dropZ )) return false;
    if (!( fabs( get2SFOSMll1() - MZ )             >   20. || dropZ )) return false;
    if (!( DPhi3lMET()                             >    2.5         )) return false;
    if (!( fabs( M3l() - MZ )                      >   10.          )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3LAR0SFOS()
{
    if (!( lepidx["TightLepton"].size() == 2 )) return false;
    if (!( pass3L0SFOS( "LooseLepton" )        )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3LAR1SFOS()
{
    if (!( lepidx["TightLepton"].size() == 2 )) return false;
    if (!( pass3L1SFOS( "LooseLepton" )        )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3LAR2SFOS()
{
    if (!( lepidx["TightLepton"].size() == 2 )) return false;
    if (!( pass3L2SFOS( "LooseLepton" )        )) return false;
    return true;
}

//______________________________________________________________________________________
//bool passTrigMM() { return wwwbaby.HLT_DoubleMu() || wwwbaby.HLT_DoubleMu_noiso(); }
//bool passTrigEM() { return wwwbaby.HLT_MuEG() || wwwbaby.HLT_MuEG_noiso(); }
//bool passTrigEE() { return wwwbaby.HLT_DoubleEl() || wwwbaby.HLT_DoubleEl_DZ() || wwwbaby.HLT_DoubleEl_DZ_2() || wwwbaby.HLT_DoubleEl_noiso(); }
bool passTrigMM() { if (wwwbaby.isData()) return wwwbaby.HLT_DoubleMu(); else return true; }
bool passTrigEM() { if (wwwbaby.isData()) return wwwbaby.HLT_MuEG(); else return true; }
bool passTrigEE() { if (wwwbaby.isData()) return wwwbaby.HLT_DoubleEl() || wwwbaby.HLT_DoubleEl_DZ() || wwwbaby.HLT_DoubleEl_DZ_2(); else return true; }

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
        else if ( lepidx["TightLepton"].size() == 2 && lepidx["LooseLepton"].size() == 3 )
            wgt *= fakefactor( lepidx["LbntLepton"][0], ffsyst );
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
    std::cout  << "M3l           " << " : " << M3l            () << std::endl;
    if (getNumSFOS() == 0)
    {
        std::cout  << "get0SFOSMll   " << " : " << get0SFOSMll    () << std::endl;
        std::cout  << "get0SFOSMee   " << " : " << get0SFOSMee    () << std::endl;
    }
    else if (getNumSFOS() == 1)
    {
        std::cout  << "get1SFOSMll   " << " : " << get1SFOSMll    () << std::endl;
    }
    else if (getNumSFOS() == 2)
    {
        std::cout  << "get2SFOSMll0  " << " : " << get2SFOSMll0   () << std::endl;
        std::cout  << "get2SFOSMll1  " << " : " << get2SFOSMll1   () << std::endl;
    }
    std::cout  << "met           " << " : " << wwwbaby.met_pt () << std::endl;
    std::cout  << "ntrkiso       " << " : " << wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() << std::endl;
    std::cout  << "grl           " << " : " << wwwbaby.evt_passgoodrunlist() << std::endl;
    std::cout  << "firstvertex   " << " : " << wwwbaby.firstgoodvertex() << std::endl;
    std::cout  << "passTrigEE    " << " : " << passTrigEE() << std::endl;
    std::cout  << "passTrigEM    " << " : " << passTrigEM() << std::endl;
    std::cout  << "passTrigMM    " << " : " << passTrigMM() << std::endl;
    std::cout  << "LepFlavProductOfSSfor3L" << " : " << LepFlavProductOfSSfor3L() << std::endl;
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
    idx["TightLepton"].clear();
    idx["LooseLepton"].clear();
    idx["VetoLepton"].clear();
    idx["LbntLepton"].clear();
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
    for ( unsigned int ilep = 0; ilep < wwwbaby.lep_pdgId().size(); ++ilep )
    {
//        if ( isTightLepton( ilep ) ) idx["TightLepton"].push_back( ilep );
//        if ( isLooseLepton( ilep ) ) idx["LooseLepton"].push_back( ilep );
//        if ( isVetoLepton( ilep ) )  idx["VetoLepton"].push_back( ilep );
//        if ( isLbntLepton( ilep ) )  idx["LbntLepton"].push_back( ilep );
        if ( isTightLepton( ilep ) )
        {
            idx["TightLepton"].push_back( ilep );
            idx["LooseLepton"].push_back( ilep );
            idx["VetoLepton"].push_back( ilep );
        }
        else if ( isLooseLepton( ilep ) )
        {
            idx["LbntLepton"].push_back( ilep );
            idx["LooseLepton"].push_back( ilep );
            idx["VetoLepton"].push_back( ilep );
        }
        else if ( isVetoLepton( ilep ) )
        {
            idx["VetoLepton"].push_back( ilep );
        }
    }
    return idx;
}

//______________________________________________________________________________________
void loadLeptonIndices()
{
    lepidx["TightLepton"]   = wwwbaby.TightLepton();
    lepidx["LooseLepton"]   = wwwbaby.LooseLepton();
    lepidx["VetoLepton"]    = wwwbaby.VetoLepton();
    lepidx["LbntLepton"]    = wwwbaby.LbntLepton();
}

//______________________________________________________________________________________
ObjIdx getJetsIndices()
{
    ObjIdx idx = makeEmptyJetidx();
    for ( unsigned int ijet = 0; ijet < wwwbaby.jets_csv().size(); ++ijet )
    {
        if ( isGoodSSJet( ijet ) )
        {
            idx["GoodSSJet"].push_back( ijet );
            idx["Good3LJet"].push_back( ijet );
        }
        else if ( isGood3LJet( ijet ) )
        {
            idx["Good3LJet"].push_back( ijet );
        }
        if ( isLooseBJet( ijet ) )
            idx["LooseBJet"].push_back( ijet );
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
    if (!(       wwwbaby.lep_isTriggerSafe_v1()[ilep]        )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0     )) return false;
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
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0     )) return false;
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
bool isGoodSSJet( unsigned int& ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].Pt()      > 30.  )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].Eta() )   <  2.5 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isGood3LJet( unsigned int& ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].Pt()      > 30.  )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].Eta() )   <  5.0 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseBJet( unsigned int& ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].Pt()      > 20.     )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].Eta() )   <  2.4    )) return false;
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
float MET()
{
    return wwwbaby.met_pt();
}

//______________________________________________________________________________________
float MllSS()
{
    if ( lepidx["VetoLepton"].size() < 3 )
        return -999;
    std::pair<int, int> sspair = getSSPair();
    return ( wwwbaby.lep_p4()[lepidx["VetoLepton"][sspair.first]]
             + wwwbaby.lep_p4()[lepidx["VetoLepton"][sspair.second]] ).mass();
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
    if ( lepidx["VetoLepton"].size() < 2 )
        return -999;
    return fabs( ROOT::Math::VectorUtil::DeltaPhi(
            wwwbaby.lep_p4()[lepidx["VetoLepton"][0]],
            wwwbaby.lep_p4()[lepidx["VetoLepton"][1]] ));
}

//______________________________________________________________________________________
float DEtall()
{
    if ( lepidx["VetoLepton"].size() < 2 )
        return -999;
    return fabs( wwwbaby.lep_p4()[lepidx["VetoLepton"][0]].eta()
             - wwwbaby.lep_p4()[lepidx["VetoLepton"][1]].eta() );
}

//______________________________________________________________________________________
float Pt3l()
{
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
             + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).pt();
}

//______________________________________________________________________________________
float M3l()
{
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
             + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).M();
}

//______________________________________________________________________________________
float DPhi3lMET()
{
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    TLorentzVector met;
    met.SetPtEtaPhiM( wwwbaby.met_pt(), 0, wwwbaby.met_phi(), 0 );
    LorentzVector metlv;
    metlv.SetPxPyPzE( met.Px(), met.Py(), met.Pz(), met.E() );
    return fabs( ROOT::Math::VectorUtil::DeltaPhi(
                wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
                + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]],
                met ));
}

//______________________________________________________________________________________
float MTmaxSS()
{
    if ( lepidx["VetoLepton"].size() < 3 )
        return -999;
    std::pair<int, int> sspair = getSSPair();
    float mt0 = MTidx( sspair.first, "VetoLepton" );
    float mt1 = MTidx( sspair.second, "VetoLepton" );
    return mt0 > mt1 ? mt0 : mt1;
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
std::pair<int, int> getSSPair()
{
    if ( lepidx["VetoLepton"].size() < 3 )
        return std::pair<int, int>(-1, -1);
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["VetoLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["VetoLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["VetoLepton"][2]];
    if ( ((pdgid0>0) - (pdgid0<0)) == ((pdgid1>0) - (pdgid1<0)) )
        return std::pair<int, int>(0, 1);
    else if ( ((pdgid0>0)-(pdgid0<0)) == ((pdgid2>0)-(pdgid2<0)) )
        return std::pair<int, int>(0, 2);
    else if ( ((pdgid1>0)-(pdgid1<0)) == ((pdgid2>0)-(pdgid2<0)) )
        return std::pair<int, int>(1, 2);
    else
        return std::pair<int, int>(-1, -1);
}

//______________________________________________________________________________________
int LepFlavProductOfSSfor3L()
{
    if ( lepidx["VetoLepton"].size() < 3 )
        return 0;
    std::pair<int, int> sspair = getSSPair();
    return wwwbaby.lep_pdgId()[lepidx["VetoLepton"][sspair.first]]
           * wwwbaby.lep_pdgId()[lepidx["VetoLepton"][sspair.second]];
}

//______________________________________________________________________________________
int LepFlavProduct()
{
    if ( lepidx["VetoLepton"].size() < 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["VetoLepton"][0]]
           * wwwbaby.lep_pdgId()[lepidx["VetoLepton"][1]];
}

//______________________________________________________________________________________
bool isLeadPlus()
{
    if ( lepidx["VetoLepton"].size() < 1 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["VetoLepton"][0]] < 0;
}

//______________________________________________________________________________________
bool isSubleadPlus()
{
    if ( lepidx["VetoLepton"].size() < 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["VetoLepton"][1]] < 0;
}

//______________________________________________________________________________________
bool isLeadMinus()
{
    if ( lepidx["VetoLepton"].size() < 1 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["VetoLepton"][0]] > 0;
}

//______________________________________________________________________________________
bool isSubleadMinus()
{
    if ( lepidx["VetoLepton"].size() < 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["VetoLepton"][1]] > 0;
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
    for ( auto& ilep : lepidx["LooseLepton"] )
        charge += wwwbaby.lep_charge()[ilep];
    return charge;
}

//______________________________________________________________________________________
int getNumSFOS()
{
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    /* Loops through pairs of entries in the lep_pdgId vector
       and counts how many have opposite value*/
    int num_SFOS = 0;
    for ( int i = 0; i < ( int ) lepidx["LooseLepton"].size(); i++ )
    {
        for ( int j = i + 1; j < ( int ) lepidx["LooseLepton"].size(); j++ )
        {
            if ( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][i]] ==
                    -( wwwbaby.lep_pdgId()[lepidx["LooseLepton"][j]] ) )
                num_SFOS++;
        }
    }
    return num_SFOS;
}

//______________________________________________________________________________________
float get0SFOSMll()
{
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][2]];
    if ( pdgid0 == pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
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
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][2]];
    if ( pdgid0 == pdgid1 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 && abs( pdgid1 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
    return -999;
}

//______________________________________________________________________________________
float get1SFOSMll()
{
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
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
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
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
    if ( lepidx["LooseLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["LooseLepton"][2]];
    if ( pdgid2 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid0 )
        return ( wwwbaby.lep_p4()[lepidx["LooseLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["LooseLepton"][1]] ).mass();
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
    if ( lepidx["VetoLepton"].size() < 2 )
        return -999;
    if ( jetidx["GoodSSJet"].size() < 2 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx["VetoLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["VetoLepton"][1]]
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
    return bkgCategory_v2();
}

//______________________________________________________________________________________
TString bkgCategory_v1()
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
TString bkgCategory_v2()
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

    if (lepidx["LooseLepton"].size() == 2) //have two SS leptons
    {
        int gentype = gentype_v2(lepidx["LooseLepton"][0], lepidx["LooseLepton"][1], -1);
        TString sn;
        if      (gentype == 0) { sn = "trueSS"; }
        else if (gentype == 2) { sn = "chargeflips"; }
        else if (gentype == 3) { sn = "SSLL"; }
        else if (gentype == 4) { sn = "fakes"; }
        else if (gentype == 5) { sn = "photonfakes"; }
        else                   { sn = "others"; }
        bkg_category = sn;
        return sn;
    }
    else if (lepidx["LooseLepton"].size() >= 3) //have three leptons for 3l channel
    {
        int gentype = gentype_v2(lepidx["LooseLepton"][0],lepidx["LooseLepton"][1], lepidx["LooseLepton"][2]);
        TString sn2;
        if      (gentype == 0) { sn2 = "trueWWW"; }
        else if (gentype == 1) { sn2 = "true3L"; } //in the end will be taken from data-driven estimate
        else if (gentype == 2) { sn2 = "chargeflips"; }
        else if (gentype == 3) { sn2 = "3lLL"; } //in the end will be taken from data-driven estimate
        else if (gentype == 4) { sn2 = "fakes"; } //in the end will be taken from data-driven estimate
        else if (gentype == 5) { sn2 = "photonfakes"; }
        else                   { sn2 = "others"; }//this gets never filled - is just a failsafe category
        bkg_category = sn2;
        return sn2;
    }
    else
        return "others";
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

//______________________________________________________________________________________
int gentype_v2(unsigned lep1_index, unsigned lep2_index, int lep3_index)
{
    bool gammafake = false;
    bool jetfake   = false;
    unsigned int ngenlep = wwwbaby.ngenLepFromTau() + wwwbaby.ngenLep();
    unsigned int nW(0), nZ(0);
    bool lep1_real = wwwbaby.lep_motherIdSS().at(lep1_index) > 0;
    bool lep2_real = wwwbaby.lep_motherIdSS().at(lep2_index) > 0;
    bool lep3_real = false;
    if (lep3_index > 0) { lep3_real = wwwbaby.lep_motherIdSS().at(lep3_index) > 0; }
    vector<int> reallepindex;
    for (unsigned int lepindex = 0; lepindex < wwwbaby.lep_p4().size(); ++lepindex)
    {
        if (wwwbaby.lep_motherIdSS().at(lepindex) > 0) { reallepindex.push_back(lepindex); }
        else if (wwwbaby.lep_motherIdSS().at(lepindex) == -3) { gammafake = true; }
        else                                           { jetfake = true; }
        if (wwwbaby.lep_isFromW().at(lepindex)) { nW++; }
        if (wwwbaby.lep_isFromZ().at(lepindex)) { nZ++; }
    }
    //found two real leptons
    if (lep3_index < 0)
    {
        bool ischargeflip = false;
        bool isSS = false;
        if (lep1_real && lep2_real)
        {
            int ilep1 =   wwwbaby.lep_genPart_index().at(lep1_index);
            int ilep2 =   wwwbaby.lep_genPart_index().at(lep2_index);
            bool lep1_chargeflip  = wwwbaby.genPart_charge().at(ilep1) != wwwbaby.lep_charge().at(lep1_index);
            bool lep2_chargeflip  = wwwbaby.genPart_charge().at(ilep2) != wwwbaby.lep_charge().at(lep2_index);
            if (!lep1_chargeflip && !lep2_chargeflip && nW == 2) { return 0; } // true SS
            else if (!lep1_chargeflip && !lep2_chargeflip) { isSS = true; } // true SS - but could be still lost lepton WZ
            if (lep1_chargeflip || lep2_chargeflip)   { ischargeflip = true; } // charge flip
        }
        if (ngenlep > 2 || reallepindex.size() > 2 || (nW > 0 && nZ > 0)) { return 3; } // lostlep
        if ((ngenlep < 2 || !lep1_real || !lep2_real) &&    jetfake) { return 4; } // jetfake - if double fake with one jet fake and one gamma fake call it jet fake
        if ((ngenlep < 2 || !lep1_real || !lep2_real) &&  gammafake) { return 5; } // gammafake
        if ((ngenlep < 2 || !lep1_real || !lep2_real) && !gammafake) { return 4; } // call all without gamma fake jetfake - safety cut
        if (isSS) { return 0; }
        if (ischargeflip) { return 2; }
        cout << "This event was not classified - 2 lepton event - v2" << endl;
        return 1;
    }
    else
    {
        //found three real leptons
        bool ischargeflip = false;
        bool isthreelep = false;
        if (lep1_real && lep2_real && lep3_real)
        {
            int ilep1 =   wwwbaby.lep_genPart_index().at(lep1_index);
            int ilep2 =   wwwbaby.lep_genPart_index().at(lep2_index);
            int ilep3 =   wwwbaby.lep_genPart_index().at(lep3_index);
            bool lep1_chargeflip  = wwwbaby.genPart_charge().at(ilep1) != wwwbaby.lep_charge().at(lep1_index);
            bool lep2_chargeflip  = wwwbaby.genPart_charge().at(ilep2) != wwwbaby.lep_charge().at(lep2_index);
            bool lep3_chargeflip  = wwwbaby.genPart_charge().at(ilep3) != wwwbaby.lep_charge().at(lep3_index);
            if (!lep1_chargeflip && !lep2_chargeflip && !lep3_chargeflip && nW == 3) { return 0; } // true WWW
            else if (!lep1_chargeflip && !lep2_chargeflip && !lep3_chargeflip) { isthreelep = true; } // true 3l, but could be lost lepton ZZ
            if (lep1_chargeflip || lep2_chargeflip || lep3_chargeflip)   { ischargeflip = true; } // charge flip
        }
        if (ngenlep > 3 || reallepindex.size() > 3 || (nW >= 2 && nZ >= 1) || (nZ >= 3)) { return 3; } // lostlep (2 lep from W and 2 from Z, or 4 from Z)
        //there is the case of having WZZ with two lost leptons --> ngenlep>3 - correctly put has lostlep
        if ((ngenlep < 3 || !lep1_real || !lep2_real || !lep3_real) &&    jetfake) { return 4; } // jetfake
        if ((ngenlep < 3 || !lep1_real || !lep2_real || !lep3_real) &&  gammafake) { return 5; } // gammafake
        if ((ngenlep < 3 || !lep1_real || !lep2_real || !lep3_real) && !gammafake) { return 4; } // jetfake
        if (isthreelep) { return 1; }
        if (ischargeflip) { return 2; }
        cout << "This event was not classified - 3 lepton event - v2" << endl;
        return 0;
    }
}

//______________________________________________________________________________________
float lepsf_EGammaTightPOG_EGammaVVV(float pt, float eta)
{
    if (pt < 30.0 && fabs(eta) < 1.0) return 0.980376899242;
    if (pt < 40.0 && fabs(eta) < 1.0) return 0.981454670429;
    if (pt < 50.0 && fabs(eta) < 1.0) return 0.980952382088;
    if (pt < 2000.0 && fabs(eta) < 1.0) return 0.976341068745;
    if (pt < 30.0 && fabs(eta) < 1.4442) return 0.945449829102;
    if (pt < 40.0 && fabs(eta) < 1.4442) return 0.94911891222;
    if (pt < 50.0 && fabs(eta) < 1.4442) return 0.95356208086;
    if (pt < 2000.0 && fabs(eta) < 1.4442) return 0.949157774448;
    if (pt < 30.0 && fabs(eta) < 1.56) return 0.938438951969;
    if (pt < 40.0 && fabs(eta) < 1.56) return 0.941834270954;
    if (pt < 50.0 && fabs(eta) < 1.56) return 0.947787821293;
    if (pt < 2000.0 && fabs(eta) < 1.56) return 0.952665507793;
    if (pt < 30.0 && fabs(eta) < 2.0) return 0.963943779469;
    if (pt < 40.0 && fabs(eta) < 2.0) return 0.971829116344;
    if (pt < 50.0 && fabs(eta) < 2.0) return 0.975824117661;
    if (pt < 2000.0 && fabs(eta) < 2.0) return 0.980078279972;
    if (pt < 30.0 && fabs(eta) < 2.5) return 0.984203636646;
    if (pt < 40.0 && fabs(eta) < 2.5) return 0.998542606831;
    if (pt < 50.0 && fabs(eta) < 2.5) return 1.00013875961;
    if (pt < 2000.0 && fabs(eta) < 2.5) return 1.00487089157;
    return 1;
}

//______________________________________________________________________________________
float lepsf_EGammaVVV_EGammaVVVEle12(float pt, float eta)
{
    if (pt < 30.0 && fabs(eta) < 1.0) return 1.01243853569;
    if (pt < 40.0 && fabs(eta) < 1.0) return 1.01314127445;
    if (pt < 50.0 && fabs(eta) < 1.0) return 1.01340222359;
    if (pt < 2000.0 && fabs(eta) < 1.0) return 1.01799964905;
    if (pt < 30.0 && fabs(eta) < 1.4442) return 1.01877808571;
    if (pt < 40.0 && fabs(eta) < 1.4442) return 1.01980721951;
    if (pt < 50.0 && fabs(eta) < 1.4442) return 1.01803350449;
    if (pt < 2000.0 && fabs(eta) < 1.4442) return 1.01880669594;
    if (pt < 30.0 && fabs(eta) < 1.56) return 1.01981687546;
    if (pt < 40.0 && fabs(eta) < 1.56) return 1.01999163628;
    if (pt < 50.0 && fabs(eta) < 1.56) return 1.01699542999;
    if (pt < 2000.0 && fabs(eta) < 1.56) return 1.01749408245;
    if (pt < 30.0 && fabs(eta) < 2.0) return 1.02057921886;
    if (pt < 40.0 && fabs(eta) < 2.0) return 1.02028226852;
    if (pt < 50.0 && fabs(eta) < 2.0) return 1.01790201664;
    if (pt < 2000.0 && fabs(eta) < 2.0) return 1.02002465725;
    if (pt < 30.0 && fabs(eta) < 2.5) return 0.978350162506;
    if (pt < 40.0 && fabs(eta) < 2.5) return 1.00657308102;
    if (pt < 50.0 && fabs(eta) < 2.5) return 1.00991177559;
    if (pt < 2000.0 && fabs(eta) < 2.5) return 1.01343548298;
    return 1;
}

//______________________________________________________________________________________
float lepsf_EGammaVVV_EGammaVVVEleLead(float pt, float eta)
{
    if (pt < 30.0 && fabs(eta) < 1.0) return 0.973178207874;
    if (pt < 40.0 && fabs(eta) < 1.0) return 1.01313531399;
    if (pt < 50.0 && fabs(eta) < 1.0) return 1.01339864731;
    if (pt < 2000.0 && fabs(eta) < 1.0) return 1.01799964905;
    if (pt < 30.0 && fabs(eta) < 1.4442) return 0.969570457935;
    if (pt < 40.0 && fabs(eta) < 1.4442) return 1.01976215839;
    if (pt < 50.0 && fabs(eta) < 1.4442) return 1.01803350449;
    if (pt < 2000.0 && fabs(eta) < 1.4442) return 1.01880669594;
    if (pt < 30.0 && fabs(eta) < 1.56) return 0.974450111389;
    if (pt < 40.0 && fabs(eta) < 1.56) return 1.01983559132;
    if (pt < 50.0 && fabs(eta) < 1.56) return 1.01696801186;
    if (pt < 2000.0 && fabs(eta) < 1.56) return 1.01749408245;
    if (pt < 30.0 && fabs(eta) < 2.0) return 0.980933725834;
    if (pt < 40.0 && fabs(eta) < 2.0) return 1.02049899101;
    if (pt < 50.0 && fabs(eta) < 2.0) return 1.01789569855;
    if (pt < 2000.0 && fabs(eta) < 2.0) return 1.01624071598;
    if (pt < 30.0 && fabs(eta) < 2.5) return 0.932764172554;
    if (pt < 40.0 && fabs(eta) < 2.5) return 1.00653338432;
    if (pt < 50.0 && fabs(eta) < 2.5) return 1.00987935066;
    if (pt < 2000.0 && fabs(eta) < 2.5) return 1.01343548298;
    return 1;
}

//______________________________________________________________________________________
float lepsf_MuMediumPOG_MuTightVVV(float pt, float eta)
{
    if (pt < 25.0 && fabs(eta) < 0.9) return 0.972987830639;
    if (pt < 30.0 && fabs(eta) < 0.9) return 0.980633199215;
    if (pt < 40.0 && fabs(eta) < 0.9) return 0.986259222031;
    if (pt < 50.0 && fabs(eta) < 0.9) return 0.989767193794;
    if (pt < 60.0 && fabs(eta) < 0.9) return 0.990144252777;
    if (pt < 120.0 && fabs(eta) < 0.9) return 0.992544472218;
    if (pt < 200.0 && fabs(eta) < 0.9) return 1.00607073307;
    if (pt < 1200.0 && fabs(eta) < 0.9) return 1.00708353519;
    if (pt < 25.0 && fabs(eta) < 1.2) return 0.982158899307;
    if (pt < 30.0 && fabs(eta) < 1.2) return 0.99262791872;
    if (pt < 40.0 && fabs(eta) < 1.2) return 0.994507312775;
    if (pt < 50.0 && fabs(eta) < 1.2) return 0.994486868382;
    if (pt < 60.0 && fabs(eta) < 1.2) return 0.995231330395;
    if (pt < 120.0 && fabs(eta) < 1.2) return 0.996162414551;
    if (pt < 200.0 && fabs(eta) < 1.2) return 0.998247504234;
    if (pt < 1200.0 && fabs(eta) < 1.2) return 1.00434410572;
    if (pt < 25.0 && fabs(eta) < 2.1) return 0.992021918297;
    if (pt < 30.0 && fabs(eta) < 2.1) return 0.999400854111;
    if (pt < 40.0 && fabs(eta) < 2.1) return 0.999857544899;
    if (pt < 50.0 && fabs(eta) < 2.1) return 0.999814629555;
    if (pt < 60.0 && fabs(eta) < 2.1) return 1.00330460072;
    if (pt < 120.0 && fabs(eta) < 2.1) return 1.00755393505;
    if (pt < 200.0 && fabs(eta) < 2.1) return 1.01337754726;
    if (pt < 1200.0 && fabs(eta) < 2.1) return 1.02070629597;
    if (pt < 25.0 && fabs(eta) < 2.4) return 0.969161331654;
    if (pt < 30.0 && fabs(eta) < 2.4) return 0.975010037422;
    if (pt < 40.0 && fabs(eta) < 2.4) return 0.990625023842;
    if (pt < 50.0 && fabs(eta) < 2.4) return 0.994753956795;
    if (pt < 60.0 && fabs(eta) < 2.4) return 0.992071866989;
    if (pt < 120.0 && fabs(eta) < 2.4) return 0.994121670723;
    if (pt < 200.0 && fabs(eta) < 2.4) return 0.994795560837;
    if (pt < 1200.0 && fabs(eta) < 2.4) return 1.01750481129;
    return 1;
}

//______________________________________________________________________________________
float lepsf_MuTightVVV_MuTightVVVMu17(float pt, float eta)
{
    if (pt < 25.0 && fabs(eta) < 0.9) return 0.990568935871;
    if (pt < 30.0 && fabs(eta) < 0.9) return 0.989429354668;
    if (pt < 40.0 && fabs(eta) < 0.9) return 0.989176571369;
    if (pt < 50.0 && fabs(eta) < 0.9) return 0.988689124584;
    if (pt < 60.0 && fabs(eta) < 0.9) return 0.988338470459;
    if (pt < 120.0 && fabs(eta) < 0.9) return 0.986210823059;
    if (pt < 200.0 && fabs(eta) < 0.9) return 0.986462116241;
    if (pt < 1200.0 && fabs(eta) < 0.9) return 0.992106556892;
    if (pt < 25.0 && fabs(eta) < 1.2) return 0.979321241379;
    if (pt < 30.0 && fabs(eta) < 1.2) return 0.979366958141;
    if (pt < 40.0 && fabs(eta) < 1.2) return 0.984407544136;
    if (pt < 50.0 && fabs(eta) < 1.2) return 0.98391777277;
    if (pt < 60.0 && fabs(eta) < 1.2) return 0.982034981251;
    if (pt < 120.0 && fabs(eta) < 1.2) return 0.981674432755;
    if (pt < 200.0 && fabs(eta) < 1.2) return 0.967723488808;
    if (pt < 1200.0 && fabs(eta) < 1.2) return 0.97421169281;
    if (pt < 25.0 && fabs(eta) < 2.1) return 0.988342583179;
    if (pt < 30.0 && fabs(eta) < 2.1) return 0.997145831585;
    if (pt < 40.0 && fabs(eta) < 2.1) return 0.997771382332;
    if (pt < 50.0 && fabs(eta) < 2.1) return 0.999460279942;
    if (pt < 60.0 && fabs(eta) < 2.1) return 0.996708571911;
    if (pt < 120.0 && fabs(eta) < 2.1) return 0.980419993401;
    if (pt < 200.0 && fabs(eta) < 2.1) return 0.921213328838;
    if (pt < 1200.0 && fabs(eta) < 2.1) return 0.875039458275;
    if (pt < 25.0 && fabs(eta) < 2.4) return 0.902593195438;
    if (pt < 30.0 && fabs(eta) < 2.4) return 0.927266776562;
    if (pt < 40.0 && fabs(eta) < 2.4) return 0.943869292736;
    if (pt < 50.0 && fabs(eta) < 2.4) return 0.958442807198;
    if (pt < 60.0 && fabs(eta) < 2.4) return 0.96275305748;
    if (pt < 120.0 && fabs(eta) < 2.4) return 0.966543674469;
    if (pt < 200.0 && fabs(eta) < 2.4) return 0.951920926571;
    if (pt < 1200.0 && fabs(eta) < 2.4) return 0.870018422604;
    return 1;
}

//______________________________________________________________________________________
float lepsf_MuTightVVV_MuTightVVVMu8(float pt, float eta)
{
    if (pt < 25.0 && fabs(eta) < 0.9) return 0.999541461468;
    if (pt < 30.0 && fabs(eta) < 0.9) return 0.997860491276;
    if (pt < 40.0 && fabs(eta) < 0.9) return 0.997504532337;
    if (pt < 50.0 && fabs(eta) < 0.9) return 0.996784865856;
    if (pt < 60.0 && fabs(eta) < 0.9) return 0.996682703495;
    if (pt < 120.0 && fabs(eta) < 0.9) return 0.994662582874;
    if (pt < 200.0 && fabs(eta) < 0.9) return 0.995876908302;
    if (pt < 1200.0 && fabs(eta) < 0.9) return 1.00186300278;
    if (pt < 25.0 && fabs(eta) < 1.2) return 1.00136160851;
    if (pt < 30.0 && fabs(eta) < 1.2) return 1.00054740906;
    if (pt < 40.0 && fabs(eta) < 1.2) return 1.00397896767;
    if (pt < 50.0 && fabs(eta) < 1.2) return 1.00303781033;
    if (pt < 60.0 && fabs(eta) < 1.2) return 1.00098395348;
    if (pt < 120.0 && fabs(eta) < 1.2) return 1.00059711933;
    if (pt < 200.0 && fabs(eta) < 1.2) return 0.994088232517;
    if (pt < 1200.0 && fabs(eta) < 1.2) return 0.998805344105;
    if (pt < 25.0 && fabs(eta) < 2.1) return 1.01416134834;
    if (pt < 30.0 && fabs(eta) < 2.1) return 1.01810801029;
    if (pt < 40.0 && fabs(eta) < 2.1) return 1.01771843433;
    if (pt < 50.0 && fabs(eta) < 2.1) return 1.01809024811;
    if (pt < 60.0 && fabs(eta) < 2.1) return 1.01902592182;
    if (pt < 120.0 && fabs(eta) < 2.1) return 1.01623928547;
    if (pt < 200.0 && fabs(eta) < 2.1) return 0.983312547207;
    if (pt < 1200.0 && fabs(eta) < 2.1) return 0.948991417885;
    if (pt < 25.0 && fabs(eta) < 2.4) return 1.00432026386;
    if (pt < 30.0 && fabs(eta) < 2.4) return 1.00704360008;
    if (pt < 40.0 && fabs(eta) < 2.4) return 1.00508832932;
    if (pt < 50.0 && fabs(eta) < 2.4) return 1.00839221478;
    if (pt < 60.0 && fabs(eta) < 2.4) return 1.00767874718;
    if (pt < 120.0 && fabs(eta) < 2.4) return 1.00606381893;
    if (pt < 200.0 && fabs(eta) < 2.4) return 0.988984107971;
    if (pt < 1200.0 && fabs(eta) < 2.4) return 0.922808051109;
    return 1;
}
