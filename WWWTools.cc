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

        return fakerate_el_data(abseta, corrpt, syst);

        // fSumw[4][3]=0.170493 , x=85 , y=1.9895 , error=0.0370097
        // fSumw[4][2]=0.101855 , x=85 , y=1.1395 , error=0.0601733
        // fSumw[4][1]=0.199399 , x=85 , y=0.4    , error=0.0770095

        // fSumw[3][3]=0.1736   , x=40 , y=1.9895 , error=0.014733
        // fSumw[3][2]=0.164039 , x=40 , y=1.1395 , error=0.0213855
        // fSumw[3][1]=0.169051 , x=40 , y=0.4    , error=0.0207381

        // fSumw[2][3]=0.217486 , x=25 , y=1.9895 , error=0.0248149
        // fSumw[2][2]=0.191887 , x=25 , y=1.1395 , error=0.0232803
        // fSumw[2][1]=0.195776 , x=25 , y=0.4    , error=0.0282233

//        if ( corrpt >= 50. )
//        {
//            if      ( abseta >= 1.479 ) return 0.170493 + syst * 0.0370097;
//            else if ( abseta >= 0.8 )   return 0.101855 + syst * 0.0601733;
//            else                        return 0.199399 + syst * 0.0770095;
//        }
//        else if ( corrpt >= 30. )
//        {
//            if      ( abseta >= 1.479 ) return 0.1736   + syst * 0.014733 ;
//            else if ( abseta >= 0.8 )   return 0.164039 + syst * 0.0213855;
//            else                        return 0.169051 + syst * 0.0207381;
//        }
//        else if ( corrpt >= 20. )
//        {
//            if      ( abseta >= 1.479 ) return 0.217486 + syst * 0.0248149;
//            else if ( abseta >= 0.8 )   return 0.191887 + syst * 0.0232803;
//            else                        return 0.195776 + syst * 0.0282233;
//        }
//        else
//            return 0;
    }
    else if ( abs( wwwbaby.lep_pdgId()[idx] ) == 13 )
    {
        float relIso = wwwbaby.lep_relIso03EAv2()[idx];
        float coneptcorr = std::max( 0., relIso - 0.06 );
        float corrpt = wwwbaby.lep_p4()[idx].pt() * ( 1 + coneptcorr );
        float abseta = fabs( wwwbaby.lep_p4()[idx].eta() );

        return fakerate_mu_data(abseta, corrpt, syst);
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


//        if ( corrpt >= 50. )
//        {
//            if      ( abseta >= 2.1 ) return 0.0583592 + syst * 0.0227116 ;
//            else if ( abseta >= 1.2 ) return 0.0554221 + syst * 0.0171941 ;
//            else                      return 0.0317323 + syst * 0.00831419;
//        }
//        else if ( corrpt >= 30. )
//        {
//            if      ( abseta >= 2.1 ) return 0.0674571 + syst * 0.0151419 ;
//            else if ( abseta >= 1.2 ) return 0.0673591 + syst * 0.00792283;
//            else                      return 0.045541  + syst * 0.00546544;
//        }
//        else if ( corrpt >= 20. )
//        {
//            if      ( abseta >= 2.1 ) return 0.124618  + syst * 0.0126553 ;
//            else if ( abseta >= 1.2 ) return 0.109371  + syst * 0.00733365;
//            else                      return 0.0592033 + syst * 0.00397246;
//        }
//        else
//            return 0;

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

float lepsf_EGammaTightID(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 20.0 && fabs(eta) < -2.0) return 0.806557357311 + isyst * 0.0180196290542;
    if ((pt) < 35.0 && fabs(eta) < -2.0) return 0.88245934248 + isyst * 0.0101198259436;
    if ((pt) < 50.0 && fabs(eta) < -2.0) return 0.918840587139 + isyst * 0.00903908258991;
    if ((pt) < 90.0 && fabs(eta) < -2.0) return 0.940397322178 + isyst * 0.00700136278822;
    if ((pt) < 150.0 && fabs(eta) < -2.0) return 1.05108559132 + isyst * 0.0222742310261;
    if ((pt) < 500.0 && fabs(eta) < -2.0) return 1.05108559132 + isyst * 0.0224317920679;
    if ((pt) < 20.0 && fabs(eta) < -1.566) return 0.828571438789 + isyst * 0.0182790616195;
    if ((pt) < 35.0 && fabs(eta) < -1.566) return 0.927487373352 + isyst * 0.0176260918721;
    if ((pt) < 50.0 && fabs(eta) < -1.566) return 0.966987609863 + isyst * 0.00730306309988;
    if ((pt) < 90.0 && fabs(eta) < -1.566) return 0.980964481831 + isyst * 0.00561477853273;
    if ((pt) < 150.0 && fabs(eta) < -1.566) return 1.00592422485 + isyst * 0.0215595581354;
    if ((pt) < 500.0 && fabs(eta) < -1.566) return 1.00592422485 + isyst * 0.0217223036202;
    if ((pt) < 20.0 && fabs(eta) < -1.444) return 1.03286385536 + isyst * 0.106278969347;
    if ((pt) < 35.0 && fabs(eta) < -1.444) return 1.00759494305 + isyst * 0.109912699864;
    if ((pt) < 50.0 && fabs(eta) < -1.444) return 0.988195598125 + isyst * 0.0166870602707;
    if ((pt) < 90.0 && fabs(eta) < -1.444) return 0.995391726494 + isyst * 0.0240068358257;
    if ((pt) < 150.0 && fabs(eta) < -1.444) return 1.10410559177 + isyst * 0.0495688722273;
    if ((pt) < 500.0 && fabs(eta) < -1.444) return 1.10410559177 + isyst * 0.0496398733023;
    if ((pt) < 20.0 && fabs(eta) < -0.8) return 1.00775194168 + isyst * 0.0272663283234;
    if ((pt) < 35.0 && fabs(eta) < -0.8) return 0.972039461136 + isyst * 0.013293087796;
    if ((pt) < 50.0 && fabs(eta) < -0.8) return 0.974666655064 + isyst * 0.00686833231553;
    if ((pt) < 90.0 && fabs(eta) < -0.8) return 0.971674859524 + isyst * 0.0186735224927;
    if ((pt) < 150.0 && fabs(eta) < -0.8) return 0.98948597908 + isyst * 0.00867577650139;
    if ((pt) < 500.0 && fabs(eta) < -0.8) return 0.98948597908 + isyst * 0.00907265261538;
    if ((pt) < 20.0 && fabs(eta) < 0.0) return 0.94072163105 + isyst * 0.0257409592209;
    if ((pt) < 35.0 && fabs(eta) < 0.0) return 0.952998399734 + isyst * 0.0150816442813;
    if ((pt) < 50.0 && fabs(eta) < 0.0) return 0.953457474709 + isyst * 0.00472981124539;
    if ((pt) < 90.0 && fabs(eta) < 0.0) return 0.953086435795 + isyst * 0.0168857486142;
    if ((pt) < 150.0 && fabs(eta) < 0.0) return 0.975495934486 + isyst * 0.0126070521959;
    if ((pt) < 500.0 && fabs(eta) < 0.0) return 0.975495934486 + isyst * 0.0128833882441;
    if ((pt) < 20.0 && fabs(eta) < 0.8) return 0.945876300335 + isyst * 0.0257409592209;
    if ((pt) < 35.0 && fabs(eta) < 0.8) return 0.981996715069 + isyst * 0.0150816442813;
    if ((pt) < 50.0 && fabs(eta) < 0.8) return 0.979811549187 + isyst * 0.00472981124539;
    if ((pt) < 90.0 && fabs(eta) < 0.8) return 0.977528095245 + isyst * 0.0168857486142;
    if ((pt) < 150.0 && fabs(eta) < 0.8) return 1.0118906498 + isyst * 0.0127268449853;
    if ((pt) < 500.0 && fabs(eta) < 0.8) return 1.0118906498 + isyst * 0.0130006350174;
    if ((pt) < 20.0 && fabs(eta) < 1.444) return 0.989769816399 + isyst * 0.0272663283234;
    if ((pt) < 35.0 && fabs(eta) < 1.444) return 0.975206613541 + isyst * 0.013293087796;
    if ((pt) < 50.0 && fabs(eta) < 1.444) return 0.97463285923 + isyst * 0.00686833231553;
    if ((pt) < 90.0 && fabs(eta) < 1.444) return 0.97893434763 + isyst * 0.0186735224927;
    if ((pt) < 150.0 && fabs(eta) < 1.444) return 1.01072704792 + isyst * 0.00888880187216;
    if ((pt) < 500.0 && fabs(eta) < 1.444) return 1.01072704792 + isyst * 0.00927656867057;
    if ((pt) < 20.0 && fabs(eta) < 1.566) return 1.03398060799 + isyst * 0.106278969347;
    if ((pt) < 35.0 && fabs(eta) < 1.566) return 0.974874377251 + isyst * 0.109912699864;
    if ((pt) < 50.0 && fabs(eta) < 1.566) return 0.966386556625 + isyst * 0.0166870602707;
    if ((pt) < 90.0 && fabs(eta) < 1.566) return 0.979907274246 + isyst * 0.0240068358257;
    if ((pt) < 150.0 && fabs(eta) < 1.566) return 1.00716328621 + isyst * 0.0494349504371;
    if ((pt) < 500.0 && fabs(eta) < 1.566) return 1.00716328621 + isyst * 0.0495061435814;
    if ((pt) < 20.0 && fabs(eta) < 2.0) return 0.827485382557 + isyst * 0.0182790616195;
    if ((pt) < 35.0 && fabs(eta) < 2.0) return 0.908937633038 + isyst * 0.0176260918721;
    if ((pt) < 50.0 && fabs(eta) < 2.0) return 0.95782315731 + isyst * 0.00730306309988;
    if ((pt) < 90.0 && fabs(eta) < 2.0) return 0.968553483486 + isyst * 0.00561477853273;
    if ((pt) < 150.0 && fabs(eta) < 2.0) return 0.988262891769 + isyst * 0.0215595581354;
    if ((pt) < 500.0 && fabs(eta) < 2.0) return 0.988262891769 + isyst * 0.0217223036202;
    if ((pt) < 20.0 && fabs(eta) < 2.5) return 0.797014951706 + isyst * 0.0180196290542;
    if ((pt) < 35.0 && fabs(eta) < 2.5) return 0.863321781158 + isyst * 0.0101198259436;
    if ((pt) < 50.0 && fabs(eta) < 2.5) return 0.907801389694 + isyst * 0.00903908258991;
    if ((pt) < 90.0 && fabs(eta) < 2.5) return 0.937662363052 + isyst * 0.00700136278822;
    if ((pt) < 150.0 && fabs(eta) < 2.5) return 1.02135682106 + isyst * 0.0222742310261;
    if ((pt) < 500.0 && fabs(eta) < 2.5) return 1.02135682106 + isyst * 0.0224317920679;
    return 1;
}

float lepsf_EGammaReco(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 500.0 && fabs(eta) < -2.45) return 1.31760430336 + isyst * 0.0182387273858;
    if ((pt) < 500.0 && fabs(eta) < -2.4) return 1.11378002167 + isyst * 0.0110667218714;
    if ((pt) < 500.0 && fabs(eta) < -2.3) return 1.02462530136 + isyst * 0.00815834887034;
    if ((pt) < 500.0 && fabs(eta) < -2.2) return 1.013641119 + isyst * 0.00713343470588;
    if ((pt) < 500.0 && fabs(eta) < -2.0) return 1.00727653503 + isyst * 0.00420281345299;
    if ((pt) < 500.0 && fabs(eta) < -1.8) return 0.994818627834 + isyst * 0.00649266791838;
    if ((pt) < 500.0 && fabs(eta) < -1.63) return 0.994786262512 + isyst * 0.00516608222797;
    if ((pt) < 500.0 && fabs(eta) < -1.566) return 0.991631805897 + isyst * 0.00551198189805;
    if ((pt) < 500.0 && fabs(eta) < -1.444) return 0.963128507137 + isyst * 0.026030162586;
    if ((pt) < 500.0 && fabs(eta) < -1.2) return 0.989701330662 + isyst * 0.00359897172237;
    if ((pt) < 500.0 && fabs(eta) < -1.0) return 0.985655725002 + isyst * 0.00506368027448;
    if ((pt) < 500.0 && fabs(eta) < -0.6) return 0.981595098972 + isyst * 0.00331156908452;
    if ((pt) < 500.0 && fabs(eta) < -0.4) return 0.984678268433 + isyst * 0.00612870275792;
    if ((pt) < 500.0 && fabs(eta) < -0.2) return 0.981613874435 + isyst * 0.0063584778335;
    if ((pt) < 500.0 && fabs(eta) < 0.0) return 0.980432569981 + isyst * 0.00530156031977;
    if ((pt) < 500.0 && fabs(eta) < 0.2) return 0.984552025795 + isyst * 0.00530156031977;
    if ((pt) < 500.0 && fabs(eta) < 0.4) return 0.988764047623 + isyst * 0.0063584778335;
    if ((pt) < 500.0 && fabs(eta) < 0.6) return 0.987742602825 + isyst * 0.00612870275792;
    if ((pt) < 500.0 && fabs(eta) < 1.0) return 0.987742602825 + isyst * 0.00331156908452;
    if ((pt) < 500.0 && fabs(eta) < 1.2) return 0.987742602825 + isyst * 0.00506368027448;
    if ((pt) < 500.0 && fabs(eta) < 1.444) return 0.98767966032 + isyst * 0.00359897172237;
    if ((pt) < 500.0 && fabs(eta) < 1.566) return 0.967597782612 + isyst * 0.026030162586;
    if ((pt) < 500.0 && fabs(eta) < 1.63) return 0.989626526833 + isyst * 0.00551198189805;
    if ((pt) < 500.0 && fabs(eta) < 1.8) return 0.992761135101 + isyst * 0.00516608222797;
    if ((pt) < 500.0 && fabs(eta) < 2.0) return 0.991761088371 + isyst * 0.00649266791838;
    if ((pt) < 500.0 && fabs(eta) < 2.2) return 0.99794024229 + isyst * 0.00420281345299;
    if ((pt) < 500.0 && fabs(eta) < 2.3) return 1.00103735924 + isyst * 0.00713343470588;
    if ((pt) < 500.0 && fabs(eta) < 2.4) return 0.989506840706 + isyst * 0.00815834887034;
    if ((pt) < 500.0 && fabs(eta) < 2.45) return 0.970518887043 + isyst * 0.0110667218714;
    if ((pt) < 500.0 && fabs(eta) < 2.5) return 0.906666696072 + isyst * 0.0182387273858;
    return 1;
}

float lepsf_MuonMediumID_PeriodBCDEF(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 25.0 && fabs(eta) < 0.9) return 0.99326056242 + isyst * 0.00299304959529;
    if ((pt) < 30.0 && fabs(eta) < 0.9) return 0.993243694305 + isyst * 0.013774752384;
    if ((pt) < 40.0 && fabs(eta) < 0.9) return 0.999212145805 + isyst * 0.0122670982741;
    if ((pt) < 50.0 && fabs(eta) < 0.9) return 0.996031761169 + isyst * 0.000219301615914;
    if ((pt) < 60.0 && fabs(eta) < 0.9) return 0.992393016815 + isyst * 0.000816908857979;
    if ((pt) < 120.0 && fabs(eta) < 0.9) return 1.00521028042 + isyst * 0.00191670688622;
    if ((pt) < 25.0 && fabs(eta) < 1.2) return 0.987515449524 + isyst * 0.00328181231155;
    if ((pt) < 30.0 && fabs(eta) < 1.2) return 0.984555661678 + isyst * 0.00144660689391;
    if ((pt) < 40.0 && fabs(eta) < 1.2) return 0.989506840706 + isyst * 0.000443537527689;
    if ((pt) < 50.0 && fabs(eta) < 1.2) return 0.991288483143 + isyst * 0.000320477322702;
    if ((pt) < 60.0 && fabs(eta) < 1.2) return 0.990756630898 + isyst * 0.00107803144403;
    if ((pt) < 120.0 && fabs(eta) < 1.2) return 0.990235745907 + isyst * 0.00256092158742;
    if ((pt) < 25.0 && fabs(eta) < 2.1) return 0.992924451828 + isyst * 0.00168451232249;
    if ((pt) < 30.0 && fabs(eta) < 2.1) return 0.991163611412 + isyst * 0.000747527518378;
    if ((pt) < 40.0 && fabs(eta) < 2.1) return 0.993620336056 + isyst * 0.000263569134308;
    if ((pt) < 50.0 && fabs(eta) < 2.1) return 0.995641350746 + isyst * 0.00898828893907;
    if ((pt) < 60.0 && fabs(eta) < 2.1) return 0.99265640974 + isyst * 0.000791923382812;
    if ((pt) < 120.0 && fabs(eta) < 2.1) return 0.996506989002 + isyst * 0.00192651638706;
    if ((pt) < 25.0 && fabs(eta) < 2.4) return 0.972319483757 + isyst * 0.00281784438243;
    if ((pt) < 30.0 && fabs(eta) < 2.4) return 0.970559298992 + isyst * 0.00149029293773;
    if ((pt) < 40.0 && fabs(eta) < 2.4) return 0.969556510448 + isyst * 0.00063092914525;
    if ((pt) < 50.0 && fabs(eta) < 2.4) return 0.971171498299 + isyst * 0.0417734207699;
    if ((pt) < 60.0 && fabs(eta) < 2.4) return 0.970107495785 + isyst * 0.00226507269067;
    if ((pt) < 120.0 && fabs(eta) < 2.4) return 0.970472455025 + isyst * 0.00527212542888;
    return 1;
}

float lepsf_MuonMediumID_PeriodGH(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 25.0 && fabs(eta) < 0.9) return 1.00160014629 + isyst * 0.00248764356599;
    if ((pt) < 30.0 && fabs(eta) < 0.9) return 0.996136546135 + isyst * 0.0110811851285;
    if ((pt) < 40.0 && fabs(eta) < 0.9) return 1.00226068497 + isyst * 0.000461913558683;
    if ((pt) < 50.0 && fabs(eta) < 0.9) return 0.999405801296 + isyst * 0.000193286700103;
    if ((pt) < 60.0 && fabs(eta) < 0.9) return 0.994069099426 + isyst * 0.000606628841237;
    if ((pt) < 120.0 && fabs(eta) < 0.9) return 1.00124228001 + isyst * 0.00145716626641;
    if ((pt) < 25.0 && fabs(eta) < 1.2) return 1.00103020668 + isyst * 0.00673470115716;
    if ((pt) < 30.0 && fabs(eta) < 1.2) return 0.99916267395 + isyst * 0.0085269708973;
    if ((pt) < 40.0 && fabs(eta) < 1.2) return 0.99788492918 + isyst * 0.000437773546648;
    if ((pt) < 50.0 && fabs(eta) < 1.2) return 0.998400390148 + isyst * 0.000315730882773;
    if ((pt) < 60.0 && fabs(eta) < 1.2) return 0.994010567665 + isyst * 0.00112190091898;
    if ((pt) < 120.0 && fabs(eta) < 1.2) return 0.999541401863 + isyst * 0.0027155142224;
    if ((pt) < 25.0 && fabs(eta) < 2.1) return 0.997374892235 + isyst * 0.0017970830944;
    if ((pt) < 30.0 && fabs(eta) < 2.1) return 0.996894299984 + isyst * 0.000833129881249;
    if ((pt) < 40.0 && fabs(eta) < 2.1) return 0.99773812294 + isyst * 0.000286574846273;
    if ((pt) < 50.0 && fabs(eta) < 2.1) return 0.996386408806 + isyst * 0.000210286445003;
    if ((pt) < 60.0 && fabs(eta) < 2.1) return 0.991427123547 + isyst * 0.000915349227664;
    if ((pt) < 120.0 && fabs(eta) < 2.1) return 0.995427966118 + isyst * 0.00230696408947;
    if ((pt) < 25.0 && fabs(eta) < 2.4) return 0.983996212482 + isyst * 0.00273030602568;
    if ((pt) < 30.0 && fabs(eta) < 2.4) return 0.97872120142 + isyst * 0.0014648021749;
    if ((pt) < 40.0 && fabs(eta) < 2.4) return 0.973877429962 + isyst * 0.000728211424762;
    if ((pt) < 50.0 && fabs(eta) < 2.4) return 0.974616527557 + isyst * 0.000733926213164;
    if ((pt) < 60.0 && fabs(eta) < 2.4) return 0.970205843449 + isyst * 0.00236731018144;
    if ((pt) < 120.0 && fabs(eta) < 2.4) return 0.971422970295 + isyst * 0.0105361766494;
    return 1;
}

float lepsf_MuonReco(float eta, float nothing, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if (fabs(eta) < 0.2 && fabs(nothing) < 1000000) return 0.996996496582 + isyst * 7.27173381037e-05;
    if (fabs(eta) < 0.4 && fabs(nothing) < 1000000) return 0.99771175854 + isyst * 8.20493871853e-05;
    if (fabs(eta) < 0.6 && fabs(nothing) < 1000000) return 0.99807762769 + isyst * 7.15881829855e-05;
    if (fabs(eta) < 0.8 && fabs(nothing) < 1000000) return 0.997803871488 + isyst * 7.4350908861e-05;
    if (fabs(eta) < 1.0 && fabs(nothing) < 1000000) return 0.99797080488 + isyst * 0.000106158649569;
    if (fabs(eta) < 1.2 && fabs(nothing) < 1000000) return 0.997147740951 + isyst * 0.000186466741555;
    if (fabs(eta) < 1.4 && fabs(nothing) < 1000000) return 0.996227437404 + isyst * 0.000181630468798;
    if (fabs(eta) < 1.6 && fabs(nothing) < 1000000) return 0.995478609224 + isyst * 0.00017062499716;
    if (fabs(eta) < 1.8 && fabs(nothing) < 1000000) return 0.995780818567 + isyst * 0.000176732421664;
    if (fabs(eta) < 2.0 && fabs(nothing) < 1000000) return 0.99389185437 + isyst * 0.000234784307893;
    if (fabs(eta) < 2.2 && fabs(nothing) < 1000000) return 0.992942696013 + isyst * 0.000330644022311;
    if (fabs(eta) < 2.4 && fabs(nothing) < 1000000) return 0.987313312455 + isyst * 0.000859346400986;
    return 1;
}

float lepsf_EGammaTightPOG_EGammaVVV(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 30.0 && fabs(eta) < 1.0) return 0.980376899242 + isyst * 0.00382222977557;
    if ((pt) < 40.0 && fabs(eta) < 1.0) return 0.981454670429 + isyst * 0.00193882357048;
    if ((pt) < 50.0 && fabs(eta) < 1.0) return 0.980952382088 + isyst * 0.00163532483243;
    if ((pt) < 2000.0 && fabs(eta) < 1.0) return 0.976341068745 + isyst * 0.00282505650216;
    if ((pt) < 30.0 && fabs(eta) < 1.4442) return 0.945449829102 + isyst * 0.00616297079806;
    if ((pt) < 40.0 && fabs(eta) < 1.4442) return 0.94911891222 + isyst * 0.00332032789849;
    if ((pt) < 50.0 && fabs(eta) < 1.4442) return 0.95356208086 + isyst * 0.00292671435923;
    if ((pt) < 2000.0 && fabs(eta) < 1.4442) return 0.949157774448 + isyst * 0.00457991622718;
    if ((pt) < 30.0 && fabs(eta) < 1.56) return 0.938438951969 + isyst * 0.0179865730823;
    if ((pt) < 40.0 && fabs(eta) < 1.56) return 0.941834270954 + isyst * 0.00994399032801;
    if ((pt) < 50.0 && fabs(eta) < 1.56) return 0.947787821293 + isyst * 0.00809667245743;
    if ((pt) < 2000.0 && fabs(eta) < 1.56) return 0.952665507793 + isyst * 0.0144871613152;
    if ((pt) < 30.0 && fabs(eta) < 2.0) return 0.963943779469 + isyst * 0.00702015595817;
    if ((pt) < 40.0 && fabs(eta) < 2.0) return 0.971829116344 + isyst * 0.00401011091053;
    if ((pt) < 50.0 && fabs(eta) < 2.0) return 0.975824117661 + isyst * 0.00357086745336;
    if ((pt) < 2000.0 && fabs(eta) < 2.0) return 0.980078279972 + isyst * 0.0062819436647;
    if ((pt) < 30.0 && fabs(eta) < 2.5) return 0.984203636646 + isyst * 0.0080419549299;
    if ((pt) < 40.0 && fabs(eta) < 2.5) return 0.998542606831 + isyst * 0.00478520093589;
    if ((pt) < 50.0 && fabs(eta) < 2.5) return 1.00013875961 + isyst * 0.0048532423643;
    if ((pt) < 2000.0 && fabs(eta) < 2.5) return 1.00487089157 + isyst * 0.00900235705955;
    return 1;
}

float lepsf_EGammaVVV_EGammaVVVEle12(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 30.0 && fabs(eta) < 1.0) return 0.938784837723 + isyst * 0.00219590521361;
    if ((pt) < 40.0 && fabs(eta) < 1.0) return 0.956399559975 + isyst * 0.00127861103515;
    if ((pt) < 50.0 && fabs(eta) < 1.0) return 0.962966561317 + isyst * 0.000316125299188;
    if ((pt) < 2000.0 && fabs(eta) < 1.0) return 0.967466413975 + isyst * 0.00155382971572;
    if ((pt) < 30.0 && fabs(eta) < 1.4442) return 0.94966250658 + isyst * 0.00368355725231;
    if ((pt) < 40.0 && fabs(eta) < 1.4442) return 0.96033769846 + isyst * 0.00197332068038;
    if ((pt) < 50.0 && fabs(eta) < 1.4442) return 0.965217530727 + isyst * 0.000220183866521;
    if ((pt) < 2000.0 && fabs(eta) < 1.4442) return 0.969266951084 + isyst * 0.00265256114334;
    if ((pt) < 30.0 && fabs(eta) < 1.56) return 0.957857549191 + isyst * 0.0115744391104;
    if ((pt) < 40.0 && fabs(eta) < 1.56) return 0.962614357471 + isyst * 0.0059518655469;
    if ((pt) < 50.0 && fabs(eta) < 1.56) return 0.964548945427 + isyst * 0.0055577892749;
    if ((pt) < 2000.0 && fabs(eta) < 1.56) return 0.962951600552 + isyst * 0.00883849432663;
    if ((pt) < 30.0 && fabs(eta) < 2.0) return 0.958219587803 + isyst * 0.00428889366748;
    if ((pt) < 40.0 && fabs(eta) < 2.0) return 0.96611481905 + isyst * 0.00169359181774;
    if ((pt) < 50.0 && fabs(eta) < 2.0) return 0.971387922764 + isyst * 0.00208029227053;
    if ((pt) < 2000.0 && fabs(eta) < 2.0) return 0.974559724331 + isyst * 1.84578668625e-05;
    if ((pt) < 30.0 && fabs(eta) < 2.5) return 0.910284698009 + isyst * 0.00477846357209;
    if ((pt) < 40.0 && fabs(eta) < 2.5) return 0.947666943073 + isyst * 0.00277335535449;
    if ((pt) < 50.0 && fabs(eta) < 2.5) return 0.955771327019 + isyst * 0.00289434846743;
    if ((pt) < 2000.0 && fabs(eta) < 2.5) return 0.962957322598 + isyst * 0.00532928187002;
    return 1;
}

float lepsf_EGammaVVV_EGammaVVVEleLead(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 30.0 && fabs(eta) < 1.0) return 0.902300655842 + isyst * 0.00213228896962;
    if ((pt) < 40.0 && fabs(eta) < 1.0) return 0.956393837929 + isyst * 0.00130303273575;
    if ((pt) < 50.0 && fabs(eta) < 1.0) return 0.962963223457 + isyst * 0.000868854723185;
    if ((pt) < 2000.0 && fabs(eta) < 1.0) return 0.967466413975 + isyst * 0.00155382971572;
    if ((pt) < 30.0 && fabs(eta) < 1.4442) return 0.903473079205 + isyst * 0.0035448160487;
    if ((pt) < 40.0 && fabs(eta) < 1.4442) return 0.960295259953 + isyst * 0.00197019492648;
    if ((pt) < 50.0 && fabs(eta) < 1.4442) return 0.965217530727 + isyst * 0.000220183866521;
    if ((pt) < 2000.0 && fabs(eta) < 1.4442) return 0.969266951084 + isyst * 0.00265256114334;
    if ((pt) < 30.0 && fabs(eta) < 1.56) return 0.914313673973 + isyst * 0.0113559311569;
    if ((pt) < 40.0 && fabs(eta) < 1.56) return 0.962467074394 + isyst * 0.00602194977727;
    if ((pt) < 50.0 && fabs(eta) < 1.56) return 0.964522898197 + isyst * 0.00558842549423;
    if ((pt) < 2000.0 && fabs(eta) < 1.56) return 0.962951600552 + isyst * 0.00883849432663;
    if ((pt) < 30.0 && fabs(eta) < 2.0) return 0.920680880547 + isyst * 0.00415683064104;
    if ((pt) < 40.0 && fabs(eta) < 2.0) return 0.966320037842 + isyst * 0.00231408855426;
    if ((pt) < 50.0 && fabs(eta) < 2.0) return 0.971381902695 + isyst * 0.00210118701395;
    if ((pt) < 2000.0 && fabs(eta) < 2.0) return 0.970944464207 + isyst * 1.80491030395e-05;
    if ((pt) < 30.0 && fabs(eta) < 2.5) return 0.867845416069 + isyst * 0.00461339431639;
    if ((pt) < 40.0 && fabs(eta) < 2.5) return 0.947629570961 + isyst * 0.00278310457802;
    if ((pt) < 50.0 && fabs(eta) < 2.5) return 0.955740690231 + isyst * 0.00208964389036;
    if ((pt) < 2000.0 && fabs(eta) < 2.5) return 0.962957322598 + isyst * 0.00532928187002;
    return 1;
}

float lepsf_MuMediumPOG_MuTightVVV(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 25.0 && fabs(eta) < 0.9) return 0.972987830639 + isyst * 0.0038003901295;
    if ((pt) < 30.0 && fabs(eta) < 0.9) return 0.980633199215 + isyst * 0.00260697870217;
    if ((pt) < 40.0 && fabs(eta) < 0.9) return 0.986259222031 + isyst * 0.00133423658569;
    if ((pt) < 50.0 && fabs(eta) < 0.9) return 0.989767193794 + isyst * 0.00124355006164;
    if ((pt) < 60.0 && fabs(eta) < 0.9) return 0.990144252777 + isyst * 0.00253156629733;
    if ((pt) < 120.0 && fabs(eta) < 0.9) return 0.992544472218 + isyst * 0.00331781595267;
    if ((pt) < 200.0 && fabs(eta) < 0.9) return 1.00607073307 + isyst * 0.0149094139022;
    if ((pt) < 1200.0 && fabs(eta) < 0.9) return 1.00708353519 + isyst * 0.0335569509038;
    if ((pt) < 25.0 && fabs(eta) < 1.2) return 0.982158899307 + isyst * 0.00740547808971;
    if ((pt) < 30.0 && fabs(eta) < 1.2) return 0.99262791872 + isyst * 0.00511323464457;
    if ((pt) < 40.0 && fabs(eta) < 1.2) return 0.994507312775 + isyst * 0.00258035645753;
    if ((pt) < 50.0 && fabs(eta) < 1.2) return 0.994486868382 + isyst * 0.00198296279702;
    if ((pt) < 60.0 && fabs(eta) < 1.2) return 0.995231330395 + isyst * 0.00452383601426;
    if ((pt) < 120.0 && fabs(eta) < 1.2) return 0.996162414551 + isyst * 0.00598093468074;
    if ((pt) < 200.0 && fabs(eta) < 1.2) return 0.998247504234 + isyst * 0.0277592866625;
    if ((pt) < 1200.0 && fabs(eta) < 1.2) return 1.00434410572 + isyst * 0.0576382502384;
    if ((pt) < 25.0 && fabs(eta) < 2.1) return 0.992021918297 + isyst * 0.00412251574796;
    if ((pt) < 30.0 && fabs(eta) < 2.1) return 0.999400854111 + isyst * 0.00313541238303;
    if ((pt) < 40.0 && fabs(eta) < 2.1) return 0.999857544899 + isyst * 0.00174861771899;
    if ((pt) < 50.0 && fabs(eta) < 2.1) return 0.999814629555 + isyst * 0.00147727821042;
    if ((pt) < 60.0 && fabs(eta) < 2.1) return 1.00330460072 + isyst * 0.00294343817454;
    if ((pt) < 120.0 && fabs(eta) < 2.1) return 1.00755393505 + isyst * 0.00406308933367;
    if ((pt) < 200.0 && fabs(eta) < 2.1) return 1.01337754726 + isyst * 0.0191509926106;
    if ((pt) < 1200.0 && fabs(eta) < 2.1) return 1.02070629597 + isyst * 0.0465072322835;
    if ((pt) < 25.0 && fabs(eta) < 2.4) return 0.969161331654 + isyst * 0.00705886719455;
    if ((pt) < 30.0 && fabs(eta) < 2.4) return 0.975010037422 + isyst * 0.00640140774693;
    if ((pt) < 40.0 && fabs(eta) < 2.4) return 0.990625023842 + isyst * 0.00328420129413;
    if ((pt) < 50.0 && fabs(eta) < 2.4) return 0.994753956795 + isyst * 0.00324093861002;
    if ((pt) < 60.0 && fabs(eta) < 2.4) return 0.992071866989 + isyst * 0.00636933685301;
    if ((pt) < 120.0 && fabs(eta) < 2.4) return 0.994121670723 + isyst * 0.0088715023098;
    if ((pt) < 200.0 && fabs(eta) < 2.4) return 0.994795560837 + isyst * 0.0453345991072;
    if ((pt) < 1200.0 && fabs(eta) < 2.4) return 1.01750481129 + isyst * 0.499230818045;
    return 1;
}

float lepsf_MuTightVVV_MuTightVVVMu8(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 25.0 && fabs(eta) < 0.9) return 0.933527767658 + isyst * 0.00252632647732;
    if ((pt) < 30.0 && fabs(eta) < 0.9) return 0.933456897736 + isyst * 0.00152941158904;
    if ((pt) < 40.0 && fabs(eta) < 0.9) return 0.934123814106 + isyst * 5.7151570156e-05;
    if ((pt) < 50.0 && fabs(eta) < 0.9) return 0.933176159859 + isyst * 0.000533302705251;
    if ((pt) < 60.0 && fabs(eta) < 0.9) return 0.931691288948 + isyst * 0.00134394291184;
    if ((pt) < 120.0 && fabs(eta) < 0.9) return 0.929507613182 + isyst * 0.00175162047452;
    if ((pt) < 200.0 && fabs(eta) < 0.9) return 0.924368977547 + isyst * 0.00746803644143;
    if ((pt) < 1200.0 && fabs(eta) < 0.9) return 0.910157442093 + isyst * 0.0167671085273;
    if ((pt) < 25.0 && fabs(eta) < 1.2) return 0.96690428257 + isyst * 0.00463995182995;
    if ((pt) < 30.0 && fabs(eta) < 1.2) return 0.968839883804 + isyst * 0.000296189435009;
    if ((pt) < 40.0 && fabs(eta) < 1.2) return 0.968176484108 + isyst * 0.00147301227187;
    if ((pt) < 50.0 && fabs(eta) < 1.2) return 0.967356085777 + isyst * 0.0012252910266;
    if ((pt) < 60.0 && fabs(eta) < 1.2) return 0.966672122478 + isyst * 0.00253712385559;
    if ((pt) < 120.0 && fabs(eta) < 1.2) return 0.965913832188 + isyst * 0.00324173670128;
    if ((pt) < 200.0 && fabs(eta) < 1.2) return 0.962421000004 + isyst * 0.0145629366561;
    if ((pt) < 1200.0 && fabs(eta) < 1.2) return 0.961861371994 + isyst * 0.0236944487724;
    if ((pt) < 25.0 && fabs(eta) < 2.1) return 0.961845993996 + isyst * 0.00034315722064;
    if ((pt) < 30.0 && fabs(eta) < 2.1) return 0.962982833385 + isyst * 0.00180503994187;
    if ((pt) < 40.0 && fabs(eta) < 2.1) return 0.962533831596 + isyst * 0.000139178337978;
    if ((pt) < 50.0 && fabs(eta) < 2.1) return 0.963077664375 + isyst * 0.0008053328847;
    if ((pt) < 60.0 && fabs(eta) < 2.1) return 0.963417887688 + isyst * 0.00160263939342;
    if ((pt) < 120.0 && fabs(eta) < 2.1) return 0.962129116058 + isyst * 0.000363735247895;
    if ((pt) < 200.0 && fabs(eta) < 2.1) return 0.935746133327 + isyst * 0.0082700842047;
    if ((pt) < 1200.0 && fabs(eta) < 2.1) return 0.89957177639 + isyst * 0.0245530433342;
    if ((pt) < 25.0 && fabs(eta) < 2.4) return 0.944805383682 + isyst * 0.00573946175886;
    if ((pt) < 30.0 && fabs(eta) < 2.4) return 0.948861658573 + isyst * 0.0148002473209;
    if ((pt) < 40.0 && fabs(eta) < 2.4) return 0.949180960655 + isyst * 0.00205940031199;
    if ((pt) < 50.0 && fabs(eta) < 2.4) return 0.949977576733 + isyst * 0.00186501847954;
    if ((pt) < 60.0 && fabs(eta) < 2.4) return 0.947516798973 + isyst * 0.00371860765466;
    if ((pt) < 120.0 && fabs(eta) < 2.4) return 0.948109149933 + isyst * 0.00521480251546;
    if ((pt) < 200.0 && fabs(eta) < 2.4) return 0.926202654839 + isyst * 0.0195444637705;
    if ((pt) < 1200.0 && fabs(eta) < 2.4) return 0.865132570267 + isyst * 0.532614728793;
    return 1;
}

float lepsf_MuTightVVV_MuTightVVVMu17(float pt, float eta, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if ((pt) < 25.0 && fabs(eta) < 0.9) return 0.922722399235 + isyst * 0.00250782468408;
    if ((pt) < 30.0 && fabs(eta) < 0.9) return 0.923307836056 + isyst * 0.000366351007182;
    if ((pt) < 40.0 && fabs(eta) < 0.9) return 0.924030900002 + isyst * 0.000730909831923;
    if ((pt) < 50.0 && fabs(eta) < 0.9) return 0.923241138458 + isyst * 0.000768592171557;
    if ((pt) < 60.0 && fabs(eta) < 0.9) return 0.921180009842 + isyst * 0.00080159636303;
    if ((pt) < 120.0 && fabs(eta) < 0.9) return 0.918352901936 + isyst * 0.00185961188485;
    if ((pt) < 200.0 && fabs(eta) < 0.9) return 0.910844266415 + isyst * 0.00739344058835;
    if ((pt) < 1200.0 && fabs(eta) < 0.9) return 0.894463002682 + isyst * 0.016525458852;
    if ((pt) < 25.0 && fabs(eta) < 1.2) return 0.942260980606 + isyst * 0.00458041035964;
    if ((pt) < 30.0 && fabs(eta) < 1.2) return 0.945485949516 + isyst * 0.00346783935787;
    if ((pt) < 40.0 && fabs(eta) < 1.2) return 0.946832001209 + isyst * 0.00141910674422;
    if ((pt) < 50.0 && fabs(eta) < 1.2) return 0.946245908737 + isyst * 0.000673513472744;
    if ((pt) < 60.0 && fabs(eta) < 1.2) return 0.945586264133 + isyst * 0.00240884834276;
    if ((pt) < 120.0 && fabs(eta) < 1.2) return 0.945138394833 + isyst * 0.00318644606928;
    if ((pt) < 200.0 && fabs(eta) < 1.2) return 0.933813333511 + isyst * 0.014169270918;
    if ((pt) < 1200.0 && fabs(eta) < 1.2) return 0.931742787361 + isyst * 0.0243677246898;
    if ((pt) < 25.0 && fabs(eta) < 2.1) return 0.934685409069 + isyst * 0.00235550744943;
    if ((pt) < 30.0 && fabs(eta) < 2.1) return 0.940949559212 + isyst * 0.0017736671693;
    if ((pt) < 40.0 && fabs(eta) < 2.1) return 0.941943407059 + isyst * 0.000142803887399;
    if ((pt) < 50.0 && fabs(eta) < 2.1) return 0.943855106831 + isyst * 0.000241549432239;
    if ((pt) < 60.0 && fabs(eta) < 2.1) return 0.94069480896 + isyst * 0.00167134522243;
    if ((pt) < 120.0 && fabs(eta) < 2.1) return 0.926266431808 + isyst * 0.00212335629489;
    if ((pt) < 200.0 && fabs(eta) < 2.1) return 0.873467862606 + isyst * 0.00932832611851;
    if ((pt) < 1200.0 && fabs(eta) < 2.1) return 0.823928773403 + isyst * 0.0230251207295;
    if ((pt) < 25.0 && fabs(eta) < 2.4) return 0.837569534779 + isyst * 0.00488000441924;
    if ((pt) < 30.0 && fabs(eta) < 2.4) return 0.866438448429 + isyst * 0.00319838676876;
    if ((pt) < 40.0 && fabs(eta) < 2.4) return 0.885224223137 + isyst * 0.00194770815806;
    if ((pt) < 50.0 && fabs(eta) < 2.4) return 0.898686885834 + isyst * 0.00179089574024;
    if ((pt) < 60.0 && fabs(eta) < 2.4) return 0.901396632195 + isyst * 0.00355735499739;
    if ((pt) < 120.0 && fabs(eta) < 2.4) return 0.906494557858 + isyst * 0.00504503993108;
    if ((pt) < 200.0 && fabs(eta) < 2.4) return 0.884702563286 + isyst * 0.0182467487213;
    if ((pt) < 1200.0 && fabs(eta) < 2.4) return 0.80960047245 + isyst * 0.0950186278378;
    return 1;
}

float fakerate_mu_data(float eta, float conecorrpt, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if (fabs(eta) < 1.2 && (conecorrpt) < 20.0) return 0.18162853637 + isyst * 0.00419139178595;
    if (fabs(eta) < 2.1 && (conecorrpt) < 20.0) return 0.228556342818 + isyst * 0.00643408235926;
    if (fabs(eta) < 2.4 && (conecorrpt) < 20.0) return 0.284501936043 + isyst * 0.0138531422818;
    if (fabs(eta) < 1.2 && (conecorrpt) < 30.0) return 0.0593080016201 + isyst * 0.00429931883058;
    if (fabs(eta) < 2.1 && (conecorrpt) < 30.0) return 0.0812573899031 + isyst * 0.00690133585711;
    if (fabs(eta) < 2.4 && (conecorrpt) < 30.0) return 0.120471940997 + isyst * 0.0154600921933;
    if (fabs(eta) < 1.2 && (conecorrpt) < 50.0) return 0.0502796544514 + isyst * 0.00337692471616;
    if (fabs(eta) < 2.1 && (conecorrpt) < 50.0) return 0.0755524005305 + isyst * 0.00454099610916;
    if (fabs(eta) < 2.4 && (conecorrpt) < 50.0) return 0.102110577821 + isyst * 0.0116371933399;
    if (fabs(eta) < 1.2 && (conecorrpt) < 10000) return 0.0504603068665 + isyst * 0.0223031276198;
    if (fabs(eta) < 2.1 && (conecorrpt) < 10000) return 0.0699800116561 + isyst * 0.0249675245212;
    if (fabs(eta) < 2.4 && (conecorrpt) < 10000) return 0.0699800116561 + isyst * 0.0249675245212;
    return 1;
}

float fakerate_el_data(float eta, float conecorrpt, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if (fabs(eta) < 0.8 && (conecorrpt) < 20.0) return 0.180862610834 + isyst * 0.011780565611;
    if (fabs(eta) < 1.479 && (conecorrpt) < 20.0) return 0.151654211681 + isyst * 0.00866843935974;
    if (fabs(eta) < 2.5 && (conecorrpt) < 20.0) return 0.166856924948 + isyst * 0.0086750832925;
    if (fabs(eta) < 0.8 && (conecorrpt) < 30.0) return 0.133459584075 + isyst * 0.0146960861119;
    if (fabs(eta) < 1.479 && (conecorrpt) < 30.0) return 0.108987284034 + isyst * 0.0112594175795;
    if (fabs(eta) < 2.5 && (conecorrpt) < 30.0) return 0.127048887898 + isyst * 0.0116391553635;
    if (fabs(eta) < 0.8 && (conecorrpt) < 50.0) return 0.115077619378 + isyst * 0.0116121608794;
    if (fabs(eta) < 1.479 && (conecorrpt) < 50.0) return 0.0997162232625 + isyst * 0.0099509146269;
    if (fabs(eta) < 2.5 && (conecorrpt) < 50.0) return 0.106395503788 + isyst * 0.00678263788235;
    if (fabs(eta) < 0.8 && (conecorrpt) < 10000) return 0.0793935587448 + isyst * 0.0340196244843;
    if (fabs(eta) < 1.479 && (conecorrpt) < 10000) return 0.0794240499611 + isyst * 0.0289512048486;
    if (fabs(eta) < 2.5 && (conecorrpt) < 10000) return 0.0995944944707 + isyst * 0.0122391050471;
    return 1;
}

float fakerate_mu_qcd(float eta, float conecorrpt, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if (fabs(eta) < 1.2 && (conecorrpt) < 20.0) return 0.192454320571 + isyst * 0.00435651562103;
    if (fabs(eta) < 2.1 && (conecorrpt) < 20.0) return 0.25212783352 + isyst * 0.00673523309854;
    if (fabs(eta) < 2.4 && (conecorrpt) < 20.0) return 0.299959325364 + isyst * 0.0139496576752;
    if (fabs(eta) < 1.2 && (conecorrpt) < 30.0) return 0.0592032661673 + isyst * 0.00397245526891;
    if (fabs(eta) < 2.1 && (conecorrpt) < 30.0) return 0.109370893274 + isyst * 0.00733365179432;
    if (fabs(eta) < 2.4 && (conecorrpt) < 30.0) return 0.124618041094 + isyst * 0.0126552722177;
    if (fabs(eta) < 1.2 && (conecorrpt) < 50.0) return 0.0455410037998 + isyst * 0.00546544224618;
    if (fabs(eta) < 2.1 && (conecorrpt) < 50.0) return 0.0673590936047 + isyst * 0.00792283130924;
    if (fabs(eta) < 2.4 && (conecorrpt) < 50.0) return 0.067457120855 + isyst * 0.0151418937682;
    if (fabs(eta) < 1.2 && (conecorrpt) < 10000) return 0.0317322626443 + isyst * 0.00831418957152;
    if (fabs(eta) < 2.1 && (conecorrpt) < 10000) return 0.0554220997027 + isyst * 0.0171941451909;
    if (fabs(eta) < 2.4 && (conecorrpt) < 10000) return 0.0583591616584 + isyst * 0.0227115872818;
    return 1;
}

float fakerate_el_qcd(float eta, float conecorrpt, int isyst)
{
    if (isyst != 1 && isyst != -1 && isyst != 0)
        printf(Form("WARNING - in function=%s, isyst=%d is not recommended!\n", __FUNCTION__, isyst));
    if (fabs(eta) < 0.8 && (conecorrpt) < 20.0) return 0.230089425906 + isyst * 0.013635548062;
    if (fabs(eta) < 1.479 && (conecorrpt) < 20.0) return 0.187918774855 + isyst * 0.0301966750315;
    if (fabs(eta) < 2.5 && (conecorrpt) < 20.0) return 0.181326427551 + isyst * 0.0229038017929;
    if (fabs(eta) < 0.8 && (conecorrpt) < 30.0) return 0.173060172437 + isyst * 0.0540543730729;
    if (fabs(eta) < 1.479 && (conecorrpt) < 30.0) return 0.127839209154 + isyst * 0.0578526851684;
    if (fabs(eta) < 2.5 && (conecorrpt) < 30.0) return 0.106022264846 + isyst * 0.0311610831766;
    if (fabs(eta) < 0.8 && (conecorrpt) < 50.0) return 0.0515026065715 + isyst * 0.0163129257811;
    if (fabs(eta) < 1.479 && (conecorrpt) < 50.0) return 0.141963824743 + isyst * 0.0705629823265;
    if (fabs(eta) < 2.5 && (conecorrpt) < 50.0) return 0.0987933061792 + isyst * 0.022152759532;
    if (fabs(eta) < 0.8 && (conecorrpt) < 10000) return 0.0854687145457 + isyst * 0.0367682367283;
    if (fabs(eta) < 1.479 && (conecorrpt) < 10000) return 0.101560108062 + isyst * 0.0386388851491;
    if (fabs(eta) < 2.5 && (conecorrpt) < 10000) return 0.0892398734513 + isyst * 0.0199718234062;
    return 1;
}


