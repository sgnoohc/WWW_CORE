#include "WWWTools.h"

unsigned int objidx_set_to_eventid;
ObjIdx lepidx;
ObjIdx jetidx;

//______________________________________________________________________________________
bool passSSpresel( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1 ) || ( lepidx["TightLepton"].size()                      ==  2  ) )) return false;
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==  2  ) )) return false;
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()   >  30. ) )) return false;
    if (!( ( Nm1idx == 4 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()   >  30. ) )) return false;
    if (!( ( Nm1idx == 5 ) || ( jetidx["GoodSSJet"]    .size()                    >=  2  ) )) return false;
    if (!( ( Nm1idx == 6 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==  0  ) )) return false;
    if (!( ( Nm1idx == 7 ) || ( isSS()                                                   ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSScommon( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSSpresel()                       ) )) return false;
    if (!( ( Nm1idx == 1 ) || ( MjjW()                     > 60.  &&
                                MjjW()                     < 100.    ) )) return false;
    if (!( ( Nm1idx == 2 ) || ( MjjLead()                  < 400.    ) )) return false;
    if (!( ( Nm1idx == 3 ) || ( DEtajjLead()               < 1.5     ) )) return false;
    if (!( ( Nm1idx == 4 ) || ( jetidx["LooseBJet"].size() == 0      ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSMM( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSScommon() ) )) return false;
    if (!( ( Nm1idx == 1 ) || ( isSSMM()       ) )) return false;
    if (!( ( Nm1idx == 2 ) || ( MllSS() > 40.  ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSEM( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSScommon()         ) )) return false;
    if (!( ( Nm1idx == 1 ) || ( isSSEM()               ) )) return false;
    if (!( ( Nm1idx == 2 ) || ( MllSS()          > 30. ) )) return false;
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.met_pt() > 40. ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool passSSEE( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSScommon()             ) )) return false;
    if (!( ( Nm1idx == 1 ) || ( isSSEE()                   ) )) return false;
    if (!( ( Nm1idx == 2 ) || ( MllSS()          >  40.    ) )) return false;
    if (!( ( Nm1idx == 3 ) || ( MllSS()          <  80. ||
                                MllSS()          > 100.    ) )) return false;
    if (!( ( Nm1idx == 4 ) || ( wwwbaby.met_pt() >  40.    ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3Lpresel( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1 ) || ( lepidx["TightLepton"].size()                      ==  3  ) )) return false;
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==  3  ) )) return false;
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()   >  20. ) )) return false;
    if (!( ( Nm1idx == 4 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()   >  20. ) )) return false;
    if (!( ( Nm1idx == 5 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][2]].pt()   >  20. ) )) return false;
    if (!( ( Nm1idx == 6 ) || ( abs( totalCharge() )                              ==  1  ) )) return false;
    if (!( ( Nm1idx == 7 ) || ( jetidx["Good3LJet"].size()                        <=  1  ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3Lcommon( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( pass3Lpresel()                   ) )) return false;
    if (!( ( Nm1idx == 1 ) || ( Pt3l() > 60.                     ) )) return false;
    if (!( ( Nm1idx == 2 ) || ( jetidx["LooseBJet"].size() ==  0 ) )) return false;
    return true;
}

//______________________________________________________________________________________
bool pass3L0SFOS( int Nm1idx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( pass3Lcommon()            ) )) return false;
    if (!( ( Nm1idx == 1 ) || ( getNumSFOS()      ==   0  ) )) return false;
    return true;
}

//    if (!( ( Nm1idx == 3 ) || ( wwwbaby.met_pt()  >   50. ) )) return false;
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
}

//______________________________________________________________________________________
bool isObjectSelected()
{
    if ( objidx_set_to_eventid == wwwbaby.evt() )
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
        if ( isTightLepton( ilep ) ) idx["TightLepton"].push_back( ilep );
        if ( isTightMuon  ( ilep ) ) idx["TightMuon"]  .push_back( ilep );
        if ( isTightElec  ( ilep ) ) idx["TightElec"]  .push_back( ilep );
        if ( isLooseLepton( ilep ) ) idx["LooseLepton"].push_back( ilep );
        if ( isLooseMuon  ( ilep ) ) idx["LooseMuon"]  .push_back( ilep );
        if ( isLooseElec  ( ilep ) ) idx["LooseElec"]  .push_back( ilep );
        if ( isLbntLepton ( ilep ) ) idx["LbntLepton"] .push_back( ilep );
        if ( isLbntMuon   ( ilep ) ) idx["LbntMuon"]   .push_back( ilep );
        if ( isLbntElec   ( ilep ) ) idx["LbntElec"]   .push_back( ilep );
    }
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
    for ( unsigned int ijet = 0; ijet < wwwbaby.jets_p4().size(); ++ijet )
    {
        for ( unsigned int jjet = ijet; jjet < wwwbaby.jets_p4().size(); ++jjet )
        {
            if ( ijet == jjet )
                continue;
            LorentzVector ijet_p4 = wwwbaby.jets_p4()[ijet];
            LorentzVector jjet_p4 = wwwbaby.jets_p4()[jjet];
            float dr = ROOT::Math::VectorUtil::DeltaR( ijet_p4, jjet_p4 );
            if ( dr < mindr )
            {
                mindr = dr;
                wjetpair = make_pair( ijet, jjet );
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
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 13     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isTightElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 11     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0     )) return false;
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
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 13     )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseElec( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_fo_noiso()[ilep] )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 11        )) return false;
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
float MllSS()
{
    if ( lepidx["TightLepton"].size() != 2 )
        return -999;
    return ( wwwbaby.lep_p4()[jetidx["TightLepton"][0]]
             + wwwbaby.lep_p4()[jetidx["TightLepton"][1]] ).mass();
}

//______________________________________________________________________________________
float Pt3l()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    return ( wwwbaby.lep_p4()[jetidx["TightLepton"][0]]
             + wwwbaby.lep_p4()[jetidx["TightLepton"][1]]
             + wwwbaby.lep_p4()[jetidx["TightLepton"][2]] ).pt();
}

//______________________________________________________________________________________
int LepFlavProduct()
{
    if ( lepidx["GoodSSJet"].size() != 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["GoodSSJet"][0]]
           * wwwbaby.lep_pdgId()[lepidx["GoodSSJet"][1]];
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
    for ( auto& ilep : lepidx["TightLepton"] )
        charge += wwwbaby.lep_charge()[ilep];
    return charge;
}

//______________________________________________________________________________________
int getNumSFOS()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    /* Loops through pairs of entries in the lep_pdgId vector
       and counts how many have opposite value*/
    int num_SFOS = 0;
    for ( int i = 0; i < ( int ) lepidx["TightLepton"].size(); i++ )
    {
        for ( int j = i + 1; j < ( int ) lepidx["TightLepton"].size(); j++ )
        {
            if ( wwwbaby.lep_pdgId()[lepidx["TightLepton"][i]] ==
                    -( wwwbaby.lep_pdgId()[lepidx["TightLepton"][j]] ) )
                num_SFOS++;
        }
    }
    return num_SFOS;
}

//______________________________________________________________________________________
float get0SFOSMll()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][2]];
    if ( pdgid0 == pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get0SFOSMee()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][2]];
    if ( pdgid0 == pdgid1 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass();
    else if ( pdgid0 == pdgid2 && abs( pdgid0 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    else if ( pdgid1 == pdgid2 && abs( pdgid1 ) == 11 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get1SFOSMll()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get2SFOSMll0()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][2]];
    if ( pdgid0 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

//______________________________________________________________________________________
float get2SFOSMll1()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    int pdgid0 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]];
    int pdgid1 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]];
    int pdgid2 = wwwbaby.lep_pdgId()[lepidx["TightLepton"][2]];
    if ( pdgid2 == -pdgid1 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass();
    else if ( pdgid0 == -pdgid2 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    else if ( pdgid1 == -pdgid0 )
        return ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]]
                 + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).mass();
    std::cout <<
        "Warning: Shouldn't be here if function call are at the right places."
        << std::endl;
    return -999;
}

