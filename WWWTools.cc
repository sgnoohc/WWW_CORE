#include "WWWTools.h"

unsigned int objidx_set_to_eventid;
ObjIdx lepidx;
ObjIdx jetidx;

#define MZ 91.1876
#define LUMI 35.9

//______________________________________________________________________________________
bool passSSpresel( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1 ) || ( lepidx["TightLepton"].size()                      ==  2  ) )) { cutidx = -1; return false; }
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==  2  ) )) { cutidx = -2; return false; }
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()   >  30. ) )) { cutidx = -3; return false; }
    if (!( ( Nm1idx == 4 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()   >  30. ) )) { cutidx = -4; return false; }
    if (!( ( Nm1idx == 5 ) || ( jetidx["GoodSSJet"]    .size()                    >=  2  ) )) { cutidx = -5; return false; }
    if (!( ( Nm1idx == 6 ) || ( wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto() ==  0  ) )) { cutidx = -6; return false; }
    if (!( ( Nm1idx == 7 ) || ( isSS()                                                   ) )) { cutidx = -7; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool passSScommon( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSSpresel()                       ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( MjjW()                     > 60.  &&
                                MjjW()                     < 100.    ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( MjjLead()                  < 400.    ) )) { cutidx = -2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( DEtajjLead()               < 1.5     ) )) { cutidx = -3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( jetidx["LooseBJet"].size() == 0      ) )) { cutidx = -4 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool passSSMM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSScommon() ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( isSSMM()       ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( MllSS() > 40.  ) )) { cutidx = -2 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool passSSEM( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSScommon()         ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( isSSEM()               ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( MllSS()          > 30. ) )) { cutidx = -2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.met_pt() > 40. ) )) { cutidx = -3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( getMTmax()       > 90. ) )) { cutidx = -4 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool passSSEE( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( passSScommon()             ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( isSSEE()                   ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( MllSS()          >  40.    ) )) { cutidx = -2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( MllSS()          <  80. ||
                                MllSS()          > 100.    ) )) { cutidx = -3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( wwwbaby.met_pt() >  40.    ) )) { cutidx = -4 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool pass3Lpresel( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 1 ) || ( lepidx["TightLepton"].size()                      ==  3  ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.nlep_VVV_cutbased_veto()                  ==  3  ) )) { cutidx = -2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]].pt()   >  20. ) )) { cutidx = -3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][1]].pt()   >  20. ) )) { cutidx = -4 ; return false; }
    if (!( ( Nm1idx == 5 ) || ( wwwbaby.lep_p4()[lepidx["TightLepton"][2]].pt()   >  20. ) )) { cutidx = -5 ; return false; }
    if (!( ( Nm1idx == 6 ) || ( abs( totalCharge() )                              ==  1  ) )) { cutidx = -6 ; return false; }
    if (!( ( Nm1idx == 7 ) || ( jetidx["Good3LJet"].size()                        <=  1  ) )) { cutidx = -7 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool pass3Lcommon( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( pass3Lpresel()                   ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( Pt3l() > 60.                     ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( jetidx["LooseBJet"].size() ==  0 ) )) { cutidx = -2 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool pass3L0SFOS( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( pass3Lcommon()                      ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( getNumSFOS()               ==   0   ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( get0SFOSMll()              >   20.  ) )) { cutidx = -2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( fabs( get0SFOSMee() - MZ ) >   15.  ) )) { cutidx = -3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( DPhi3lMET()                >    2.7 ) )) { cutidx = -4 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool pass3L1SFOS( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( pass3Lcommon()                      ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( getNumSFOS()               ==   1   ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.met_pt()           >   45.  ) )) { cutidx = -2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( get1SFOSMll()              <   55.||
                                get1SFOSMll()              >  100.  ) )) { cutidx = -3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( DPhi3lMET()                >    2.5 ) )) { cutidx = -4 ; return false; }
    cutidx = 1;
    return true;
}

//______________________________________________________________________________________
bool pass3L2SFOS( int Nm1idx, int& cutidx )
{
    setObjectIndices();
    if (!( ( Nm1idx == 0 ) || ( pass3Lcommon()                       ) )) { cutidx = 0 ; return false; }
    if (!( ( Nm1idx == 1 ) || ( getNumSFOS()                ==   2   ) )) { cutidx = -1 ; return false; }
    if (!( ( Nm1idx == 2 ) || ( wwwbaby.met_pt()            >   55.  ) )) { cutidx = -2 ; return false; }
    if (!( ( Nm1idx == 3 ) || ( fabs( get2SFOSMll0() - MZ ) >   20.  ) )) { cutidx = -3 ; return false; }
    if (!( ( Nm1idx == 4 ) || ( fabs( get2SFOSMll1() - MZ ) >   20.  ) )) { cutidx = -4 ; return false; }
    if (!( ( Nm1idx == 5 ) || ( DPhi3lMET()                 >    2.5 ) )) { cutidx = -5 ; return false; }
    cutidx = 1;
    return true;
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
    std::cout  << "MllSS         " << " : " << MllSS          () << std::endl;
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
    return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["TightLepton"][1]] ).mass();
}

//______________________________________________________________________________________
float Pt3l()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    return ( wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
             + wwwbaby.lep_p4()[lepidx["TightLepton"][1]]
             + wwwbaby.lep_p4()[lepidx["TightLepton"][2]] ).pt();
}

//______________________________________________________________________________________
float DPhi3lMET()
{
    if ( lepidx["TightLepton"].size() != 3 )
        return -999;
    TLorentzVector met;
    met.SetPtEtaPhiM( wwwbaby.met_pt(), 0, wwwbaby.met_phi(), 0 );
    LorentzVector metlv;
    metlv.SetPxPyPzE( met.Px(), met.Py(), met.Pz(), met.E() );
    return fabs( ROOT::Math::VectorUtil::DeltaPhi(
                wwwbaby.lep_p4()[lepidx["TightLepton"][0]]
                + wwwbaby.lep_p4()[lepidx["TightLepton"][1]]
                + wwwbaby.lep_p4()[lepidx["TightLepton"][2]],
                met ));
}

//______________________________________________________________________________________
float getMTmax()
{
    if ( lepidx["TightLepton"].size() < 2 )
        return -999;
    // Set MET
    TLorentzVector met;
    met.SetPtEtaPhiM( wwwbaby.met_pt(), 0, wwwbaby.met_phi(), 0 );
    LorentzVector metlv;
    metlv.SetPxPyPzE( met.Px(), met.Py(), met.Pz(), met.E() );
    // Get Lepton p4
    LorentzVector p4;
    p4 = wwwbaby.lep_p4()[lepidx["TightLepton"][0]];
    float MT0 = sqrt( 2 * p4.pt() * metlv.pt() * ( 1 - cos( ROOT::Math::VectorUtil::DeltaPhi( metlv, p4 ) ) ) );
    p4 = wwwbaby.lep_p4()[lepidx["TightLepton"][1]];
    float MT1 = sqrt( 2 * p4.pt() * metlv.pt() * ( 1 - cos( ROOT::Math::VectorUtil::DeltaPhi( metlv, p4 ) ) ) );
    return MT0 > MT1 ? MT0 : MT1;
}

//______________________________________________________________________________________
int LepFlavProduct()
{
    if ( lepidx["TightLepton"].size() != 2 )
        return false;
    return wwwbaby.lep_pdgId()[lepidx["TightLepton"][0]]
           * wwwbaby.lep_pdgId()[lepidx["TightLepton"][1]];
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

