#include "WWWTools.h"

//______________________________________________________________________________________
bool isTightSSLepton( int ilep )
{
    return isTightSSMuon( ilep ) || isTightSSElectron( ilep );
}

//______________________________________________________________________________________
bool isTightSSMuon( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 13     )) return false;
    if (!(       wwwbaby.lep_p4()[ilep].pt()       >  30.    )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isTightSSElectron( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_tight()[ilep] )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 11     )) return false;
    if (!(       wwwbaby.lep_p4()[ilep].pt()       >  30.    )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5   )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6   )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.06  )) return false;
    if (!( fabs( wwwbaby.lep_ip3d()[ilep]       )  <   0.015 )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0     )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseSSLepton( int ilep )
{
    return isLooseSSMuon( ilep ) || isLooseSSElectron( ilep );
}

//______________________________________________________________________________________
bool isLooseSSMuon( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_fo()[ilep]    )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 13     )) return false;
    if (!(       wwwbaby.lep_p4()[ilep].pt()       >  30.    )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4   )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLooseSSElectron( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_fo_noiso()[ilep] )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 11        )) return false;
    if (!(       wwwbaby.lep_p4()[ilep].pt()       >  30.       )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5      )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6      )) return false;
    if (!(       wwwbaby.lep_3ch_agree()[ilep]     !=  0        )) return false;
    return true;
}

//______________________________________________________________________________________
bool isLbntSSLepton( int ilep )
{
    return isLbntSSMuon( ilep ) || isLbntSSElectron( ilep );
}

//______________________________________________________________________________________
bool isLbntSSMuon( int ilep )
{
    return isLooseSSMuon( ilep ) && !isTightSSMuon( ilep );
}

//______________________________________________________________________________________
bool isLbntSSElectron( int ilep )
{
    return isLooseSSElectron( ilep ) && !isTightSSElectron( ilep );
}

//______________________________________________________________________________________
bool isVetoLepton( int ilep )
{
    return isVetoMuon( ilep ) || isVetoElectron( ilep );
}

//______________________________________________________________________________________
bool isVetoMuon( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_veto_noiso()[ilep] )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 13        )) return false;
    if (!(       wwwbaby.lep_p4()[ilep].pt()       >   5.       )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.4      )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.4      )) return false;
    return true;
}

//______________________________________________________________________________________
bool isVetoElectron( int ilep )
{
    if (!(       wwwbaby.lep_pass_VVV_cutbased_veto_noiso()[ilep] )) return false;
    if (!(       wwwbaby.lep_pdgId()[ilep]         == 11        )) return false;
    if (!(       wwwbaby.lep_p4()[ilep].pt()       >   5.       )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   2.5      )) return false;
    if (!( fabs( wwwbaby.lep_p4()[ilep].eta()   )  <   1.4 ||
           fabs( wwwbaby.lep_p4()[ilep].eta()   )  >   1.6      )) return false;
    if (!(       wwwbaby.lep_relIso03EA()[ilep]    <   0.4      )) return false;
    return true;
}

//______________________________________________________________________________________
bool isGoodSSJet( int ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].pt()      > 30.  )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].eta() )   <  2.5 )) return false;
    return true;
}

//______________________________________________________________________________________
bool isGoodLooseBJet( int ijet )
{
    if (!(       wwwbaby.jets_p4()[ijet].pt()      > 20.     )) return false;
    if (!( fabs( wwwbaby.jets_p4()[ijet].eta() )   <  2.4    )) return false;
    if (!(       wwwbaby.jets_csv()[ijet]          >  0.5426 )) return false;
    return true;
}
