#include <vector>
#include <algorithm>
#include <iterator>

#include <TMath.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMatrixDSym.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>

#include "MvaVariablesEventClassification.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
//#include "TTH/CommonClassifier/interface/DLBDTVars.h"
#include "DLBDTVars.h"
#include "MvaReaderTopJets.h"
#include "MvaVariablesTopJets.h"
#include "higgsUtils.h"
#include "analysisStructs.h"



MvaVariablesEventClassification::MvaVariablesEventClassification():
MvaVariablesBase(),
multiplicity_jets_(MvaVariableInt(name_multiplicity_jets_)),
btagDiscriminatorAverage_tagged_(MvaVariableFloat(name_btagDiscriminatorAverage_tagged_)),
btagDiscriminatorAverage_untagged_(MvaVariableFloat(name_btagDiscriminatorAverage_untagged_)),
minDeltaR_jet_jet_(MvaVariableFloat(name_minDeltaR_jet_jet_)),
minDeltaR_tag_tag_(MvaVariableFloat(name_minDeltaR_tag_tag_)),
avgDeltaR_jet_jet_(MvaVariableFloat(name_avgDeltaR_jet_jet_)),
avgDeltaR_jet_tag_(MvaVariableFloat(name_avgDeltaR_jet_tag_)),
avgDeltaR_tag_tag_(MvaVariableFloat(name_avgDeltaR_tag_tag_)),
ptSum_jets_leptons_(MvaVariableFloat(name_ptSum_jets_leptons_)),
multiplicity_higgsLikeDijet15_(MvaVariableInt(name_multiplicity_higgsLikeDijet15_)),
mass_higgsLikeDijet_(MvaVariableFloat(name_mass_higgsLikeDijet_)),
mass_higgsLikeDijet2_(MvaVariableFloat(name_mass_higgsLikeDijet2_)),
mass_jet_jet_min_deltaR_(MvaVariableFloat(name_mass_jet_jet_min_deltaR_)),
mass_tag_tag_min_deltaR_(MvaVariableFloat(name_mass_tag_tag_min_deltaR_)),
mass_jet_tag_min_deltaR_(MvaVariableFloat(name_mass_jet_tag_min_deltaR_)),
mass_tag_tag_max_mass_(MvaVariableFloat(name_mass_tag_tag_max_mass_)),
median_mass_jet_jet_(MvaVariableFloat(name_median_mass_jet_jet_)),
maxDeltaEta_jet_jet_(MvaVariableFloat(name_maxDeltaEta_jet_jet_)),
maxDeltaEta_tag_tag_(MvaVariableFloat(name_maxDeltaEta_tag_tag_)),
HT_jets_(MvaVariableFloat(name_HT_jets_)),
HT_tags_(MvaVariableFloat(name_HT_tags_)),
pT_jet_jet_min_deltaR_(MvaVariableFloat(name_pT_jet_jet_min_deltaR_)),
pT_jet_tag_min_deltaR_(MvaVariableFloat(name_pT_jet_tag_min_deltaR_)),
pT_tag_tag_min_deltaR_(MvaVariableFloat(name_pT_tag_tag_min_deltaR_)),
mass_jet_jet_jet_max_pT_(MvaVariableFloat(name_mass_jet_jet_jet_max_pT_)),
mass_jet_tag_tag_max_pT_(MvaVariableFloat(name_mass_jet_tag_tag_max_pT_)),
centrality_jets_leps_(MvaVariableFloat(name_centrality_jets_leps_)),
centrality_tags_(MvaVariableFloat(name_centrality_tags_)),
twist_jet_jet_max_mass_(MvaVariableFloat(name_twist_jet_jet_max_mass_)),
twist_jet_tag_max_mass_(MvaVariableFloat(name_twist_jet_tag_max_mass_)),
twist_tag_tag_max_mass_(MvaVariableFloat(name_twist_tag_tag_max_mass_)),
twist_tag_tag_min_deltaR_(MvaVariableFloat(name_twist_tag_tag_min_deltaR_)),
sphericity_jet_(MvaVariableFloat(name_sphericity_jet_)),
aplanarity_jet_(MvaVariableFloat(name_aplanarity_jet_)),
circularity_jet_(MvaVariableFloat(name_circularity_jet_)),
isotropy_jet_(MvaVariableFloat(name_isotropy_jet_)),
C_jet_(MvaVariableFloat(name_C_jet_)),
D_jet_(MvaVariableFloat(name_D_jet_)),
transSphericity_jet_(MvaVariableFloat(name_transSphericity_jet_)),
sphericity_tag_(MvaVariableFloat(name_sphericity_tag_)),
aplanarity_tag_(MvaVariableFloat(name_aplanarity_tag_)),
circularity_tag_(MvaVariableFloat(name_circularity_tag_)),
isotropy_tag_(MvaVariableFloat(name_isotropy_tag_)),
C_tag_(MvaVariableFloat(name_C_tag_)),
D_tag_(MvaVariableFloat(name_D_tag_)),
transSphericity_tag_(MvaVariableFloat(name_transSphericity_tag_)),
H0_jet_(MvaVariableFloat(name_H0_jet_)),
H1_jet_(MvaVariableFloat(name_H1_jet_)),
H2_jet_(MvaVariableFloat(name_H2_jet_)),
H3_jet_(MvaVariableFloat(name_H3_jet_)),
H4_jet_(MvaVariableFloat(name_H4_jet_)),
R1_jet_(MvaVariableFloat(name_R1_jet_)),
R2_jet_(MvaVariableFloat(name_R2_jet_)),
R3_jet_(MvaVariableFloat(name_R3_jet_)),
R4_jet_(MvaVariableFloat(name_R4_jet_)),
H0_tag_(MvaVariableFloat(name_H0_tag_)),
H1_tag_(MvaVariableFloat(name_H1_tag_)),
H2_tag_(MvaVariableFloat(name_H2_tag_)),
H3_tag_(MvaVariableFloat(name_H3_tag_)),
H4_tag_(MvaVariableFloat(name_H4_tag_)),
R1_tag_(MvaVariableFloat(name_R1_tag_)),
R2_tag_(MvaVariableFloat(name_R2_tag_)),
R3_tag_(MvaVariableFloat(name_R3_tag_)),
R4_tag_(MvaVariableFloat(name_R4_tag_)),
mass_bb_(MvaVariableFloat(name_mass_bb_))
{}



MvaVariablesEventClassification::MvaVariablesEventClassification(
    const EventMetadata& eventMetadata,
    const int multiplicity_jets,
    const double& btagDiscriminatorAverage_tagged, const double& btagDiscriminatorAverage_untagged,
    const double& minDeltaR_jet_jet, const double& minDeltaR_tag_tag,
    const double& avgDeltaR_jet_jet, const double& avgDeltaR_jet_tag, const double& avgDeltaR_tag_tag,
    const double& ptSum_jets_leptons,
    const int multiplicity_higgsLikeDijet15,
    const double& mass_higgsLikeDijet, const double& mass_higgsLikeDijet2,
    const double& mass_jet_jet_min_deltaR, const double& mass_tag_tag_min_deltaR, const double& mass_jet_tag_min_deltaR,
    const double& mass_tag_tag_max_mass,
    const double& median_mass_jet_jet,
    const double& maxDeltaEta_jet_jet, const double& maxDeltaEta_tag_tag,
    const double& HT_jets, const double& HT_tags,
    const double& pT_jet_jet_min_deltaR, const double& pT_jet_tag_min_deltaR, const double& pT_tag_tag_min_deltaR,
    const double& mass_jet_jet_jet_max_pT, const double& mass_jet_tag_tag_max_pT,
    const double& centrality_jets_leps, const double& centrality_tags,
    const double& twist_jet_jet_max_mass, const double& twist_jet_tag_max_mass, const double& twist_tag_tag_max_mass,
    const double& twist_tag_tag_min_deltaR,
    const double& sphericity_jet, const double& aplanarity_jet, const double& circularity_jet, const double& isotropy_jet,
    const double& C_jet, const double& D_jet, const double& transSphericity_jet,
    const double& sphericity_tag, const double& aplanarity_tag, const double& circularity_tag, const double& isotropy_tag,
    const double& C_tag, const double& D_tag, const double& transSphericity_tag,
    const double& H0_jet, const double& H1_jet, const double& H2_jet, const double& H3_jet, const double& H4_jet,
    const double& R1_jet, const double& R2_jet, const double& R3_jet, const double& R4_jet,
    const double& H0_tag, const double& H1_tag, const double& H2_tag, const double& H3_tag, const double& H4_tag,
    const double& R1_tag, const double& R2_tag, const double& R3_tag, const double& R4_tag,
    const double& mass_bb,
    const double& eventWeight):
MvaVariablesBase(eventMetadata, eventWeight),
multiplicity_jets_(MvaVariableInt(name_multiplicity_jets_)),
btagDiscriminatorAverage_tagged_(MvaVariableFloat(name_btagDiscriminatorAverage_tagged_)),
btagDiscriminatorAverage_untagged_(MvaVariableFloat(name_btagDiscriminatorAverage_untagged_)),
minDeltaR_jet_jet_(MvaVariableFloat(name_minDeltaR_jet_jet_)),
minDeltaR_tag_tag_(MvaVariableFloat(name_minDeltaR_tag_tag_)),
avgDeltaR_jet_jet_(MvaVariableFloat(name_avgDeltaR_jet_jet_)),
avgDeltaR_jet_tag_(MvaVariableFloat(name_avgDeltaR_jet_tag_)),
avgDeltaR_tag_tag_(MvaVariableFloat(name_avgDeltaR_tag_tag_)),
ptSum_jets_leptons_(MvaVariableFloat(name_ptSum_jets_leptons_)),
multiplicity_higgsLikeDijet15_(MvaVariableInt(name_multiplicity_higgsLikeDijet15_)),
mass_higgsLikeDijet_(MvaVariableFloat(name_mass_higgsLikeDijet_)),
mass_higgsLikeDijet2_(MvaVariableFloat(name_mass_higgsLikeDijet2_)),
mass_jet_jet_min_deltaR_(MvaVariableFloat(name_mass_jet_jet_min_deltaR_)),
mass_tag_tag_min_deltaR_(MvaVariableFloat(name_mass_tag_tag_min_deltaR_)),
mass_jet_tag_min_deltaR_(MvaVariableFloat(name_mass_jet_tag_min_deltaR_)),
mass_tag_tag_max_mass_(MvaVariableFloat(name_mass_tag_tag_max_mass_)),
median_mass_jet_jet_(MvaVariableFloat(name_median_mass_jet_jet_)),
maxDeltaEta_jet_jet_(MvaVariableFloat(name_maxDeltaEta_jet_jet_)),
maxDeltaEta_tag_tag_(MvaVariableFloat(name_maxDeltaEta_tag_tag_)),
HT_jets_(MvaVariableFloat(name_HT_jets_)),
HT_tags_(MvaVariableFloat(name_HT_tags_)),
pT_jet_jet_min_deltaR_(MvaVariableFloat(name_pT_jet_jet_min_deltaR_)),
pT_jet_tag_min_deltaR_(MvaVariableFloat(name_pT_jet_tag_min_deltaR_)),
pT_tag_tag_min_deltaR_(MvaVariableFloat(name_pT_tag_tag_min_deltaR_)),
mass_jet_jet_jet_max_pT_(MvaVariableFloat(name_mass_jet_jet_jet_max_pT_)),
mass_jet_tag_tag_max_pT_(MvaVariableFloat(name_mass_jet_tag_tag_max_pT_)),
centrality_jets_leps_(MvaVariableFloat(name_centrality_jets_leps_)),
centrality_tags_(MvaVariableFloat(name_centrality_tags_)),
twist_jet_jet_max_mass_(MvaVariableFloat(name_twist_jet_jet_max_mass_)),
twist_jet_tag_max_mass_(MvaVariableFloat(name_twist_jet_tag_max_mass_)),
twist_tag_tag_max_mass_(MvaVariableFloat(name_twist_tag_tag_max_mass_)),
twist_tag_tag_min_deltaR_(MvaVariableFloat(name_twist_tag_tag_min_deltaR_)),
sphericity_jet_(MvaVariableFloat(name_sphericity_jet_)),
aplanarity_jet_(MvaVariableFloat(name_aplanarity_jet_)),
circularity_jet_(MvaVariableFloat(name_circularity_jet_)),
isotropy_jet_(MvaVariableFloat(name_isotropy_jet_)),
C_jet_(MvaVariableFloat(name_C_jet_)),
D_jet_(MvaVariableFloat(name_D_jet_)),
transSphericity_jet_(MvaVariableFloat(name_transSphericity_jet_)),
sphericity_tag_(MvaVariableFloat(name_sphericity_tag_)),
aplanarity_tag_(MvaVariableFloat(name_aplanarity_tag_)),
circularity_tag_(MvaVariableFloat(name_circularity_tag_)),
isotropy_tag_(MvaVariableFloat(name_isotropy_tag_)),
C_tag_(MvaVariableFloat(name_C_tag_)),
D_tag_(MvaVariableFloat(name_D_tag_)),
transSphericity_tag_(MvaVariableFloat(name_transSphericity_tag_)),
H0_jet_(MvaVariableFloat(name_H0_jet_)),
H1_jet_(MvaVariableFloat(name_H1_jet_)),
H2_jet_(MvaVariableFloat(name_H2_jet_)),
H3_jet_(MvaVariableFloat(name_H3_jet_)),
H4_jet_(MvaVariableFloat(name_H4_jet_)),
R1_jet_(MvaVariableFloat(name_R1_jet_)),
R2_jet_(MvaVariableFloat(name_R2_jet_)),
R3_jet_(MvaVariableFloat(name_R3_jet_)),
R4_jet_(MvaVariableFloat(name_R4_jet_)),
H0_tag_(MvaVariableFloat(name_H0_tag_)),
H1_tag_(MvaVariableFloat(name_H1_tag_)),
H2_tag_(MvaVariableFloat(name_H2_tag_)),
H3_tag_(MvaVariableFloat(name_H3_tag_)),
H4_tag_(MvaVariableFloat(name_H4_tag_)),
R1_tag_(MvaVariableFloat(name_R1_tag_)),
R2_tag_(MvaVariableFloat(name_R2_tag_)),
R3_tag_(MvaVariableFloat(name_R3_tag_)),
R4_tag_(MvaVariableFloat(name_R4_tag_)),
mass_bb_(MvaVariableFloat(name_mass_bb_))
{
    // Fill the variables for MVA TTree
    multiplicity_jets_.setValue(multiplicity_jets);
    btagDiscriminatorAverage_tagged_.setValue(btagDiscriminatorAverage_tagged);
    btagDiscriminatorAverage_untagged_.setValue(btagDiscriminatorAverage_untagged);
    minDeltaR_jet_jet_.setValue(minDeltaR_jet_jet);
    minDeltaR_tag_tag_.setValue(minDeltaR_tag_tag);
    avgDeltaR_jet_jet_.setValue(avgDeltaR_jet_jet);
    avgDeltaR_jet_tag_.setValue(avgDeltaR_jet_tag);
    avgDeltaR_tag_tag_.setValue(avgDeltaR_tag_tag);
    ptSum_jets_leptons_.setValue(ptSum_jets_leptons);
    multiplicity_higgsLikeDijet15_.setValue(multiplicity_higgsLikeDijet15);
    mass_higgsLikeDijet_.setValue(mass_higgsLikeDijet);
    mass_higgsLikeDijet2_.setValue(mass_higgsLikeDijet2);
    mass_jet_jet_min_deltaR_.setValue(mass_jet_jet_min_deltaR);
    mass_jet_tag_min_deltaR_.setValue(mass_jet_tag_min_deltaR);
    mass_tag_tag_min_deltaR_.setValue(mass_tag_tag_min_deltaR);
    mass_tag_tag_max_mass_.setValue(mass_tag_tag_max_mass);
    median_mass_jet_jet_.setValue(median_mass_jet_jet);
    maxDeltaEta_jet_jet_.setValue(maxDeltaEta_jet_jet);
    maxDeltaEta_tag_tag_.setValue(maxDeltaEta_tag_tag);
    HT_jets_.setValue(HT_jets);
    HT_tags_.setValue(HT_tags);
    pT_jet_jet_min_deltaR_.setValue(pT_jet_jet_min_deltaR);
    pT_jet_tag_min_deltaR_.setValue(pT_jet_tag_min_deltaR);
    pT_tag_tag_min_deltaR_.setValue(pT_tag_tag_min_deltaR);
    mass_jet_jet_jet_max_pT_.setValue(mass_jet_jet_jet_max_pT);
    mass_jet_tag_tag_max_pT_.setValue(mass_jet_tag_tag_max_pT);
    centrality_jets_leps_.setValue(centrality_jets_leps);
    centrality_tags_.setValue(centrality_tags);
    twist_jet_jet_max_mass_.setValue(twist_jet_jet_max_mass);
    twist_jet_tag_max_mass_.setValue(twist_jet_tag_max_mass);
    twist_tag_tag_max_mass_.setValue(twist_tag_tag_max_mass);
    twist_tag_tag_min_deltaR_.setValue(twist_tag_tag_min_deltaR);

    // Event shape variables
    sphericity_jet_.setValue(sphericity_jet);
    aplanarity_jet_.setValue(aplanarity_jet);
    circularity_jet_.setValue(circularity_jet);
    isotropy_jet_.setValue(isotropy_jet);
    C_jet_.setValue(C_jet);
    D_jet_.setValue(D_jet);
    transSphericity_jet_.setValue(transSphericity_jet);

    sphericity_tag_.setValue(sphericity_tag);
    aplanarity_tag_.setValue(aplanarity_tag);
    circularity_tag_.setValue(circularity_tag);
    isotropy_tag_.setValue(isotropy_tag);
    C_tag_.setValue(C_tag);
    D_tag_.setValue(D_tag);
    transSphericity_tag_.setValue(transSphericity_tag);

    H0_jet_.setValue(H0_jet);
    H1_jet_.setValue(H1_jet);
    H2_jet_.setValue(H2_jet);
    H3_jet_.setValue(H3_jet);
    H4_jet_.setValue(H4_jet);
    R1_jet_.setValue(R1_jet);
    R2_jet_.setValue(R2_jet);
    R3_jet_.setValue(R3_jet);
    R4_jet_.setValue(R4_jet);

    H0_tag_.setValue(H0_tag);
    H1_tag_.setValue(H1_tag);
    H2_tag_.setValue(H2_tag);
    H3_tag_.setValue(H3_tag);
    H4_tag_.setValue(H4_tag);
    R1_tag_.setValue(R1_tag);
    R2_tag_.setValue(R2_tag);
    R3_tag_.setValue(R3_tag);
    R4_tag_.setValue(R4_tag);

    mass_bb_.setValue(mass_bb);
}



MvaVariablesEventClassification::MvaVariablesEventClassification(const EventMetadata& eventMetadata, const double& eventWeight,
                                                                 const DLBDTMvaVariablesEventClassification& dlBdtMvaVariables):
MvaVariablesBase(eventMetadata, eventWeight),
multiplicity_jets_(MvaVariableInt(dlBdtMvaVariables.multiplicity_jets_.name().data())),
btagDiscriminatorAverage_tagged_(MvaVariableFloat(dlBdtMvaVariables.btagDiscriminatorAverage_tagged_.name().data())),
btagDiscriminatorAverage_untagged_(MvaVariableFloat(dlBdtMvaVariables.btagDiscriminatorAverage_untagged_.name().data())),
minDeltaR_jet_jet_(MvaVariableFloat(dlBdtMvaVariables.minDeltaR_jet_jet_.name().data())),
minDeltaR_tag_tag_(MvaVariableFloat(dlBdtMvaVariables.minDeltaR_tag_tag_.name().data())),
avgDeltaR_jet_jet_(MvaVariableFloat(dlBdtMvaVariables.avgDeltaR_jet_jet_.name().data())),
avgDeltaR_jet_tag_(MvaVariableFloat(dlBdtMvaVariables.avgDeltaR_jet_tag_.name().data())),
avgDeltaR_tag_tag_(MvaVariableFloat(dlBdtMvaVariables.avgDeltaR_tag_tag_.name().data())),
ptSum_jets_leptons_(MvaVariableFloat(dlBdtMvaVariables.ptSum_jets_leptons_.name().data())),
multiplicity_higgsLikeDijet15_(MvaVariableInt(dlBdtMvaVariables.multiplicity_higgsLikeDijet15_.name().data())),
mass_higgsLikeDijet_(MvaVariableFloat(dlBdtMvaVariables.mass_higgsLikeDijet_.name().data())),
mass_higgsLikeDijet2_(MvaVariableFloat(dlBdtMvaVariables.mass_higgsLikeDijet2_.name().data())),
mass_jet_jet_min_deltaR_(MvaVariableFloat(dlBdtMvaVariables.mass_jet_jet_min_deltaR_.name().data())),
mass_tag_tag_min_deltaR_(MvaVariableFloat(dlBdtMvaVariables.mass_tag_tag_min_deltaR_.name().data())),
mass_jet_tag_min_deltaR_(MvaVariableFloat(dlBdtMvaVariables.mass_jet_tag_min_deltaR_.name().data())),
mass_tag_tag_max_mass_(MvaVariableFloat(dlBdtMvaVariables.mass_tag_tag_max_mass_.name().data())),
median_mass_jet_jet_(MvaVariableFloat(dlBdtMvaVariables.median_mass_jet_jet_.name().data())),
maxDeltaEta_jet_jet_(MvaVariableFloat(dlBdtMvaVariables.maxDeltaEta_jet_jet_.name().data())),
maxDeltaEta_tag_tag_(MvaVariableFloat(dlBdtMvaVariables.maxDeltaEta_tag_tag_.name().data())),
HT_jets_(MvaVariableFloat(dlBdtMvaVariables.HT_jets_.name().data())),
HT_tags_(MvaVariableFloat(dlBdtMvaVariables.HT_tags_.name().data())),
pT_jet_jet_min_deltaR_(MvaVariableFloat(dlBdtMvaVariables.pT_jet_jet_min_deltaR_.name().data())),
pT_jet_tag_min_deltaR_(MvaVariableFloat(dlBdtMvaVariables.pT_jet_tag_min_deltaR_.name().data())),
pT_tag_tag_min_deltaR_(MvaVariableFloat(dlBdtMvaVariables.pT_tag_tag_min_deltaR_.name().data())),
mass_jet_jet_jet_max_pT_(MvaVariableFloat(dlBdtMvaVariables.mass_jet_jet_jet_max_pT_.name().data())),
mass_jet_tag_tag_max_pT_(MvaVariableFloat(dlBdtMvaVariables.mass_jet_tag_tag_max_pT_.name().data())),
centrality_jets_leps_(MvaVariableFloat(dlBdtMvaVariables.centrality_jets_leps_.name().data())),
centrality_tags_(MvaVariableFloat(dlBdtMvaVariables.centrality_tags_.name().data())),
twist_jet_jet_max_mass_(MvaVariableFloat(dlBdtMvaVariables.twist_jet_jet_max_mass_.name().data())),
twist_jet_tag_max_mass_(MvaVariableFloat(dlBdtMvaVariables.twist_jet_tag_max_mass_.name().data())),
twist_tag_tag_max_mass_(MvaVariableFloat(dlBdtMvaVariables.twist_tag_tag_max_mass_.name().data())),
twist_tag_tag_min_deltaR_(MvaVariableFloat(dlBdtMvaVariables.twist_tag_tag_min_deltaR_.name().data())),
sphericity_jet_(MvaVariableFloat(dlBdtMvaVariables.sphericity_jet_.name().data())),
aplanarity_jet_(MvaVariableFloat(dlBdtMvaVariables.aplanarity_jet_.name().data())),
circularity_jet_(MvaVariableFloat(dlBdtMvaVariables.circularity_jet_.name().data())),
isotropy_jet_(MvaVariableFloat(dlBdtMvaVariables.isotropy_jet_.name().data())),
C_jet_(MvaVariableFloat(dlBdtMvaVariables.C_jet_.name().data())),
D_jet_(MvaVariableFloat(dlBdtMvaVariables.D_jet_.name().data())),
transSphericity_jet_(MvaVariableFloat(dlBdtMvaVariables.transSphericity_jet_.name().data())),
sphericity_tag_(MvaVariableFloat(dlBdtMvaVariables.sphericity_tag_.name().data())),
aplanarity_tag_(MvaVariableFloat(dlBdtMvaVariables.aplanarity_tag_.name().data())),
circularity_tag_(MvaVariableFloat(dlBdtMvaVariables.circularity_tag_.name().data())),
isotropy_tag_(MvaVariableFloat(dlBdtMvaVariables.isotropy_tag_.name().data())),
C_tag_(MvaVariableFloat(dlBdtMvaVariables.C_tag_.name().data())),
D_tag_(MvaVariableFloat(dlBdtMvaVariables.D_tag_.name().data())),
transSphericity_tag_(MvaVariableFloat(dlBdtMvaVariables.transSphericity_tag_.name().data())),
H0_jet_(MvaVariableFloat(dlBdtMvaVariables.H0_jet_.name().data())),
H1_jet_(MvaVariableFloat(dlBdtMvaVariables.H1_jet_.name().data())),
H2_jet_(MvaVariableFloat(dlBdtMvaVariables.H2_jet_.name().data())),
H3_jet_(MvaVariableFloat(dlBdtMvaVariables.H3_jet_.name().data())),
H4_jet_(MvaVariableFloat(dlBdtMvaVariables.H4_jet_.name().data())),
R1_jet_(MvaVariableFloat(dlBdtMvaVariables.R1_jet_.name().data())),
R2_jet_(MvaVariableFloat(dlBdtMvaVariables.R2_jet_.name().data())),
R3_jet_(MvaVariableFloat(dlBdtMvaVariables.R3_jet_.name().data())),
R4_jet_(MvaVariableFloat(dlBdtMvaVariables.R4_jet_.name().data())),
H0_tag_(MvaVariableFloat(dlBdtMvaVariables.H0_tag_.name().data())),
H1_tag_(MvaVariableFloat(dlBdtMvaVariables.H1_tag_.name().data())),
H2_tag_(MvaVariableFloat(dlBdtMvaVariables.H2_tag_.name().data())),
H3_tag_(MvaVariableFloat(dlBdtMvaVariables.H3_tag_.name().data())),
H4_tag_(MvaVariableFloat(dlBdtMvaVariables.H4_tag_.name().data())),
R1_tag_(MvaVariableFloat(dlBdtMvaVariables.R1_tag_.name().data())),
R2_tag_(MvaVariableFloat(dlBdtMvaVariables.R2_tag_.name().data())),
R3_tag_(MvaVariableFloat(dlBdtMvaVariables.R3_tag_.name().data())),
R4_tag_(MvaVariableFloat(dlBdtMvaVariables.R4_tag_.name().data())),
mass_bb_(MvaVariableFloat(dlBdtMvaVariables.mass_bb_.name().data()))
{
    // Fill the variables for MVA TTree
//    eventNumber_.setValue(dlBdtMvaVariables.eventNumber_.value_);
    multiplicity_jets_.setValue(dlBdtMvaVariables.multiplicity_jets_.value_);
    btagDiscriminatorAverage_tagged_.setValue(dlBdtMvaVariables.btagDiscriminatorAverage_tagged_.value_);
    btagDiscriminatorAverage_untagged_.setValue(dlBdtMvaVariables.btagDiscriminatorAverage_untagged_.value_);
    minDeltaR_jet_jet_.setValue(dlBdtMvaVariables.minDeltaR_jet_jet_.value_);
    minDeltaR_tag_tag_.setValue(dlBdtMvaVariables.minDeltaR_tag_tag_.value_);
    avgDeltaR_jet_jet_.setValue(dlBdtMvaVariables.avgDeltaR_jet_jet_.value_);
    avgDeltaR_jet_tag_.setValue(dlBdtMvaVariables.avgDeltaR_jet_tag_.value_);
    avgDeltaR_tag_tag_.setValue(dlBdtMvaVariables.avgDeltaR_tag_tag_.value_);
    ptSum_jets_leptons_.setValue(dlBdtMvaVariables.ptSum_jets_leptons_.value_);
    multiplicity_higgsLikeDijet15_.setValue(dlBdtMvaVariables.multiplicity_higgsLikeDijet15_.value_);
    mass_higgsLikeDijet_.setValue(dlBdtMvaVariables.mass_higgsLikeDijet_.value_);
    mass_higgsLikeDijet2_.setValue(dlBdtMvaVariables.mass_higgsLikeDijet2_.value_);
    mass_jet_jet_min_deltaR_.setValue(dlBdtMvaVariables.mass_jet_jet_min_deltaR_.value_);
    mass_jet_tag_min_deltaR_.setValue(dlBdtMvaVariables.mass_jet_tag_min_deltaR_.value_);
    mass_tag_tag_min_deltaR_.setValue(dlBdtMvaVariables.mass_tag_tag_min_deltaR_.value_);
    mass_tag_tag_max_mass_.setValue(dlBdtMvaVariables.mass_tag_tag_max_mass_.value_);
    median_mass_jet_jet_.setValue(dlBdtMvaVariables.median_mass_jet_jet_.value_);
    maxDeltaEta_jet_jet_.setValue(dlBdtMvaVariables.maxDeltaEta_jet_jet_.value_);
    maxDeltaEta_tag_tag_.setValue(dlBdtMvaVariables.maxDeltaEta_tag_tag_.value_);
    HT_jets_.setValue(dlBdtMvaVariables.HT_jets_.value_);
    HT_tags_.setValue(dlBdtMvaVariables.HT_tags_.value_);
    pT_jet_jet_min_deltaR_.setValue(dlBdtMvaVariables.pT_jet_jet_min_deltaR_.value_);
    pT_jet_tag_min_deltaR_.setValue(dlBdtMvaVariables.pT_jet_tag_min_deltaR_.value_);
    pT_tag_tag_min_deltaR_.setValue(dlBdtMvaVariables.pT_tag_tag_min_deltaR_.value_);
    mass_jet_jet_jet_max_pT_.setValue(dlBdtMvaVariables.mass_jet_jet_jet_max_pT_.value_);
    mass_jet_tag_tag_max_pT_.setValue(dlBdtMvaVariables.mass_jet_tag_tag_max_pT_.value_);
    centrality_jets_leps_.setValue(dlBdtMvaVariables.centrality_jets_leps_.value_);
    centrality_tags_.setValue(dlBdtMvaVariables.centrality_tags_.value_);
    twist_jet_jet_max_mass_.setValue(dlBdtMvaVariables.twist_jet_jet_max_mass_.value_);
    twist_jet_tag_max_mass_.setValue(dlBdtMvaVariables.twist_jet_tag_max_mass_.value_);
    twist_tag_tag_max_mass_.setValue(dlBdtMvaVariables.twist_tag_tag_max_mass_.value_);
    twist_tag_tag_min_deltaR_.setValue(dlBdtMvaVariables.twist_tag_tag_min_deltaR_.value_);

    // Event shape variables
    sphericity_jet_.setValue(dlBdtMvaVariables.sphericity_jet_.value_);
    aplanarity_jet_.setValue(dlBdtMvaVariables.aplanarity_jet_.value_);
    circularity_jet_.setValue(dlBdtMvaVariables.circularity_jet_.value_);
    isotropy_jet_.setValue(dlBdtMvaVariables.isotropy_jet_.value_);
    C_jet_.setValue(dlBdtMvaVariables.C_jet_.value_);
    D_jet_.setValue(dlBdtMvaVariables.D_jet_.value_);
    transSphericity_jet_.setValue(dlBdtMvaVariables.transSphericity_jet_.value_);

    sphericity_tag_.setValue(dlBdtMvaVariables.sphericity_tag_.value_);
    aplanarity_tag_.setValue(dlBdtMvaVariables.aplanarity_tag_.value_);
    circularity_tag_.setValue(dlBdtMvaVariables.circularity_tag_.value_);
    isotropy_tag_.setValue(dlBdtMvaVariables.isotropy_tag_.value_);
    C_tag_.setValue(dlBdtMvaVariables.C_tag_.value_);
    D_tag_.setValue(dlBdtMvaVariables.D_tag_.value_);
    transSphericity_tag_.setValue(dlBdtMvaVariables.transSphericity_tag_.value_);

    H0_jet_.setValue(dlBdtMvaVariables.H0_jet_.value_);
    H1_jet_.setValue(dlBdtMvaVariables.H1_jet_.value_);
    H2_jet_.setValue(dlBdtMvaVariables.H2_jet_.value_);
    H3_jet_.setValue(dlBdtMvaVariables.H3_jet_.value_);
    H4_jet_.setValue(dlBdtMvaVariables.H4_jet_.value_);
    R1_jet_.setValue(dlBdtMvaVariables.R1_jet_.value_);
    R2_jet_.setValue(dlBdtMvaVariables.R2_jet_.value_);
    R3_jet_.setValue(dlBdtMvaVariables.R3_jet_.value_);
    R4_jet_.setValue(dlBdtMvaVariables.R4_jet_.value_);

    H0_tag_.setValue(dlBdtMvaVariables.H0_tag_.value_);
    H1_tag_.setValue(dlBdtMvaVariables.H1_tag_.value_);
    H2_tag_.setValue(dlBdtMvaVariables.H2_tag_.value_);
    H3_tag_.setValue(dlBdtMvaVariables.H3_tag_.value_);
    H4_tag_.setValue(dlBdtMvaVariables.H4_tag_.value_);
    R1_tag_.setValue(dlBdtMvaVariables.R1_tag_.value_);
    R2_tag_.setValue(dlBdtMvaVariables.R2_tag_.value_);
    R3_tag_.setValue(dlBdtMvaVariables.R3_tag_.value_);
    R4_tag_.setValue(dlBdtMvaVariables.R4_tag_.value_);

    mass_bb_.setValue(dlBdtMvaVariables.mass_bb_.value_);
}

MvaVariablesEventClassification:: TopPairVariable:: TopPairVariable(std::string weight_topSystemBDT) : topSystemWeight_(new MvaReaderTopJets("BDT top system identification")), weight_topSystemBDT_(weight_topSystemBDT){
    topSystemWeight_->book(weight_topSystemBDT_);
}

std::pair<int,int> MvaVariablesEventClassification::TopPairVariable::jetPairsFromMVA(const EventMetadata& eventMetadata,
                                                                                       const tth::RecoObjectIndices& recoObjectIndices,
                                                                                       const tth::GenObjectIndices& genObjectIndices,
                                                                                       const RecoObjects& recoObjects,
                                                                                       const double weight) {
    
    //Setting up the MVA input #################################
    //Loop over all jet combinations and get MVA input variables
    std::vector<MvaVariablesBase*> v_mvaVariables = MvaVariablesTopJets::fillVariables(eventMetadata, recoObjectIndices, genObjectIndices, recoObjects, weight);
    // Getting the MVA weights from weights file as vector, one entry per jet pair
    std::vector<float> v_mvaWeightsCorrect;
    if(topSystemWeight_) v_mvaWeightsCorrect = topSystemWeight_->mvaWeights(v_mvaVariables);
    MvaVariablesTopJets::clearVariables(v_mvaVariables);
    
    std::pair<int, int> topPair;

    if(!topSystemWeight_) return topPair;
    
    double maxWeight = -999.;
    int maxIndex = -1;
    
    // Loop over jet pairs
    const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
    for(size_t i=0; i<jetIndexPairs.size(); ++i) {
        if (v_mvaWeightsCorrect.at(i) > maxWeight) {
            maxWeight = v_mvaWeightsCorrect.at(i);
            maxIndex = i;
            topPair = jetIndexPairs.at(maxIndex);
        }
    } 
   return topPair;
}

MvaVariablesEventClassification* MvaVariablesEventClassification::fillVariables(const EventMetadata& eventMetadata,
                                                                                const tth::RecoObjectIndices& recoObjectIndices,
                                                                                const RecoObjects& recoObjects,
                                                                                const tth::GenObjectIndices& genObjectIndices,
                                                                                const double& eventWeight)
{
    using ROOT::Math::VectorUtil::DeltaR;
    
    // Access relevant objects and indices
    const std::vector<double>& jetBtags(*recoObjects.jetBtags_);
    const VLV& leptons(*recoObjects.allLeptons_);
    const VLV& jets(*recoObjects.jets_);
    
    // Identify the most likely pair to stem from tt 
    std::pair<int,int> topPair;
    
    std::string weight_topSystemBDT = "/nfs/dust/cms/user/chrisjcc/Variables2/Area1/analysisArea/xmlWeights/correct_step7_cate5_cate6_cate7_c1.weights.xml"; 
    
    if(recoObjectIndices.jetIndices_.size()>1) topPair = TopPairVariable::Instance(weight_topSystemBDT).jetPairsFromMVA(eventMetadata, recoObjectIndices, genObjectIndices, recoObjects, eventWeight); 
    
    // Calculate several jet-dependent quantities
    double btagDiscriminatorSumTagged(0.);
    double btagDiscriminatorSumUntagged(0.);
    double minDeltaRJetJet(999.);
    double ptSumJets(0.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
        // Calculate the btag-discriminator averages, setting values<0. to 0.
        const double& btagDiscriminator(jetBtags.at(*i_index));
        const double btagDiscriminatorPositive(btagDiscriminator>=0. ? btagDiscriminator : 0.);
        // Avoid b-tag values where the algorithm did not work, giving values>1., setting them to 1.
        const double btagDiscriminatorInRange(btagDiscriminatorPositive<=1. ? btagDiscriminatorPositive : 1.);
        if(std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
            btagDiscriminatorSumTagged += btagDiscriminatorInRange;
        }
        else{
            btagDiscriminatorSumUntagged += btagDiscriminatorInRange;
        }
        
        // Minimum deltaR between jets
        for(auto j_index = i_index+1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
            const double deltaR = DeltaR(jets.at(*i_index), jets.at(*j_index));
            if(deltaR < minDeltaRJetJet) minDeltaRJetJet = deltaR;
        }
        
        // Scalar sum of pt of all jets
        ptSumJets += jets.at(*i_index).pt();
    }
    const int numberOfJets(recoObjectIndices.jetIndices_.size());
    const int numberOfTaggedJets(recoObjectIndices.bjetIndices_.size());
    const int numberOfUntaggedJets(numberOfJets - numberOfTaggedJets);
    const double btagDiscriminatorAverage_tagged = numberOfTaggedJets>0 ? btagDiscriminatorSumTagged/static_cast<double>(numberOfTaggedJets) : 0.;
    const double btagDiscriminatorAverage_untagged = numberOfUntaggedJets>0 ? btagDiscriminatorSumUntagged/static_cast<double>(numberOfUntaggedJets) : 0.;
    const double ptSumJetsLeptons = ptSumJets + leptons.at(recoObjectIndices.leptonIndex_).pt() + leptons.at(recoObjectIndices.antiLeptonIndex_).pt();
    
    // Calculate several dijet dependent quantities
    int numberOfHiggsLikeDijet15(0);
    double higgsLikeDijetMass(-999.);
    double higgsLikeDijetMass2(-999.);
    for(const auto& indexPair : recoObjectIndices.jetIndexPairs_){
        const bool hasBtag = (std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), indexPair.first) != recoObjectIndices.bjetIndices_.end()) || 
                             (std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), indexPair.second) != recoObjectIndices.bjetIndices_.end());
        constexpr double higgsMass(125.);
        const double dijetMass = (jets.at(indexPair.first) + jets.at(indexPair.second)).M();
        if(std::abs(higgsMass - dijetMass) < std::abs(higgsMass - higgsLikeDijetMass)) higgsLikeDijetMass = dijetMass;
        if(hasBtag){
            if(std::abs(higgsMass - dijetMass) < std::abs(higgsMass - higgsLikeDijetMass2)) higgsLikeDijetMass2 = dijetMass;
            if(std::abs(dijetMass - higgsMass) < 15.) ++numberOfHiggsLikeDijet15;
        }
    }
    

    // Calculate several new jet-dependent quantities                                                                                                                          
    double minDeltaRTagTag(999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag1;

      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag1      = jets.at(*i_index);
      }
      else
        continue;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag2;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
          tag2      = jets.at(*j_index);
        }
        else
          continue;

        const double deltaR = DeltaR(tag1, tag2);

        if(deltaR < minDeltaRTagTag) {
          minDeltaRTagTag = deltaR;
        }
      }
    }
    
    std::vector<double> deltaRjetjetCollection;

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet1  = jets.at(*i_index);

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV jet2 = jets.at(*j_index);

        const double deltaR = DeltaR(jet1, jet2);

        deltaRjetjetCollection.push_back(deltaR);
      }
    }
    double sumDeltaR_jet_jet = std::accumulate(std::begin(deltaRjetjetCollection), std::end(deltaRjetjetCollection), 0.0);
    double avgDeltaRJetJet =  sumDeltaR_jet_jet/deltaRjetjetCollection.size();

    std::vector<double> deltaRjettagCollection;

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){

      LV jet;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()
           || std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
          jet      = jets.at(*i_index);
          tag      = jets.at(*j_index);
        }
        else
          continue;
        const double deltaR = DeltaR(jet, tag);
        deltaRjettagCollection.push_back(deltaR);
      }
    }
    double sumDeltaR_jet_tag = std::accumulate(std::begin(deltaRjettagCollection), std::end(deltaRjettagCollection), 0.0);
    double avgDeltaRJetTag =  sumDeltaR_jet_tag/deltaRjettagCollection.size();

    std::vector<double> deltaRtagtagCollection;

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag1;

      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag1      = jets.at(*i_index);
      }
      else
        continue;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag2;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
          tag2      = jets.at(*j_index);
        }
        else
          continue;

        const double deltaR = DeltaR(tag1, tag2);

        deltaRtagtagCollection.push_back(deltaR);
      }
    }
    double sumDeltaR_tag_tag = std::accumulate(std::begin(deltaRtagtagCollection), std::end(deltaRtagtagCollection), 0.0);
    double avgDeltaRTagTag =  sumDeltaR_tag_tag/deltaRtagtagCollection.size();


    double maxDeltaEta_jet_jet(-999.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet1 = jets.at(*i_index);

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV jet2  = jets.at(*j_index);

        const double DeltaEta = fabs(jet1.Eta() - jet2.Eta());

        if(DeltaEta > maxDeltaEta_jet_jet) {
          maxDeltaEta_jet_jet = DeltaEta;
        }
      }
    }
    
    double maxDeltaEta_tag_tag(-999.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag1;

      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag1      = jets.at(*i_index);
      }
      else
        continue;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag2;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
          tag2      = jets.at(*j_index);
        }
        else
          continue;

        const double DeltaEta = fabs(tag1.Eta() - tag2.Eta());

        if(DeltaEta > maxDeltaEta_tag_tag) {
          maxDeltaEta_tag_tag = DeltaEta;
        }
      }
    }


    double sumJetPt(0.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet = jets.at(*i_index);
      sumJetPt += jet.pt();
    }


    double sumTagPt(0.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag;

      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag      = jets.at(*i_index);
      }
      else
        continue;
      sumTagPt += tag.pt();
    }


    double minDeltaRjetjet(999.);
    double pT_jet_jet_min_deltaR(-999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet1 = jets.at(*i_index);

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV jet2 = jets.at(*j_index);

        const double deltaR = DeltaR(jet1, jet2);

        if(deltaR < minDeltaRjetjet) {
          minDeltaRjetjet = deltaR;
          //pT_jet_jet_min_deltaR = fabs(jet1.pt() + jet2.pt());                                                                                              
          pT_jet_jet_min_deltaR = (jet1 + jet2).pt();
        }
      }
    }


    double minDeltaRjettag(999.);
    double pT_jet_tag_min_deltaR(-999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()
           || std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
          jet = jets.at(*i_index);
          tag      = jets.at(*j_index);
        }
        else
          continue;

        const double deltaR = DeltaR(jet, tag);

        if(deltaR < minDeltaRjettag) {
          minDeltaRjettag = deltaR;
          //pTjettag_min_deltaR = fabs(jet.pt() + tag.pt());                                                                                                  
          pT_jet_tag_min_deltaR = (jet + tag).pt();
        }
      }
    }


    double minDeltaRtagtag(999.);
    double pT_tag_tag_min_deltaR(-999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag1;

      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag1      = jets.at(*i_index);
      }
      else
        continue;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag2;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
          tag2      = jets.at(*j_index);
        }
        else
          continue;

        const double deltaR = DeltaR(tag1, tag2);

        if(deltaR < minDeltaRtagtag) {
          minDeltaRtagtag = deltaR;
          //pT_tag_tag_min_deltaR = fabs(tag1.pt() + tag2.pt());                                                                                              
          pT_tag_tag_min_deltaR = (tag1 + tag2).pt();
        }
      }
    }


    double mass_jet_jet_jet_max_pT(-999.);
    double max_sum_jet_pT(-999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet1  = jets.at(*i_index);

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV jet2  = jets.at(*j_index);

        for(auto k_index = j_index + 1; k_index != recoObjectIndices.jetIndices_.end(); ++k_index){

          LV jet3  = jets.at(*k_index);

          //double sum_jet_pT = jet1.pt()+jet2.pt()+jet3.pt();                                                                                                
          double sum_jet_pT = (jet1+jet2+jet3).pt();

          if(sum_jet_pT > max_sum_jet_pT) {
            max_sum_jet_pT = sum_jet_pT;
            mass_jet_jet_jet_max_pT = (jet1+jet2+jet3).M();
          }
        }
      }
    }


    double mass_jet_tag_tag_max_pT(-999.);
    double max_sum_tag_pT(-999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet1;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV jet2;

        for(auto k_index = j_index + 1; k_index != recoObjectIndices.jetIndices_.end(); ++k_index){
          LV jet3;

          if((std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()
              && std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end())
             || (std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()
                 && std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *k_index) != recoObjectIndices.bjetIndices_.end())
             || (std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()
                 && std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *k_index) != recoObjectIndices.bjetIndices_.end()))
            {
              jet1     = jets.at(*i_index);
              jet2     = jets.at(*j_index);
              jet3     = jets.at(*k_index);
            }
          else
            continue;

          //double sum_tag_pT = jet1.pt()+jet2.pt()+jet3.pt();                                                                                                
          double sum_tag_pT = (jet1+jet2+jet3).pt();

          if(sum_tag_pT > max_sum_tag_pT) {
            max_sum_tag_pT = sum_tag_pT;
            mass_jet_tag_tag_max_pT = (jet1+jet2+jet3).M();
          }
        }
      }
    }

    
    double mass_jet_jet_min_deltaR(-999.);
    double minDeltaR_JetJet(999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet1 = jets.at(*i_index);

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV jet2 = jets.at(*j_index);

        const double deltaR = DeltaR(jet1, jet2);

        if(deltaR < minDeltaR_JetJet) {
          minDeltaR_JetJet = deltaR;
          mass_jet_jet_min_deltaR = (jet1+jet2).M();
        }
      }
    }

    

    double mass_jet_tag_min_deltaR(-999.);
    double minDeltaR_JetTag(999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet;
      LV tag;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        //  LV tag;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()
           || std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){

          jet      = jets.at(*i_index);
          tag      = jets.at(*j_index);
        }
        else
          continue;
      
        const double deltaR = DeltaR(jet, tag);
        
        if(deltaR < minDeltaR_JetTag) {
          minDeltaR_JetTag = deltaR;
          mass_jet_tag_min_deltaR = (jet+tag).M();
        }
      }
    }
   
    
    double mass_tag_tag_min_deltaR(-999.);
    double minDeltaR_TagTag(999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag1;
      
      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag1      = jets.at(*i_index);
      }
      else
        continue;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag2;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
          tag2      = jets.at(*j_index);
        }
        else
          continue;

        const double deltaR = DeltaR(tag1, tag2);

        if(deltaR < minDeltaR_TagTag) {
          minDeltaR_TagTag = deltaR;
          mass_tag_tag_min_deltaR = (tag1+tag2).M();
        }
      }
    }
    


    double mass_tag_tag_max_mass(-999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag1;

      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag1      = jets.at(*i_index);
      }
      else
        continue;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        LV tag2;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
          tag2      = jets.at(*j_index);
        }
        else
          continue;

        const double mass = (tag1+tag2).M();

        if(mass > mass_tag_tag_max_mass) {
          mass_tag_tag_max_mass = mass;
        }
      }
    }


    double median_mass_jet_jet(-999.);
    std::vector<double> massJetJetCollection;

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){

      LV jet1 = jets.at(*i_index);

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){

        LV jet2 = jets.at(*j_index);

        const double mass = (jet1 + jet2).M();

        massJetJetCollection.push_back(mass);
      }
    }

    std::sort(massJetJetCollection.begin(),massJetJetCollection.end());
    size_t n = massJetJetCollection.size();

    if (n != 0) {

      if((int)n%2 == 0)
        median_mass_jet_jet = (massJetJetCollection.at(n/2) + massJetJetCollection.at(n/2-1))/2;
      else
        median_mass_jet_jet = massJetJetCollection.at((n-1)/2);
    }




    // Centrality calculations
    double centrality_jets_leps(-999.);
    double centrality_tags(-999.0);
    double sumJetPT(.0);
    double sumJetE(.0);
    double sumTagPT(.0);
    double sumTagE(.0);    

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){                                   
      
      LV jet = jets.at(*i_index); 
      sumJetPT += jet.pt();
      sumJetE += jet.E();
      
      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){         
        sumTagPT += jet.pt();
        sumTagE += jet.E();
      }
    }

    LV lepton = leptons.at(recoObjectIndices.leptonIndex_);
    LV antilepton = leptons.at(recoObjectIndices.antiLeptonIndex_);

    centrality_jets_leps = (sumJetPT + lepton.pt() + antilepton.pt())/(sumJetE + lepton.E() + antilepton.E());
    centrality_tags = sumTagPT/sumTagE;



    // Twist between jet pair 
    double max_mass_jet_jet(0.);
    double twist_jet_jet_max_mass(-999.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      
      LV jet1 = jets.at(*i_index);

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){

        LV jet2 = jets.at(*j_index); 

        const double mass = (jet1 + jet2).M();

        if(mass > max_mass_jet_jet) {
          max_mass_jet_jet = mass;
          twist_jet_jet_max_mass = TMath::ATan(ROOT::Math::VectorUtil::DeltaPhi(jet1, jet2)/(jet1.Eta() - jet2.Eta()));
        }
        
      }
    }


    // Twist between jet and b-tagged jet pair
    double max_mass_jet_tag(0.);
    double twist_jet_tag_max_mass(-999.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){

      LV jet1;

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){

        LV jet2;

        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()
           || std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){

          jet1 = jets.at(*i_index);
          jet2 = jets.at(*j_index);

          const double mass = (jet1 + jet2).M();

          if(mass > max_mass_jet_tag) {
            max_mass_jet_tag = mass;
            twist_jet_tag_max_mass = TMath::ATan(ROOT::Math::VectorUtil::DeltaPhi(jet1, jet2)/(jet1.Eta() - jet2.Eta()));
          }

        }
      }
    }


    // Twist between b-tagged jet pair 
    double max_mass_tag_tag(0.);
    double twist_tag_tag_max_mass(-999.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        
        LV tag1;
        LV tag2;
        
        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
              
          if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
                    
            tag1 = jets.at(*i_index);
            tag2 = jets.at(*j_index);
                    
            const double mass = (tag1 + tag2).M();

            if(mass > max_mass_tag_tag) {
              max_mass_tag_tag = mass;
              twist_tag_tag_max_mass = TMath::ATan(ROOT::Math::VectorUtil::DeltaPhi(tag1, tag2)/(tag1.Eta() - tag2.Eta()));
            }
                    
          }  
        }
      }
    }


    // Twist between b-tagged jet pair 
    double min_deltaR_tag_tag(999.);
    double twist_tag_tag_min_deltaR(-999.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){

      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        
        LV tag1;
        LV tag2;
        
        if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
            
          if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end()){
                
            tag1 = jets.at(*i_index);
            tag2 = jets.at(*j_index);
                
            const double deltaR = DeltaR(tag1, tag2);
                
            if(deltaR < min_deltaR_tag_tag) {
              min_deltaR_tag_tag = deltaR;
              twist_tag_tag_min_deltaR = TMath::ATan(ROOT::Math::VectorUtil::DeltaPhi(tag1, tag2)/(tag1.Eta() - tag2.Eta()));
            }

          }  
        }
      }
    }


    // Event shape variable for jets in the event
    std::vector<LV> recoJetCollection;

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV jet = jets.at(*i_index);

      recoJetCollection.push_back(jet);
    }
   
    EventShapeVariables eventshape_jets(recoJetCollection);
    
    // Spherecity eigenvalue varaibles jets 
    double sphericity_jet  = eventshape_jets.sphericity();
    double aplanarity_jet  = eventshape_jets.aplanarity();
    double circularity_jet = eventshape_jets.circularity();
    double isotropy_jet    = eventshape_jets.isotropy();
    double C_jet           = eventshape_jets.C();
    double D_jet           = eventshape_jets.D();
    double transSphericity_jet = eventshape_jets.transSphericity();

    // Fox Wolfram moments variables
    double H0_jet = eventshape_jets.H(0);
    double H1_jet = eventshape_jets.H(1);
    double H2_jet = eventshape_jets.H(2);
    double H3_jet = eventshape_jets.H(3);
    double H4_jet = eventshape_jets.H(4);

    double R1_jet = eventshape_jets.R(1);
    double R2_jet = eventshape_jets.R(2);
    double R3_jet = eventshape_jets.R(3);
    double R4_jet = eventshape_jets.R(4);

 
    // Event shape variables for b-tag jets in the event
    std::vector<LV> recoBJetCollection;

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      LV tag;

      // Select b-tagged jets
      if(std::find(recoObjectIndices.bjetIndices_.begin(),recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
        tag  = jets.at(*i_index);
      }
      else
        continue;

      recoBJetCollection.push_back(tag);
    }

    EventShapeVariables eventshape_tags(recoBJetCollection);        

    // Sphericity associated variables 
    double sphericity_tag  = eventshape_tags.sphericity();
    double aplanarity_tag  = eventshape_tags.aplanarity();
    double circularity_tag = eventshape_tags.circularity();
    double isotropy_tag    = eventshape_tags.isotropy();
    double C_tag           = eventshape_tags.C();
    double D_tag           = eventshape_tags.D();
    double transSphericity_tag = eventshape_tags.transSphericity();

    // Fox Wolfram moments associated variables
    double H0_tag = eventshape_tags.H(0);
    double H1_tag = eventshape_tags.H(1);
    double H2_tag = eventshape_tags.H(2);
    double H3_tag = eventshape_tags.H(3);
    double H4_tag = eventshape_tags.H(4);

    double R1_tag = eventshape_tags.R(1);
    double R2_tag = eventshape_tags.R(2);
    double R3_tag = eventshape_tags.R(3);
    double R4_tag = eventshape_tags.R(4);


    std::pair<int, int> p_mass_jj(-999,-999);
    double mass_jj(-999.);

    int i_index_maxCSV(-999);
    int j_index_maxCSV(-999);
    double maxCSV1(-999.);
    double maxCSV2(-999.);

    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
      //const bool btagged1 = std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end();
      
      for(auto j_index = i_index + 1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
        //const bool btagged2 = std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), *j_index) != recoObjectIndices.bjetIndices_.end();

        // Veto jet if it's been tagged as most likely coming from the top-antitop system
        if((*i_index == topPair.first || *i_index == topPair.second) 
           || (*j_index == topPair.first || *j_index == topPair.second))
          continue;
        
        // Require at least one b-tagged jet else continue (study requiring at least one b-tagged jet)
        //if(!btagged1 && !btagged2) continue;
        
        if (jetBtags.at(*i_index) > maxCSV1)
          i_index_maxCSV = *i_index;
        if (jetBtags.at(*j_index) > maxCSV2)
          j_index_maxCSV = *j_index;
        
        p_mass_jj =  std::make_pair(i_index_maxCSV, j_index_maxCSV);
        
      }
    }

    mass_jj = (p_mass_jj != std::make_pair(-999,-999)) ? (jets.at(p_mass_jj.first) + jets.at(p_mass_jj.second)).M() : -999.;


    return new MvaVariablesEventClassification(eventMetadata,
                                               numberOfJets,
                                               btagDiscriminatorAverage_tagged, btagDiscriminatorAverage_untagged,
                                               minDeltaRJetJet, minDeltaRTagTag,
                                               avgDeltaRJetJet, avgDeltaRJetTag, avgDeltaRTagTag,
                                               ptSumJetsLeptons,
                                               numberOfHiggsLikeDijet15,
                                               higgsLikeDijetMass, higgsLikeDijetMass2,
                                               mass_jet_jet_min_deltaR, mass_tag_tag_min_deltaR, mass_jet_tag_min_deltaR,
                                               mass_tag_tag_max_mass,
                                               median_mass_jet_jet,
                                               maxDeltaEta_jet_jet, maxDeltaEta_tag_tag,
                                               sumJetPt, sumTagPt,
                                               pT_jet_jet_min_deltaR, pT_jet_tag_min_deltaR, pT_tag_tag_min_deltaR,
                                               mass_jet_jet_jet_max_pT, mass_jet_tag_tag_max_pT,
                                               centrality_jets_leps, centrality_tags, 
                                               twist_jet_jet_max_mass, twist_jet_tag_max_mass, twist_tag_tag_max_mass,
                                               twist_tag_tag_min_deltaR,
                                               sphericity_jet,  aplanarity_jet,  circularity_jet,  isotropy_jet,
                                               C_jet,  D_jet, transSphericity_jet,
                                               sphericity_tag,  aplanarity_tag,  circularity_tag,  isotropy_tag,
                                               C_tag,  D_tag,  transSphericity_tag,
                                               H0_jet,  H1_jet,  H2_jet,  H3_jet,  H4_jet,
                                               R1_jet,  R2_jet,  R3_jet,  R4_jet,
                                               H0_tag,  H1_tag,  H2_tag,  H3_tag,  H4_tag,
                                               R1_tag,  R2_tag,  R3_tag,  R4_tag,
                                               mass_jj,
                                               eventWeight);
}


DLBDTMvaVariablesEventClassification* temp::dlBdtFillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                               const RecoObjects& recoObjects)
{
    // Access relevant objects and indices
    const VLV& leptons(*recoObjects.allLeptons_);
    const std::vector<int>& leptonPdgIds(*recoObjects.lepPdgId_);
    const VLV& jets(*recoObjects.jets_);
    const std::vector<double>& jetBtags(*recoObjects.jetBtags_);
    
    // Objects in input format for CommonClassifier
    // FIXME: Make btagWP there a single value per class
    std::vector<TLorentzVector> v_lepton;
    std::vector<double> v_leptonCharge;
    std::vector<TLorentzVector> v_jet;
    std::vector<double> v_jetBtag;
    constexpr double btagWP(0.800);
    
    for(const int index : recoObjectIndices.allLeptonIndices_){
        v_lepton.push_back(common::LVtoTLV(leptons.at(index)));
        v_leptonCharge.push_back(leptonPdgIds.at(index)>0 ? -1. : 1.);
    }
    
    for(const int index : recoObjectIndices.jetIndices_){
        v_jet.push_back(common::LVtoTLV(jets.at(index)));
        v_jetBtag.push_back(jetBtags.at(index));
    }
    
    
    return DLBDTMvaVariablesEventClassification::fillVariables(v_lepton, v_leptonCharge, v_jet, v_jetBtag, btagWP);
}



MvaVariablesEventClassification* temp::fillVariables(const EventMetadata& eventMetadata,
                                                     const tth::RecoObjectIndices& recoObjectIndices,
                                                     const RecoObjects& recoObjects,
                                                     const double& eventWeight)
{
    DLBDTMvaVariablesEventClassification* dlBdtMvaVariables = dlBdtFillVariables(recoObjectIndices, recoObjects);
    DLBDTMvaVariablesEventClassification dlBdtMvaVariablesTemp = *dlBdtMvaVariables;
    delete dlBdtMvaVariables;
    return new MvaVariablesEventClassification(eventMetadata, eventWeight, dlBdtMvaVariablesTemp);
}






// ---------------------------------- Class MvaVariablesEventClassification::EventShapeVariables -------------------------------------------



/// constructor from XYZ coordinates
MvaVariablesEventClassification::EventShapeVariables::EventShapeVariables(const VLV& inputVectors):
inputVectors_(makeVecForEventShape(inputVectors))
{}



std::vector<ROOT::Math::XYZVector> MvaVariablesEventClassification::EventShapeVariables::makeVecForEventShape(const VLV& jets)
{
  std::vector<ROOT::Math::XYZVector> p;

  for(unsigned int i=0; i< jets.size(); i++){

    ROOT::Math::XYZVector Vjet;

    Vjet = ROOT::Math::XYZVector(jets[i].Px(), jets[i].Py(), jets[i].Pz());
    p.push_back(Vjet);
  }

  return p;
}



/// needs the number of steps to determine how fine the granularity of the algorithm in phi should be
double MvaVariablesEventClassification::EventShapeVariables::isotropy(const unsigned int& numberOfSteps) const
{
  const double deltaPhi=2*TMath::Pi()/numberOfSteps;
  double phi = 0, eIn =-1., eOut=-1.;
 
  for(unsigned int i=0; i<numberOfSteps; ++i){
    phi+=deltaPhi;
    double sum=0;
    
    for(unsigned int j=0; j<inputVectors_.size(); ++j){

      // sum over inner product of unit vectors and momenta
      sum+=TMath::Abs(TMath::Cos(phi)*inputVectors_[j].x()+TMath::Sin(phi)*inputVectors_[j].y());
    }
    
    if( eOut<0. || sum<eOut ) eOut=sum;
    if( eIn <0. || sum>eIn  ) eIn =sum;
  }
  return (eIn-eOut)/eIn;
}



/// the return value is 1 for spherical and 0 linear events in r-phi.
double MvaVariablesEventClassification::EventShapeVariables::circularity(const unsigned int& numberOfSteps) const
{
  const double deltaPhi=2*TMath::Pi()/numberOfSteps;
  double circularity=-1, phi=0, area = 0;
 
  for(unsigned int i=0;i<inputVectors_.size();i++) {
    area+=TMath::Sqrt(inputVectors_[i].x()*inputVectors_[i].x()+inputVectors_[i].y()*inputVectors_[i].y());
  }
  
  for(unsigned int i=0; i<numberOfSteps; ++i){
    phi+=deltaPhi;
    double sum=0, tmp=0.;
    
    for(unsigned int j=0; j<inputVectors_.size(); ++j){
      sum+=TMath::Abs(TMath::Cos(phi)*inputVectors_[j].x()+TMath::Sin(phi)*inputVectors_[j].y());
    }
    
    tmp=TMath::Pi()/2*sum/area;
    
    if(circularity<0 || tmp<circularity){
      circularity=tmp;
    }
  }
  return circularity;
}



// helper function to fill the 3 dimensional momentum tensor from the inputVecotrs where needed
TMatrixDSym MvaVariablesEventClassification::EventShapeVariables::compMomentumTensor(double r) const
{
  TMatrixDSym momentumTensor(3);
  momentumTensor.Zero();

  if (inputVectors_.size() < 2){ 
    return momentumTensor;
  }
  
   // fill momentumTensor from inputVectors
   double norm = 1.;
   
   for ( int i = 0; i < (int)inputVectors_.size(); ++i ){

     double p2 = inputVectors_[i].Dot(inputVectors_[i]);

     double pR = (r == 2.) ? p2 : TMath::Power(p2, 0.5*r);
     norm += pR;

     double pRminus2 = (r == 2.) ? 1. : TMath::Power(p2, 0.5*r - 1.);

     momentumTensor(0,0) += pRminus2*inputVectors_[i].x()*inputVectors_[i].x();
     momentumTensor(0,1) += pRminus2*inputVectors_[i].x()*inputVectors_[i].y();
     momentumTensor(0,2) += pRminus2*inputVectors_[i].x()*inputVectors_[i].z();
     momentumTensor(1,0) += pRminus2*inputVectors_[i].y()*inputVectors_[i].x();
     momentumTensor(1,1) += pRminus2*inputVectors_[i].y()*inputVectors_[i].y();
     momentumTensor(1,2) += pRminus2*inputVectors_[i].y()*inputVectors_[i].z();
     momentumTensor(2,0) += pRminus2*inputVectors_[i].z()*inputVectors_[i].x();
     momentumTensor(2,1) += pRminus2*inputVectors_[i].z()*inputVectors_[i].y();
     momentumTensor(2,2) += pRminus2*inputVectors_[i].z()*inputVectors_[i].z();
   }

   return (1./norm)*momentumTensor;
}



/// helper function to fill the 3 dimensional vector of eigen-values;
/// the largest (smallest) eigen-value is stored at index position 0 (2)
TVectorD MvaVariablesEventClassification::EventShapeVariables::compEigenValues(double r) const
{
  TVectorD eigenValues(3);
  TMatrixDSym myTensor = compMomentumTensor(r);
  
  if(myTensor.IsSymmetric()){
    if(myTensor.NonZeros() != 0) myTensor.EigenVectors(eigenValues);
   }

   return eigenValues;
}



/// Return values are 1 for spherical, 3/4 for plane and 0 for linear events
double MvaVariablesEventClassification::EventShapeVariables::sphericity(double r) const
{
  TVectorD eigenValues = compEigenValues(r);

  return 1.5*(eigenValues(1) + eigenValues(2));
}



/// Return values are 0.5 for spherical and 0 for plane and linear events
double MvaVariablesEventClassification::EventShapeVariables::aplanarity(double r) const
{
  TVectorD eigenValues = compEigenValues(r);
  
  return 1.5*eigenValues(2);
}



/// Return value is between 0 and 1 and measures the 3-jet structure of the event (C vanishes for a "perfect" 2-jet event)
double MvaVariablesEventClassification::EventShapeVariables::C(double r) const
{
  TVectorD eigenValues = compEigenValues(r);
  
   return 3.*(eigenValues(0)*eigenValues(1) + eigenValues(0)*eigenValues(2) + eigenValues(1)*eigenValues(2));
}



/// Return value is between 0 and 1 and measures the 4-jet structure of the event (D vanishes for a planar event)
double MvaVariablesEventClassification::EventShapeVariables::D(double r) const
{
  TVectorD eigenValues = compEigenValues(r);
  
  return 27.*eigenValues(0)*eigenValues(1)*eigenValues(2);
}



/// Return value is between 0 and 1, value of 0 is "pencile-like" limit, while 1 for "isotropic-like" limit  
double MvaVariablesEventClassification::EventShapeVariables::transSphericity(double r) const 
{
  
  TVectorD eigenValues = compEigenValues(r);
  
  return 2.*eigenValues(1)/(eigenValues(0) + eigenValues(1));
}



/// Calcualtes the Fox-Wolfram moment for a given order
double MvaVariablesEventClassification::EventShapeVariables::H(int i) const 
{

  FoxWolframMoments fwm(i);
  fwm.Compute(inputVectors_);
  
  return fwm.H(i);
}



// Calcualtes the 0-th Fox-Wolfram moment
double MvaVariablesEventClassification::EventShapeVariables::ZerothMoment() const 
{
  FoxWolframMoments fwm(0);
  fwm.Compute(inputVectors_);
  
  return fwm.ZerothMoment();
}



// Calculates the ratio between Fox-Wolfram moments
double MvaVariablesEventClassification::EventShapeVariables::R(int order) const
{
  
  FoxWolframMoments fwm(order);
  fwm.Compute(inputVectors_);

  return fwm.R(order);
}








// ---------------------------------- Class MvaVariablesEventClassification::FoxWolframMoments -------------------------------------------



MvaVariablesEventClassification::FoxWolframMoments::FoxWolframMoments(int maxorder):
_nmom( maxorder+1 ),
_FWarray( _nmom ),
_sumarray( _nmom )
{}



void MvaVariablesEventClassification::FoxWolframMoments::Compute(const std::vector<ROOT::Math::XYZVector>& inputVectors)
{
  // initialize
  Reset();
  
  if(inputVectors.size()== 0) 
    return;
  
  double s = 0.;
  int l;
  
  // Start a loop over the all particle candidates
  for(unsigned int i=0; i<inputVectors.size(); ++i){

    // Candidate particle's 3-momentum
    TVector3 p1(inputVectors[i].X(),inputVectors[i].Y(),inputVectors[i].Z());
    double pmag1 = p1.Mag();
    
    // loop over other particle's candidates, starting at the next one in the list
    for(unsigned int j=i; j<inputVectors.size(); ++j){
      
      // Candidate particle's 3-momentum
      TVector3 p2(inputVectors[j].X(),inputVectors[j].Y(),inputVectors[j].Z());
      double pmag2 = p2.Mag();
        
      // the cosine of the angle between the two candidate particles
      double cosPhi =  cos( p1.Angle(p2) );
      
      // the contribution of this pair of track
      // (note the factor 2 : the pair enters the sum twice)
      for( l=0; l<_nmom; l++ )
    _sumarray(l) += 2 * pmag1 * pmag2 * Legendre( l, 0, cosPhi );
    }
    
    // contribution for this moment
    for( l=0; l<_nmom; l++ )
      _sumarray(l) += pmag1*pmag1*Legendre( l, 0, 1. );
      
    // total energy
    s += p1.Mag();
    
  }
  
  // In case of
  if( s<=0. ) return;
   
  // normalize Fox Wolfram Moments
  for(int i=0; i<_nmom; i++)
    _FWarray(i) = _sumarray(i)/pow(s,2) ;
}



void MvaVariablesEventClassification::FoxWolframMoments::Reset() 
{
  for ( int i = 0; i< _nmom; ++i) 
    { 
      _FWarray(i) = 0.;
      _sumarray(i) = 0.;
    }
}



double MvaVariablesEventClassification::FoxWolframMoments::R( int order ) const
{
  if(H(0) > 0.) {
    if(order < _nmom){ 
      return (H(order)/H(0));
    }
  }
  return 0.;
}



double MvaVariablesEventClassification::FoxWolframMoments::Legendre( int l, int m, double x )
{
  assert(m >= 0.);
  assert(m <= l);
  assert(fabs(x) <= 1.);
  
  double pmm = 1.;
  
  if(m > 0)
    {
      double somx2 = sqrt((1. - x) * (1. + x));
      double fact = 1.;
      
      for(int i=0; i<m; i++)
    {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
    }
  
  if(l == m)
    return pmm;
  
  else
    {
      double pmmp1 = x * (2 * m + 1) * pmm;
   
      if(l == m + 1)
    return pmmp1;
      else
    {
      for(int ll = m+2; ll <= l; ++ll)
        {
          double pll = (x*(2 * ll - 1)*pmmp1 - (ll + m - 1)*pmm)/(ll - m);
          pmm = pmmp1;
          pmmp1 = pll;
        }
      
      return pmmp1;
    }
    }
}






