#ifndef flightDirection_h_1
#define flightDirection_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "TVector3.h"
#include "TF2.h"
#include <marlinutil/HelixClass.h>
#include "SLDCorrection.h"

int getPrimaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , bool cheatVertices , std::vector<double> &primaryVertex , bool m_displayEvent );
int getParentHadronFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , TVector3 &recoFlightDirection , std::string inputPrimaryVertex , std::string inputBuildUpVertex , std::string inputJetCollection , int vertexingScenario , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent , SLDCorrection* thisProcessor );
int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , std::string inputBuildUpVertex , bool cheatVertices , std::vector<double> &secondayVertex , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent );
int getJetAxis( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &jetAxis , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
int getLeadingParticleFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &leadingParticleFlightDirection , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
float intersectHelixHelix( EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 Momentum , std::vector<double> point , std::vector<double> &PCAatLeptonTrack , std::vector<double> &PCAatDownStreamLine );
float intersectHelixLine( EVENT::LCEvent *pLCEvent , EVENT::Track *track , TVector3 Momentum , std::vector<double> point , std::vector<double> &PCAatLeptonTrack , std::vector<double> &PCAatDownStreamLine , std::string inputPrimaryVertex );
int getTrueDownStreamVertex( EVENT::MCParticle *SLDLepton , TVector3 &trueDSVertex , TVector3 &trueDSVertexMomentum );
int getTrueVertex( EVENT::MCParticle *MotherParticle , TVector3 &trueDSVertex ,TVector3 &trueDSVertexMomentum );


#endif
