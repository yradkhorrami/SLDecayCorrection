#ifndef flightDirection_h_1
#define flightDirection_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "TVector3.h"
#include <marlinutil/HelixClass.h>


int getParentHadronFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , TVector3 &recoFlightDirection , std::string inputPrimaryVertex , std::string inputBuildUpVertex , std::string inputJetCollection , int vertexingScenario , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , float &helicesDistance );
int getPrimaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , bool cheatVertices , std::vector<double> &primaryVertex );
int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , std::string inputBuildUpVertex , bool cheatVertices , std::vector<double> &secondayVertex , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , float &helicesDistance );
int getJetAxis( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &jetAxis , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
int getLeadingParticleFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &leadingParticleFlightDirection , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
int intersectLeptonDSVertex( EVENT::ReconstructedParticle *linkedRecoLepton , EVENT::Vertex *downStreamVertex , std::vector<double> &secondayVertex , float &helicesDistance );

#endif
