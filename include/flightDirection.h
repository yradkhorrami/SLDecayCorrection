#ifndef flightDirection_h_1
#define flightDirection_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"


int getParentHadronFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , std::string inputBuildUpVertex , std::string inputJetCollection , int vertexinScenario , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
int getPrimaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , bool cheatVertices , std::vector<double> &primaryVertex );
int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , std::string inputBuildUpVertex , bool cheatVertices , std::vector<double> &secondayVertex , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
#endif
