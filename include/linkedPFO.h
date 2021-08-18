#ifndef linkedPFO_h_1
#define linkedPFO_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"

EVENT::ReconstructedParticle* getLinkedPFO( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *visibleMCP , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , bool getChargedTLV , bool getNeutralTLV );

#endif
