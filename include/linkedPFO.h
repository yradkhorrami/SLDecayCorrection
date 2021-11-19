#ifndef linkedPFO_h_1
#define linkedPFO_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "TVector3.h"
/*
EVENT::ReconstructedParticle* getLinkedPFO( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *visibleMCP , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , bool getChargedTLV , bool getNeutralTLV );

EVENT::MCParticle* getLinkedMCP( EVENT::LCEvent *pLCEvent , EVENT::ReconstructedParticle *recoParticle , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , bool getChargedMCP , bool getNeutralMCP );
std::vector<EVENT::MCParticle*> getMCParticlesWithVertex( EVENT::MCParticle *MCP , std::vector<EVENT::MCParticle*> MCParticlesWithVertex );

std::vector<EVENT::MCParticle*> getTrueChargedDecayProducts( EVENT::MCParticle *MCP , std::vector<EVENT::MCParticle*> trueChargedDecayProducts );

std::vector<EVENT::MCParticle*> getTrueNeutralDecayProducts( EVENT::MCParticle *MCP , std::vector<EVENT::MCParticle*> trueNeutralDecayProducts );

float getWidestCosAlphaOfDecayProducts( std::vector<EVENT::MCParticle*> mcpVector , TVector3 parentHadronFlightDirection );

float getWidestCosAlphaOfVertices( std::vector<EVENT::MCParticle*> mcpVector , TVector3 parentHadronFlightDirection , std::vector<float> truePrimaryVertex );
*/
#endif
