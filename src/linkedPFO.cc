#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include "linkedPFO.h"

using namespace lcio;
using namespace marlin;

/*
ReconstructedParticle* getLinkedPFO( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *visibleMCP , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , bool getChargedTLV , bool getNeutralTLV )
{
	streamlog_out(DEBUG0) << "	Look for PFO linked to MCParticle (PDG: " << visibleMCP->getPDG() << " , GenStat: " << visibleMCP->getGeneratorStatus() << ")" << std::endl;
	ReconstructedParticle* linkedPFO{};
	bool foundlinkedPFO = false;
	try
	{
		LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( recoMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( mcTruthRecoLinkCollection ) );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for PFO <-> MCParticle" << std::endl;
		return NULL;
	}
	LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( mcTruthRecoLinkCollection ) );
	const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( visibleMCP );
	const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( visibleMCP );
	streamlog_out(DEBUG0) << "	Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
	double maxweightPFOtoMCP = 0.;
	double maxweightMCPtoPFO = 0.;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
	{
		double pfo_weight = 0.0;
		if ( getChargedTLV )
		{
			pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
		}
		else if ( getNeutralTLV )
		{
			pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
		}
		streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << int( PFOweightvec.at( i_pfo ) ) % 10000 / 1000.0 << " , Cluster: " << int( PFOweightvec.at( i_pfo ) ) / 10000 / 1000.0 << ")" << std::endl;
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
		if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightMCPtoPFO = pfo_weight;
			iMCPtoPFOmax = i_pfo;
			streamlog_out(DEBUG0) << "	PFO at index: " << i_pfo << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
		}
	}
	if ( iMCPtoPFOmax != -1 )
	{
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( iMCPtoPFOmax );
		LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( recoMCTruthLinkCollection ) );
		const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( testPFO );
		const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( testPFO );
		for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
		{
			double mcp_weight = 0.0;
			if ( getChargedTLV )
			{
				mcp_weight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
			}
			else if ( getNeutralTLV )
			{
				mcp_weight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
			}
			MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
			if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightPFOtoMCP = mcp_weight;
				iPFOtoMCPmax = i_mcp;
				streamlog_out(DEBUG0) << "	MCParticle at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
			}
		}
		if ( iPFOtoMCPmax != -1 )
		{
			if ( MCPvec.at( iPFOtoMCPmax ) == visibleMCP )
			{
				streamlog_out(DEBUG0) << "	Linked PFO to MCParticle found successfully " << std::endl;
				linkedPFO = testPFO;
				foundlinkedPFO = true;
			}
		}
	}

	if( foundlinkedPFO )
	{
		streamlog_out(DEBUG0) << "	Found linked RecoLepton (px,py,pz,E): ( " << linkedPFO->getMomentum()[ 0 ] << " 	, " << linkedPFO->getMomentum()[ 1 ] << " 	, " << linkedPFO->getMomentum()[ 2 ] << " 	, " << linkedPFO->getEnergy() << " 	)" << std::endl;
		return linkedPFO;
	}
	else
	{
		streamlog_out(DEBUG0) << "	Couldn't Find a PFO linked to MCParticle" << std::endl;
		return NULL;
	}
}

MCParticle* getLinkedMCP( EVENT::LCEvent *pLCEvent , EVENT::ReconstructedParticle *recoParticle , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , bool getChargedMCP , bool getNeutralMCP )
{
	streamlog_out(DEBUG0) << "	Look for MCP linked to ReconstructedParticle (TYPE: " << recoParticle->getType() << " , nTracks: " << recoParticle->getTracks().size() << " , nClusters: " << recoParticle->getClusters().size() << ")" << std::endl;
	MCParticle* linkedMCP{};
	bool foundlinkedMCP = false;
	try
	{
		LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( recoMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( mcTruthRecoLinkCollection ) );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for PFO <-> MCParticle" << std::endl;
		return NULL;
	}
	double maxweightPFOtoMCP = 0.;
	double maxweightMCPtoPFO = 0.;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( recoMCTruthLinkCollection ) );
	const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( recoParticle );
	const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( recoParticle );
	for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
	{
		double mcp_weight = 0.0;
		if ( getChargedMCP )
		{
			mcp_weight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
		}
		else if ( getNeutralMCP )
		{
			mcp_weight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
		}
		MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
		if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightPFOtoMCP = mcp_weight;
			iPFOtoMCPmax = i_mcp;
			streamlog_out(DEBUG0) << "	MCParticle at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
		}
	}
	if ( iPFOtoMCPmax != -1 )
	{
		MCParticle* testMCP = ( MCParticle* )MCPvec.at( iPFOtoMCPmax );
		LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( mcTruthRecoLinkCollection ) );
		const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( testMCP );
		const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( testMCP );
		streamlog_out(DEBUG0) << "	Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
		for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
		{
			double pfo_weight = 0.0;
			if ( getChargedMCP )
			{
				pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
			}
			else if ( getNeutralMCP )
			{
				pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
			}
			streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << int( PFOweightvec.at( i_pfo ) ) % 10000 / 1000.0 << " , Cluster: " << int( PFOweightvec.at( i_pfo ) ) / 10000 / 1000.0 << ")" << std::endl;
			ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
			if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoPFO = pfo_weight;
				iMCPtoPFOmax = i_pfo;
				streamlog_out(DEBUG0) << "	PFO at index: " << i_pfo << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
			}
		}
		if ( PFOvec.at( iMCPtoPFOmax ) == recoParticle )
		{
			streamlog_out(DEBUG0) << "	Linked MCParticle to PFO found successfully " << std::endl;
			linkedMCP = testMCP;
			foundlinkedMCP = true;
		}
	}
	if( foundlinkedMCP )
	{
		streamlog_out(DEBUG0) << "	Found linked RecoLepton (px,py,pz,E): ( " << linkedMCP->getMomentum()[ 0 ] << " 	, " << linkedMCP->getMomentum()[ 1 ] << " 	, " << linkedMCP->getMomentum()[ 2 ] << " 	, " << linkedMCP->getEnergy() << " 	)" << std::endl;
		return linkedMCP;
	}
	else
	{
		streamlog_out(DEBUG0) << "	Couldn't Find a PFO linked to MCParticle" << std::endl;
		return NULL;
	}
}

std::vector<EVENT::MCParticle*> getMCParticlesWithVertex( EVENT::MCParticle *MCP , std::vector<EVENT::MCParticle*> MCParticlesWithVertex )
{
	int nDaughters = MCP->getDaughters().size();
	for ( int i_daughter1 = 0 ; i_daughter1 < nDaughters ; ++i_daughter1 )
	{
		MCParticle* daughter1 = MCP->getDaughters()[ i_daughter1 ];
		if ( daughter1->getGeneratorStatus() == 1 )
		{
			if (  daughter1->getCharge() > 0.5 )
			{
				for ( int i_daughter2 = 0 ; i_daughter2 < nDaughters ; ++i_daughter2 )
				{
					MCParticle* daughter2 = MCP->getDaughters()[ i_daughter2 ];
					if ( daughter2->getGeneratorStatus() == 1 && daughter2->getCharge() < -0.5 )
					{
						MCParticlesWithVertex.push_back( MCP );
					}
				}
			}
		}
		else
		{
			MCParticlesWithVertex = getMCParticlesWithVertex( daughter1 , MCParticlesWithVertex );
		}
	}
	return MCParticlesWithVertex;
}

std::vector<EVENT::MCParticle*> getTrueChargedDecayProducts( EVENT::MCParticle *MCP , std::vector<EVENT::MCParticle*> trueChargedDecayProducts )
{
	int nDaughters = MCP->getDaughters().size();
	for ( int i_daughter = 0 ; i_daughter < nDaughters ; ++i_daughter )
	{
		MCParticle* daughter = MCP->getDaughters()[ i_daughter ];
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if (  std::fabs( daughter->getCharge() ) > 0.5 ) trueChargedDecayProducts.push_back( daughter );
		}
		else
		{
			trueChargedDecayProducts = getTrueChargedDecayProducts( daughter , trueChargedDecayProducts );
		}
	}
	return trueChargedDecayProducts;
}

std::vector<EVENT::MCParticle*> getTrueNeutralDecayProducts( EVENT::MCParticle *MCP , std::vector<EVENT::MCParticle*> trueNeutralDecayProducts )
{
	int nDaughters = MCP->getDaughters().size();
	for ( int i_daughter = 0 ; i_daughter < nDaughters ; ++i_daughter )
	{
		MCParticle* daughter = MCP->getDaughters()[ i_daughter ];
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && std::fabs( daughter->getCharge() ) < 0.5 ) trueNeutralDecayProducts.push_back( daughter );
		}
		else
		{
			trueNeutralDecayProducts = getTrueNeutralDecayProducts( daughter , trueNeutralDecayProducts );
		}
	}
	return trueNeutralDecayProducts;
}

float getWidestCosAlphaOfDecayProducts( std::vector<EVENT::MCParticle*> mcpVector , TVector3 parentHadronFlightDirection )
{
	float widestCosAlpha = 1.0;
	for ( unsigned int i_mcp = 0 ; i_mcp < mcpVector.size() ; ++i_mcp )
	{
		TVector3 mcpFlightDirection( mcpVector[ i_mcp ]->getMomentum() );
		mcpFlightDirection.SetMag( 1.0 );
		if ( mcpFlightDirection.Dot( parentHadronFlightDirection ) < widestCosAlpha ) widestCosAlpha = mcpFlightDirection.Dot( parentHadronFlightDirection );
	}
	return widestCosAlpha;
}

float getWidestCosAlphaOfVertices( std::vector<EVENT::MCParticle*> mcpVector , TVector3 parentHadronFlightDirection , std::vector<float> truePrimaryVertex )
{
	float widestCosAlpha = 1.0;
	for ( unsigned int i_mcp = 0 ; i_mcp < mcpVector.size() ; ++i_mcp )
	{
		TVector3 mcpFlightDirection( mcpVector[ i_mcp ]->getEndpoint()[ 0 ] - truePrimaryVertex[ 0 ] , mcpVector[ i_mcp ]->getEndpoint()[ 1 ] - truePrimaryVertex[ 1 ] , mcpVector[ i_mcp ]->getEndpoint()[ 2 ] - truePrimaryVertex[ 2 ] );
		mcpFlightDirection.SetMag( 1.0 );
		if ( mcpFlightDirection.Dot( parentHadronFlightDirection ) < widestCosAlpha ) widestCosAlpha = mcpFlightDirection.Dot( parentHadronFlightDirection );
	}
	return widestCosAlpha;
}
*/
