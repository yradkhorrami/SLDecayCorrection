#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include "linkedPFO.h"

using namespace lcio;
using namespace marlin;

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
