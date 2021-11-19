#include "FindParticle.h"

int getNeutralMCPs( EVENT::MCParticle *parentMCP , mcpVector &neutralMCPs )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentMCP->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentMCP->getDaughters()[ i_daughter ];
		streamlog_out(DEBUG0) << "	Checking Daughter[ " << daughter->id() <<" ]: genStatus = " << daughter->getGeneratorStatus() << "	, PDG = " << daughter->getPDG() << "	, Charge = " << daughter->getCharge() << std::endl;
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && fabs( daughter->getCharge() ) <= 0.1 )
			{
				streamlog_out(DEBUG1) << "" << std::endl;
				streamlog_out(DEBUG1) << "Found One Stable Neutral Particle in Semi-Leptonic Decay Products:" << std::endl;
				streamlog_out(DEBUG1) << *daughter << std::endl;
				neutralMCPs.push_back( daughter );
			}
		}
		else
		{
			getNeutralMCPs( daughter , neutralMCPs );
		}
	}
	return neutralMCPs.size();
}

int getChargedMCPs( EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentMCP , mcpVector &chargedMCPs )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentMCP->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentMCP->getDaughters()[ i_daughter ];
		streamlog_out(DEBUG0) << "	Checking Daughter[ " << daughter->id() <<" ]: genStatus = " << daughter->getGeneratorStatus() << "	, PDG = " << daughter->getPDG() << "	, Charge = " << daughter->getCharge() << std::endl;
		if ( daughter->getGeneratorStatus() == 1 && daughter != SLDLepton )
		{
			if ( fabs( daughter->getCharge() ) > 0.1 )
			{
				streamlog_out(DEBUG1) << "" << std::endl;
				streamlog_out(DEBUG1) << "Found One Stable Charged Particle in Semi-Leptonic Decay Products:" << std::endl;
				streamlog_out(DEBUG1) << *daughter << std::endl;
				chargedMCPs.push_back( daughter );
			}
		}
		else
		{
			getChargedMCPs( SLDLepton , daughter , chargedMCPs );
		}
	}
	return chargedMCPs.size();
}

int getTrueVertices( EVENT::MCParticle *parentMCP , mcpVector &MCPsWithVertex )
{
	int nDaughters = parentMCP->getDaughters().size();
	streamlog_out(DEBUG0) << "	Checking Particle [ " << parentMCP->id() <<" ]: genStatus = " << parentMCP->getGeneratorStatus() << "	, PDG = " << parentMCP->getPDG() << "	, Charge = " << parentMCP->getCharge() << " ,	nDaughters = " << nDaughters << std::endl;
	for ( long unsigned int i_daughter1 = 0 ; i_daughter1 < ( parentMCP->getDaughters() ).size() ; ++i_daughter1 )
	{
		EVENT::MCParticle* daughter1 = parentMCP->getDaughters()[ i_daughter1 ];
		if ( daughter1->getGeneratorStatus() == 1 )
		{
			if (  daughter1->getCharge() > 0.5 )
			{
				for ( int i_daughter2 = 0 ; i_daughter2 < nDaughters ; ++i_daughter2 )
				{
					MCParticle* daughter2 = parentMCP->getDaughters()[ i_daughter2 ];
					if ( daughter2->getGeneratorStatus() == 1 && daughter2->getCharge() < -0.5 )
					{
						streamlog_out(DEBUG1) << "" << std::endl;
						streamlog_out(DEBUG1) << "		Found One True Vertex (MCP) in Semi-Leptonic Decay Products:" << std::endl;
						streamlog_out(DEBUG1) << *parentMCP << std::endl;
						bool MCPwasCounted = false;
						for ( unsigned int i = 0 ; i < MCPsWithVertex.size() ; ++i )
						{
							MCParticle* countedMCP = MCPsWithVertex[ i ];
							if ( countedMCP == parentMCP ) MCPwasCounted = true;
						}
						if ( !MCPwasCounted ) MCPsWithVertex.push_back( parentMCP );
					}
				}
			}
		}
		else
		{
			getTrueVertices( daughter1 , MCPsWithVertex );
		}
	}
	return MCPsWithVertex.size();
}

EVENT::MCParticle* getLinkedMCP( EVENT::ReconstructedParticle *recoParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedMCP , bool getNeutralMCP , float &weightPFOtoMCP , float &weightMCPtoPFO )
{
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "	Look for MCP linked to Reconstructed Particle:" << std::endl;
	streamlog_out(DEBUG1) << *recoParticle << std::endl;
	MCParticle* linkedMCP{};
	bool foundlinkedMCP = false;
	double maxweightPFOtoMCP = 0.;
	double maxweightMCPtoPFO = 0.;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( recoParticle );
	const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( recoParticle );
	for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
	{
		double mcp_weight = 0.0;
		double trackWeight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
		double clusterWeight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
		if ( getChargedMCP && !getNeutralMCP )
		{
			mcp_weight = trackWeight;
		}
		else if ( getNeutralMCP && !getChargedMCP )
		{
			mcp_weight = clusterWeight;
		}
		else
		{
			mcp_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
		}
		MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
		if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightPFOtoMCP = mcp_weight;
			iPFOtoMCPmax = i_mcp;
			streamlog_out(DEBUG0) << "	MCParticle at index: " << testMCP->id() << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
		}
	}
	if ( iPFOtoMCPmax != -1 )
	{
		MCParticle *testMCP = (MCParticle *) MCPvec.at( iPFOtoMCPmax );
		const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( testMCP );
		const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( testMCP );
		streamlog_out(DEBUG0) << "	Test Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
		for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
		{
			double pfo_weight = 0.0;
			double trackWeight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
			double clusterWeight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
			if ( getChargedMCP && !getNeutralMCP )
			{
				pfo_weight = trackWeight;
			}
			else if ( getNeutralMCP && !getChargedMCP )
			{
				pfo_weight = clusterWeight;
			}
			else
			{
				pfo_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
			}
			streamlog_out(DEBUG0) << "	Test Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << trackWeight << " , Cluster: " << clusterWeight << ")" << std::endl;
			ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
			if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoPFO = pfo_weight;
				iMCPtoPFOmax = i_pfo;
				streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
			}
		}
		if ( iMCPtoPFOmax != -1 )
		{
			if ( PFOvec.at( iMCPtoPFOmax ) == recoParticle )
			{
				linkedMCP = testMCP;
				foundlinkedMCP = true;
			}
		}
	}
	if( foundlinkedMCP )
	{
		streamlog_out(DEBUG1) << "	Linked MCParticle to PFO found successfully " << std::endl;
		streamlog_out(DEBUG1) << *linkedMCP << std::endl;
		weightPFOtoMCP = maxweightPFOtoMCP;
		weightMCPtoPFO = maxweightMCPtoPFO;
		return linkedMCP;
	}
	else
	{
		streamlog_out(DEBUG1) << "	Couldn't Find a MCParticle linked to PFO" << std::endl;
		return NULL;
	}

}

EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO )
{
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "	Look for PFO linked to visible MCParticle:" << std::endl;
	streamlog_out(DEBUG1) << *mcParticle << std::endl;
	ReconstructedParticle* linkedPFO{};
	bool foundlinkedPFO = false;
	const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( mcParticle );
	const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( mcParticle );
	streamlog_out(DEBUG0) << "	Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
	double maxweightPFOtoMCP = 0.;
	double maxweightMCPtoPFO = 0.;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
	{
		double pfo_weight = 0.0;
		double trackWeight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
		double clusterWeight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
		if ( getChargedPFO && !getNeutralPFO )
		{
			pfo_weight = trackWeight;
		}
		else if ( getNeutralPFO && !getChargedPFO )
		{
			pfo_weight = clusterWeight;
		}
		else
		{
			pfo_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
		}
		streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << trackWeight << " , Cluster: " << clusterWeight << ")" << std::endl;
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
		if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightMCPtoPFO = pfo_weight;
			iMCPtoPFOmax = i_pfo;
			streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
		}
	}
	if ( iMCPtoPFOmax != -1 )
	{
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( iMCPtoPFOmax );
		const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( testPFO );
		const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( testPFO );
		for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
		{
			double mcp_weight = 0.0;
			double trackWeight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
			double clusterWeight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
			if ( getChargedPFO && !getNeutralPFO )
			{
				mcp_weight = trackWeight;
			}
			else if ( getNeutralPFO && !getChargedPFO )
			{
				mcp_weight = clusterWeight;
			}
			else
			{
				mcp_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
			}
			MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
			if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightPFOtoMCP = mcp_weight;
				iPFOtoMCPmax = i_mcp;
				streamlog_out(DEBUG0) << "	MCParticle at index: " << testMCP->id() << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
			}
		}
		if ( iPFOtoMCPmax != -1 )
		{
			if ( MCPvec.at( iPFOtoMCPmax ) == mcParticle )
			{
				linkedPFO = testPFO;
				foundlinkedPFO = true;
			}
		}
	}

	if( foundlinkedPFO )
	{
		streamlog_out(DEBUG1) << "	Linked PFO to MCParticle found successfully " << std::endl;
		streamlog_out(DEBUG1) << *linkedPFO << std::endl;
		weightPFOtoMCP = maxweightPFOtoMCP;
		weightMCPtoPFO = maxweightMCPtoPFO;
		return linkedPFO;
	}
	else
	{
		streamlog_out(DEBUG1) << "	Couldn't Find a PFO linked to MCParticle" << std::endl;
		return NULL;
	}
}

float getWidestCosAlphaOfDecayProducts( std::vector<EVENT::MCParticle*> mcpVec , TVector3 direction )
{
	direction.SetMag( 1.0 );
	float widestCosAlpha = 1.0;
	for ( unsigned int i_mcp = 0 ; i_mcp < mcpVec.size() ; ++i_mcp )
	{
		TVector3 mcpFlightDirection( mcpVec[ i_mcp ]->getMomentum() );
		mcpFlightDirection.SetMag( 1.0 );
		if ( mcpFlightDirection.Dot( direction ) < widestCosAlpha ) widestCosAlpha = mcpFlightDirection.Dot( direction );
	}
	return widestCosAlpha;
}

float getWidestCosAlphaOfChargedPFOs( std::vector<EVENT::ReconstructedParticle*> pfoVec , TVector3 direction )
{
	direction.SetMag( 1.0 );
	float widestCosAlpha = 1.0;
	for ( unsigned int i_pfo = 0 ; i_pfo < pfoVec.size() ; ++i_pfo )
	{
		TVector3 pfoFlightDirection( pfoVec[ i_pfo ]->getMomentum() );
		pfoFlightDirection.SetMag( 1.0 );
		if ( pfoFlightDirection.Dot( direction ) < widestCosAlpha ) widestCosAlpha = pfoFlightDirection.Dot( direction );
	}
	return widestCosAlpha;
}

float getWidestCosAlphaOfVertices( std::vector<EVENT::MCParticle*> mcpVec , TVector3 direction , std::vector<float> truePrimaryVertex )
{
	float widestCosAlpha = 1.0;
	direction.SetMag( 1.0 );
	for ( unsigned int i_mcp = 0 ; i_mcp < mcpVec.size() ; ++i_mcp )
	{
		TVector3 mcpFlightDirection( mcpVec[ i_mcp ]->getEndpoint()[ 0 ] - truePrimaryVertex[ 0 ] , mcpVec[ i_mcp ]->getEndpoint()[ 1 ] - truePrimaryVertex[ 1 ] , mcpVec[ i_mcp ]->getEndpoint()[ 2 ] - truePrimaryVertex[ 2 ] );
		mcpFlightDirection.SetMag( 1.0 );
		if ( mcpFlightDirection.Dot( direction ) < widestCosAlpha ) widestCosAlpha = mcpFlightDirection.Dot( direction );
	}
	return widestCosAlpha;
}
