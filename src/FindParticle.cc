#include "FindParticle.h"

void getTrueDecayProducts( EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentMCP->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentMCP->getDaughters()[ i_daughter ];
		if ( daughter == SLDLepton || daughter == trueNeutrino ) continue;
		bool trueParticleExist = false;
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( fabs( daughter->getCharge() ) < 0.1 )
			{
				for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
				{
					if ( daughter == trueNeutralDecayProducts[ i_par ] ) trueParticleExist = true;
				}
				if ( !trueParticleExist ) trueNeutralDecayProducts.push_back( daughter );
			}
			else
			{
				for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProducts.size() ; ++i_par )
				{
					if ( daughter == trueChargedDecayProducts[ i_par ] ) trueParticleExist = true;
				}
				if ( !trueParticleExist ) trueChargedDecayProducts.push_back( daughter );
			}
		}
		else
		{
			getTrueDecayProducts( daughter , SLDLepton , trueNeutrino , trueNeutralDecayProducts , trueChargedDecayProducts );
		}
	}
}

void getTruePVADecayProducts(	EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , 
								pfo &linkedRecoLepton , float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep , 
								pfoVector &truePVANeutralDecayProducts , pfoVector &truePVAChargedDecayProducts , 
								LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentMCP->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentMCP->getDaughters()[ i_daughter ];
		ReconstructedParticle* linkedChargedPFO{};
		ReconstructedParticle* linkedNeutralPFO{};
		bool foundLinkedChargedPFO = false;
		bool foundLinkedNeutralPFO = false;
		float weightChargedPFOtoMCP = 0.0;
		float weightChargedMCPtoPFO = 0.0;
		float weightNeutralPFOtoMCP = 0.0;
		float weightNeutralMCPtoPFO = 0.0;
		bool recoParticleExist = false;
		if ( daughter == trueNeutrino )
		{
			continue;
		}
		else if ( daughter == SLDLepton )
		{
			getLinkedRecoLepton( SLDLepton , linkedRecoLepton , RecoMCParticleNav , MCParticleRecoNav , weightRecoLeptoMCLep , weightMCLeptoRecoLep );
			getTruePVADecayProducts( SLDLepton , SLDLepton , trueNeutrino , linkedRecoLepton , weightNeutralPFOtoMCP , weightNeutralMCPtoPFO , truePVANeutralDecayProducts , truePVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav );
		}
		else
		{
			linkedChargedPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , true , false , weightChargedPFOtoMCP , weightChargedMCPtoPFO , foundLinkedChargedPFO );
			linkedNeutralPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , false , true , weightNeutralPFOtoMCP , weightNeutralMCPtoPFO , foundLinkedNeutralPFO );
		}
		if ( foundLinkedChargedPFO )
		{
			for ( unsigned int i_par = 0 ; i_par < truePVAChargedDecayProducts.size() ; ++i_par )
			{
				if ( linkedChargedPFO == truePVAChargedDecayProducts[ i_par ] ) recoParticleExist = true;
			}
			for ( unsigned int i_par = 0 ; i_par < truePVANeutralDecayProducts.size() ; ++i_par )
			{
				if ( linkedChargedPFO == truePVANeutralDecayProducts[ i_par ] ) recoParticleExist = true;
			}
			if ( !recoParticleExist && linkedChargedPFO != linkedRecoLepton ) truePVAChargedDecayProducts.push_back( linkedChargedPFO );
		}
		else if ( foundLinkedNeutralPFO )
		{
			for ( unsigned int i_par = 0 ; i_par < truePVAChargedDecayProducts.size() ; ++i_par )
			{
				if ( linkedNeutralPFO == truePVAChargedDecayProducts[ i_par ] ) recoParticleExist = true;
			}
			for ( unsigned int i_par = 0 ; i_par < truePVANeutralDecayProducts.size() ; ++i_par )
			{
				if ( linkedNeutralPFO == truePVANeutralDecayProducts[ i_par ] ) recoParticleExist = true;
			}
			if ( !recoParticleExist && linkedNeutralPFO != linkedRecoLepton )
			{
				if ( ( linkedNeutralPFO->getTracks() ).size() == 0 )
				{
					truePVANeutralDecayProducts.push_back( linkedNeutralPFO );
				}
				else
				{
					truePVAChargedDecayProducts.push_back( linkedNeutralPFO );
				}
			}
		}
		else
		{
			getTruePVADecayProducts( daughter , SLDLepton , trueNeutrino , linkedRecoLepton , weightChargedPFOtoMCP , weightChargedMCPtoPFO , truePVANeutralDecayProducts , truePVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav );
		}
	}
}

void getLinkedRecoLepton(	EVENT::MCParticle *SLDLepton , pfo &linkedRecoLepton , 
							LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , 
							float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep )
{
	ReconstructedParticle* linkedRecoLep{};
	bool foundLinkedRecoLepton = false;
	linkedRecoLep = getLinkedPFO( SLDLepton , RecoMCParticleNav , MCParticleRecoNav , true , false , weightRecoLeptoMCLep , weightMCLeptoRecoLep , foundLinkedRecoLepton );
	if ( !foundLinkedRecoLepton )
	{
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( SLDLepton->getDaughters() ).size() ; ++i_daughter )
		{
			EVENT::MCParticle *daughter = SLDLepton->getDaughters()[ i_daughter ];
			getLinkedRecoLepton( daughter , linkedRecoLep , RecoMCParticleNav , MCParticleRecoNav , weightRecoLeptoMCLep , weightMCLeptoRecoLep );
		}
	}
	linkedRecoLepton = linkedRecoLep;
}

void getDecayProducts(	mcpVector trueDecayProducts , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts , 
						pfoVector &cheatedPVANeutralDecayProducts , pfoVector &cheatedPVAChargedDecayProducts , 
						LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav )
{
	for ( unsigned int i_mcp = 0 ; i_mcp < trueDecayProducts.size() ; ++i_mcp )
	{
		EVENT::MCParticle *trueDecayProduct = trueDecayProducts[ i_mcp ];
		ReconstructedParticle* linkedPFO{};
		bool foundLinkedChargedPFO = false;
		bool foundLinkedNeutralPFO = false;
		float weightChargedPFOtoMCP = 0.0;
		float weightChargedMCPtoPFO = 0.0;
		float weightNeutralPFOtoMCP = 0.0;
		float weightNeutralMCPtoPFO = 0.0;
		bool trueParticleExist = false;
		bool recoParticleExist = false;
		if ( fabs( trueDecayProduct->getCharge() ) < 0.1 )
		{
			linkedPFO = getLinkedPFO( trueDecayProduct , RecoMCParticleNav , MCParticleRecoNav , false , true , weightNeutralPFOtoMCP , weightNeutralMCPtoPFO , foundLinkedNeutralPFO );
			if ( foundLinkedNeutralPFO )
			{
				trueParticleExist = false;
				recoParticleExist = false;
				for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
				{
					if ( trueDecayProduct == trueNeutralDecayProducts[ i_par ] ) trueParticleExist = true;
				}
				if ( !trueParticleExist ) trueNeutralDecayProducts.push_back( trueDecayProduct );
				for ( unsigned int i_par = 0 ; i_par < cheatedPVANeutralDecayProducts.size() ; ++i_par )
				{
					if ( linkedPFO == cheatedPVANeutralDecayProducts[ i_par ] ) recoParticleExist = true;
				}
				if ( !recoParticleExist ) cheatedPVANeutralDecayProducts.push_back( linkedPFO );
			}
			else
			{
				investigateDecayChain(	trueDecayProduct , trueNeutralDecayProducts , trueChargedDecayProducts , cheatedPVANeutralDecayProducts , cheatedPVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav );
			}
		}
		else
		{
			linkedPFO = getLinkedPFO( trueDecayProduct , RecoMCParticleNav , MCParticleRecoNav , true , false , weightChargedPFOtoMCP , weightChargedMCPtoPFO , foundLinkedChargedPFO );
			if ( foundLinkedChargedPFO )
			{
				trueParticleExist = false;
				recoParticleExist = false;
				for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProducts.size() ; ++i_par )
				{
					if ( trueDecayProduct == trueChargedDecayProducts[ i_par ] ) trueParticleExist = true;
				}
				if ( !trueParticleExist ) trueChargedDecayProducts.push_back( trueDecayProduct );
				for ( unsigned int i_par = 0 ; i_par < cheatedPVAChargedDecayProducts.size() ; ++i_par )
				{
					if ( linkedPFO == cheatedPVAChargedDecayProducts[ i_par ] ) recoParticleExist = true;
				}
				if ( !recoParticleExist ) cheatedPVAChargedDecayProducts.push_back( linkedPFO );
			}
			else
			{
				investigateDecayChain(	trueDecayProduct , trueNeutralDecayProducts , trueChargedDecayProducts , cheatedPVANeutralDecayProducts , cheatedPVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav );
			}
		}
	}
}

void investigateDecayChain(	EVENT::MCParticle *parentMCP , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts , 
							pfoVector &cheatedPVANeutralDecayProducts , pfoVector &cheatedPVAChargedDecayProducts , 
							LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav )
{
	bool foundLinkedChargedPFO = false;
	bool foundLinkedNeutralPFO = false;
	float weightChargedPFOtoMCP = 0.0;
	float weightChargedMCPtoPFO = 0.0;
	float weightNeutralPFOtoMCP = 0.0;
	float weightNeutralMCPtoPFO = 0.0;
	bool trueParticleExist = false;
	bool recoParticleExist = false;
	int nNeutralDecayProducts = 0;
	int nChargedDecayProducts = 0;
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentMCP->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentMCP->getDaughters()[ i_daughter ];
		ReconstructedParticle* linkedPFO{};
		if ( fabs( daughter->getCharge() ) < 0.1 )
		{
			linkedPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , false , true , weightNeutralPFOtoMCP , weightNeutralMCPtoPFO , foundLinkedNeutralPFO );
			if ( foundLinkedNeutralPFO )
			{
				++nNeutralDecayProducts;
				for ( unsigned int i_par = 0 ; i_par < cheatedPVANeutralDecayProducts.size() ; ++i_par )
				{
					if ( linkedPFO == cheatedPVANeutralDecayProducts[ i_par ] ) recoParticleExist = true;
				}
				if ( !recoParticleExist ) cheatedPVANeutralDecayProducts.push_back( linkedPFO );
			}
			else
			{
				investigateDecayChain(	daughter , trueNeutralDecayProducts , trueChargedDecayProducts , cheatedPVANeutralDecayProducts , cheatedPVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav );
			}
		}
		else
		{
			linkedPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , true , false , weightChargedPFOtoMCP , weightChargedMCPtoPFO , foundLinkedChargedPFO );
			if ( foundLinkedChargedPFO )
			{
				++nChargedDecayProducts;
				for ( unsigned int i_par = 0 ; i_par < cheatedPVAChargedDecayProducts.size() ; ++i_par )
				{
					if ( linkedPFO == cheatedPVAChargedDecayProducts[ i_par ] ) recoParticleExist = true;
				}
				if ( !recoParticleExist ) cheatedPVAChargedDecayProducts.push_back( linkedPFO );
			}
			else
			{
				investigateDecayChain(	daughter , trueNeutralDecayProducts , trueChargedDecayProducts , cheatedPVANeutralDecayProducts , cheatedPVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav );
			}
		}
	}
	if ( nChargedDecayProducts == 0 && nNeutralDecayProducts == 0 )
	{
		if ( fabs( parentMCP->getCharge() ) < 0.1 )
		{
			trueParticleExist = false;
			for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
			{
				if ( parentMCP == trueNeutralDecayProducts[ i_par ] ) trueParticleExist = true;
			}
			if ( !trueParticleExist ) trueNeutralDecayProducts.push_back( parentMCP );
		}
		else
		{
			trueParticleExist = false;
			for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProducts.size() ; ++i_par )
			{
				if ( parentMCP == trueChargedDecayProducts[ i_par ] ) trueParticleExist = true;
			}
			if ( !trueParticleExist ) trueChargedDecayProducts.push_back( parentMCP );
		}
	}
	else if ( nChargedDecayProducts >= nNeutralDecayProducts )
	{
		trueParticleExist = false;
		for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProducts.size() ; ++i_par )
		{
			if ( parentMCP == trueChargedDecayProducts[ i_par ] ) trueParticleExist = true;
		}
		if ( !trueParticleExist ) trueChargedDecayProducts.push_back( parentMCP );
	}
	else
	{
		trueParticleExist = false;
		for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
		{
			if ( parentMCP == trueNeutralDecayProducts[ i_par ] ) trueParticleExist = true;
		}
		if ( !trueParticleExist ) trueNeutralDecayProducts.push_back( parentMCP );
	}
}

void getDecayProducts( 		EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , 
							mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts , 
							pfo &linkedRecoLepton , pfoVector &cheatedPVANeutralDecayProducts , pfoVector &cheatedPVAChargedDecayProducts , 
							LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , float &weightRecoLep2MCLep , float &weightMCLep2RecoLep , int generatorStatus )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentMCP->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentMCP->getDaughters()[ i_daughter ];
		if ( daughter == trueNeutrino ) continue;
		ReconstructedParticle* linkedPFO{};
		bool foundLinkedChargedPFO = false;
		bool foundLinkedNeutralPFO = false;
		float weightChargedPFOtoMCP = 0.0;
		float weightChargedMCPtoPFO = 0.0;
		float weightNeutralPFOtoMCP = 0.0;
		float weightNeutralMCPtoPFO = 0.0;
		bool trueParticleExist = false;
		bool recoParticleExist = false;
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( fabs( daughter->getCharge() ) > 0.1 )
			{
				if ( !trueParticleExist && SLDLepton != daughter ) trueChargedDecayProducts.push_back( daughter );
			}
			else
			{
				if ( !trueParticleExist ) trueNeutralDecayProducts.push_back( daughter );
			}
		}
		else
		{
			getDecayProducts( daughter , SLDLepton , trueNeutrino , trueNeutralDecayProducts , trueChargedDecayProducts , linkedRecoLepton , cheatedPVANeutralDecayProducts , cheatedPVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav , weightRecoLep2MCLep , weightMCLep2RecoLep , 1 );
		}

		if ( fabs( daughter->getCharge() ) > 0.1 )
		{
			linkedPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , true , false , weightChargedPFOtoMCP , weightChargedMCPtoPFO , foundLinkedChargedPFO );
		}
		else if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 )
		{
			linkedPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , false , true , weightNeutralPFOtoMCP , weightNeutralMCPtoPFO , foundLinkedNeutralPFO );
		}
		if ( daughter->getGeneratorStatus() == generatorStatus )
		{
			if ( foundLinkedChargedPFO )
			{
				if ( SLDLepton == daughter )
				{
					linkedRecoLepton = linkedPFO;
					weightRecoLep2MCLep = weightChargedPFOtoMCP;
					weightMCLep2RecoLep = weightChargedMCPtoPFO;
				}
				else if ( fabs( daughter->getCharge() ) > 0.1 )
				{
					for ( unsigned int i_par = 0 ; i_par < cheatedPVAChargedDecayProducts.size() ; ++i_par )
					{
						if ( linkedPFO == cheatedPVAChargedDecayProducts[ i_par ] ) recoParticleExist = true;
					}
					if ( !recoParticleExist ) cheatedPVAChargedDecayProducts.push_back( linkedPFO );
				}
			}
			else if ( foundLinkedNeutralPFO )
			{
				if ( SLDLepton == daughter )
				{
					if ( linkedRecoLepton == NULL )
					{
						linkedRecoLepton = linkedPFO;
						weightRecoLep2MCLep = weightChargedPFOtoMCP;
						weightMCLep2RecoLep = weightChargedMCPtoPFO;
					}
				}
				else
				{
					for ( unsigned int i_par = 0 ; i_par < cheatedPVANeutralDecayProducts.size() ; ++i_par )
					{
						if ( linkedPFO == cheatedPVANeutralDecayProducts[ i_par ] ) recoParticleExist = true;
					}
					if ( !recoParticleExist ) cheatedPVANeutralDecayProducts.push_back( linkedPFO );
				}
			}
		}
	}
}

EVENT::ReconstructedParticle* getLinkedRecoLepton(	EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , 
													LCRelationNavigator MCParticleRecoNav , float &weightRecoLeptoMCLep , 
													float &weightMCLeptoRecoLep , pfoVector &cheatedPVANeutralDecayProducts )
{
	ReconstructedParticle* linkedRecoLepton{};
	bool foundLinkedRecoLepton = false;
	linkedRecoLepton = getLinkedPFO( mcParticle , RecoMCParticleNav , MCParticleRecoNav , true , false , weightRecoLeptoMCLep , weightMCLeptoRecoLep , foundLinkedRecoLepton );
	if ( foundLinkedRecoLepton )
	{
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( mcParticle->getDaughters() ).size() ; ++i_daughter )
		{
			EVENT::MCParticle *daughter = mcParticle->getDaughters()[ i_daughter ];
			if ( mcParticle == daughter ) continue;
			ReconstructedParticle* linkedPFO{};
			float weightPFOtoMCP = 0.0;
			float weightMCPtoPFO = 0.0;
			bool foundLinkedPFO = false;
			bool recoParticleExist = false;
			if ( fabs( daughter->getCharge() ) < 0.1 )
			{
				linkedPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , false , true , weightPFOtoMCP , weightMCPtoPFO , foundLinkedPFO );
				if ( foundLinkedPFO )
				{
					for ( unsigned int i_par = 0 ; i_par < cheatedPVANeutralDecayProducts.size() ; ++i_par )
					{
						if ( linkedPFO == cheatedPVANeutralDecayProducts[ i_par ] ) recoParticleExist = true;
					}
					if ( !recoParticleExist ) cheatedPVANeutralDecayProducts.push_back( linkedPFO );
				}
			}
		}
	}
	else
	{
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( mcParticle->getDaughters() ).size() ; ++i_daughter )
		{
			EVENT::MCParticle *daughter = mcParticle->getDaughters()[ i_daughter ];
			if ( daughter->getPDG() == mcParticle->getPDG() )
			{
				linkedRecoLepton = getLinkedRecoLepton(	daughter , RecoMCParticleNav , MCParticleRecoNav , weightRecoLeptoMCLep , weightMCLeptoRecoLep , cheatedPVANeutralDecayProducts );
			}
			else
			{
				ReconstructedParticle* linkedPFO{};
				float weightPFOtoMCP = 0.0;
				float weightMCPtoPFO = 0.0;
				bool foundLinkedPFO = false;
				bool recoParticleExist = false;
				if ( fabs( daughter->getCharge() ) < 0.1 )
				{
					linkedPFO = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , false , true , weightPFOtoMCP , weightMCPtoPFO , foundLinkedPFO );
					if ( foundLinkedPFO )
					{
						for ( unsigned int i_par = 0 ; i_par < cheatedPVANeutralDecayProducts.size() ; ++i_par )
						{
							if ( linkedPFO == cheatedPVANeutralDecayProducts[ i_par ] ) recoParticleExist = true;
						}
						if ( !recoParticleExist ) cheatedPVANeutralDecayProducts.push_back( linkedPFO );
					}
				}
			}
		}
	}
	return linkedRecoLepton;
}

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
		if ( daughter == SLDLepton ) continue;
		if ( daughter->getGeneratorStatus() == 1 )
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
	weightPFOtoMCP = 0.0;
	weightMCPtoPFO = 0.0;
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
			pfo_weight = ( trackWeight >= clusterWeight ? trackWeight : clusterWeight );
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
	if ( getChargedPFO && maxweightMCPtoPFO < 0.8 )
	{
		streamlog_out(DEBUG1) << "	MCParticle has link weight lower than 0.8 ( " << maxweightMCPtoPFO << " ), looking for linked PFO in clusters" << std::endl;
		for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
		{
			double pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
			streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0 << " , Cluster: " << ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0 << ")" << std::endl;
			ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
			if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoPFO = pfo_weight;
				iMCPtoPFOmax = i_pfo;
				streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
			}
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

EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO , bool &foundlinkedPFO )
{
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "	Look for PFO linked to visible MCParticle:" << std::endl;
	streamlog_out(DEBUG1) << *mcParticle << std::endl;
	ReconstructedParticle* linkedPFO{};
	foundlinkedPFO = false;
	const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( mcParticle );
	const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( mcParticle );
	streamlog_out(DEBUG0) << "	Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
	weightPFOtoMCP = 0.0;
	weightMCPtoPFO = 0.0;
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
	if ( getChargedPFO && maxweightMCPtoPFO < 0.5 )
	{
		streamlog_out(DEBUG1) << "	MCParticle has link weight lower than 0.8 ( " << maxweightMCPtoPFO << " ), looking for linked PFO in clusters" << std::endl;
		maxweightMCPtoPFO = 0.0;
		for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
		{
			double pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
			streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0 << " , Cluster: " << ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0 << ")" << std::endl;
			ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
			if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoPFO = pfo_weight;
				iMCPtoPFOmax = i_pfo;
				streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
			}
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

void isMCParticleFromSLD( EVENT::MCParticle* parentHadron , EVENT::MCParticle* testMCParticle , bool &isMCPFromSLD )
{
	for ( unsigned int i_parent = 0 ; i_parent < testMCParticle->getParents().size() ; ++i_parent )
	{
		if ( parentHadron == testMCParticle->getParents()[ i_parent ] )
		{
			isMCPFromSLD = true;
		}
		else
		{
			isMCParticleFromSLD( parentHadron , testMCParticle->getParents()[ i_parent ] , isMCPFromSLD );
		}
	}
}
