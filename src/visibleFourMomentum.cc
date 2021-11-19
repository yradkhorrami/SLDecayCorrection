#include "streamlog/streamlog.h"
#include "marlin/Global.h"

#include "visibleFourMomentum.h"
//#include "linkedPFO.h"

using namespace lcio;
using namespace marlin;
/*
int getLeptonFourMomentum( LCEvent *pLCEvent , MCParticle *SLDLepton  , bool cheatLepton4momentum , TLorentzVector &fourMomentumLepton , TLorentzVector &trueFourMomentumLepton , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	ReconstructedParticle* linkedRecoLepton = NULL;
	trueFourMomentumLepton = TLorentzVector( SLDLepton->getMomentum()[ 0 ] , SLDLepton->getMomentum()[ 1 ] , SLDLepton->getMomentum()[ 2 ] , SLDLepton->getEnergy() );
	if ( cheatLepton4momentum )
	{
		fourMomentumLepton = trueFourMomentumLepton;
	}
	else
	{
		linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );

		if ( linkedRecoLepton == NULL )
		{
			return -1;
		}

		else
		{
			fourMomentumLepton = TLorentzVector( linkedRecoLepton->getMomentum()[ 0 ] , linkedRecoLepton->getMomentum()[ 1 ] , linkedRecoLepton->getMomentum()[ 2 ] , linkedRecoLepton->getEnergy() );
		}
	}
	streamlog_out(DEBUG1) << "		Lepton" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << SLDLepton->getPDG() << "	, " << SLDLepton->getMass() << "	, " << SLDLepton->getMomentum()[ 0 ] << "	, " << SLDLepton->getMomentum()[ 1 ] << "	, " << SLDLepton->getMomentum()[ 2 ] << "	, " << SLDLepton->getEnergy() << "	, " << SLDLepton->getCharge() << "	)" << std::endl;
	if ( !cheatLepton4momentum && linkedRecoLepton != NULL ) streamlog_out(DEBUG1) << "			Reco:(	" << linkedRecoLepton->getType() << "	, " << linkedRecoLepton->getMass() << "	, " << linkedRecoLepton->getMomentum()[ 0 ] << "	, " << linkedRecoLepton->getMomentum()[ 1 ] << "	, " << linkedRecoLepton->getMomentum()[ 2 ] << "	, " << linkedRecoLepton->getEnergy() << "	, " << linkedRecoLepton->getCharge() << "	)" << std::endl;
	return 1;
}

int getChargedFourMomentum( LCEvent *pLCEvent , MCParticle *SLDLepton  , bool cheatCharged4momentum , TLorentzVector &fourMomentumCharged , TLorentzVector &trueFourMomentumCharged , float &restCharge , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	mcpVector chargedMCPs;
	pfoVector chargedPFOs;
//	TLorentzVector trueFourMomentumCharged( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoFourMomentumCharged( 0.0 , 0.0 , 0.0 , 0.0 );
	float totalTrueCharge = 0.0;
	float totalRecoCharge = 0.0;
	try
	{
		EVENT::MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
		getLinkedChargedMCPs( pLCEvent , SLDLepton , parentHadron , chargedMCPs );
		getLinkedChargedPFOs( pLCEvent , SLDLepton , parentHadron , chargedPFOs , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
		for ( unsigned int i_mcp = 0 ; i_mcp < chargedMCPs.size() ; ++i_mcp )
		{
			trueFourMomentumCharged += TLorentzVector( chargedMCPs[ i_mcp ]->getMomentum()[ 0 ] , chargedMCPs[ i_mcp ]->getMomentum()[ 1 ] , chargedMCPs[ i_mcp ]->getMomentum()[ 2 ] , chargedMCPs[ i_mcp ]->getEnergy() );
			totalTrueCharge += chargedMCPs[ i_mcp ]->getCharge();
		}
		for ( unsigned int i_pfo = 0 ; i_pfo < chargedPFOs.size() ; ++i_pfo )
		{
			recoFourMomentumCharged += TLorentzVector( chargedPFOs[ i_pfo ]->getMomentum()[ 0 ] , chargedPFOs[ i_pfo ]->getMomentum()[ 1 ] , chargedPFOs[ i_pfo ]->getMomentum()[ 2 ] , chargedPFOs[ i_pfo ]->getEnergy() );
			totalRecoCharge += chargedPFOs[ i_pfo ]->getCharge();
		}
		if ( cheatCharged4momentum )
		{
			fourMomentumCharged = trueFourMomentumCharged;
			restCharge = totalTrueCharge;
		}
		else
		{
			fourMomentumCharged = recoFourMomentumCharged;
			restCharge = totalRecoCharge;
		}
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Lepton for semi-leptonic decay not found" << std::endl;
        }
	streamlog_out(DEBUG1) << "	Found " << chargedMCPs.size() << " charged MCParticles for the semi-leptonic decay" << std::endl;
	streamlog_out(DEBUG1) << "	Found " << chargedPFOs.size() << " charged PFOs for the semi-leptonic decay" << std::endl;
	streamlog_out(DEBUG1) << "		Charged" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << "qqq" << "	, " << trueFourMomentumCharged.M() << "	, " << trueFourMomentumCharged.Px() << "	, " << trueFourMomentumCharged.Py() << "	, " << trueFourMomentumCharged.Pz() << "	, " << trueFourMomentumCharged.E() << "	, " << totalTrueCharge << "	)" << std::endl;
	streamlog_out(DEBUG1) << "			Reco:(	" << "qqq" << "	, " << recoFourMomentumCharged.M() << "	, " << recoFourMomentumCharged.Px() << "	, " << recoFourMomentumCharged.Py() << "	, " << recoFourMomentumCharged.Pz() << "	, " << recoFourMomentumCharged.E() << "	, " << totalRecoCharge << "	)" << std::endl;
	return 1;
}

int getLinkedChargedMCPs( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentHadron , mcpVector &chargedMCPs )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
		streamlog_out(DEBUG0) << "	Daughter[ " << i_daughter <<" ]: genStatus = " << daughter->getGeneratorStatus() << " , PDG = " << daughter->getPDG() << " , Charge = " << daughter->getCharge() << std::endl;
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && daughter != SLDLepton && fabs( daughter->getCharge() ) >= 0.1 )
			{
				chargedMCPs.push_back( daughter );
			}
		}
		else
		{
			getLinkedChargedMCPs( pLCEvent , SLDLepton , daughter , chargedMCPs );
		}
	}
	return 1;
}

int getLinkedChargedPFOs( LCEvent *pLCEvent , MCParticle *SLDLepton , MCParticle *parentHadron , pfoVector &chargedPFOs , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
		ReconstructedParticle* linkedPFO = NULL;
		streamlog_out(DEBUG0) << "	Daughter[ " << i_daughter <<" ]: genStatus = " << daughter->getGeneratorStatus() << " , PDG = " << daughter->getPDG() << " , Charge = " << daughter->getCharge() << std::endl;
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && daughter != SLDLepton && fabs( daughter->getCharge() ) >= 0.1 )
			{
				linkedPFO = getLinkedPFO( pLCEvent , daughter , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
				if ( linkedPFO != NULL ) chargedPFOs.push_back( linkedPFO );
			}
		}
		else
		{
			getLinkedChargedPFOs( pLCEvent , SLDLepton , daughter , chargedPFOs , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
		}
	}
	return 1;
}

int getNeutralFourMomentum( LCEvent *pLCEvent , MCParticle *SLDLepton  , bool cheatNeutral4momentum , TLorentzVector &fourMomentumNeutral , TLorentzVector &trueFourMomentumNeutral , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	mcpVector neutralMCPs;
	pfoVector neutralPFOs;
//	TLorentzVector trueFourMomentumNeutral( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoFourMomentumNeutral( 0.0 , 0.0 , 0.0 , 0.0 );
	try
	{
		EVENT::MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
		getLinkedNeutralMCPs( pLCEvent , SLDLepton , parentHadron , neutralMCPs );
		getLinkedNeutralPFOs( pLCEvent , SLDLepton , parentHadron , neutralPFOs , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
		for ( unsigned int i_mcp = 0 ; i_mcp < neutralMCPs.size() ; ++i_mcp )
		{
			trueFourMomentumNeutral += TLorentzVector( neutralMCPs[ i_mcp ]->getMomentum()[ 0 ] , neutralMCPs[ i_mcp ]->getMomentum()[ 1 ] , neutralMCPs[ i_mcp ]->getMomentum()[ 2 ] , neutralMCPs[ i_mcp ]->getEnergy() );
		}
		for ( unsigned int i_pfo = 0 ; i_pfo < neutralPFOs.size() ; ++i_pfo )
		{
			recoFourMomentumNeutral += TLorentzVector( neutralPFOs[ i_pfo ]->getMomentum()[ 0 ] , neutralPFOs[ i_pfo ]->getMomentum()[ 1 ] , neutralPFOs[ i_pfo ]->getMomentum()[ 2 ] , neutralPFOs[ i_pfo ]->getEnergy() );
		}
		if ( cheatNeutral4momentum )
		{
			fourMomentumNeutral = trueFourMomentumNeutral;
		}
		else
		{
			fourMomentumNeutral = recoFourMomentumNeutral;
		}
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Lepton for semi-leptonic decay not found" << std::endl;
        }
	streamlog_out(DEBUG1) << "	Found " << neutralMCPs.size() << " neutral MCParticles for the semi-leptonic decay" << std::endl;
	streamlog_out(DEBUG1) << "	Found " << neutralPFOs.size() << " neutral PFOs for the semi-leptonic decay" << std::endl;
	streamlog_out(DEBUG1) << "		Neutral" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << "nnn" << "	, " << trueFourMomentumNeutral.M() << "	, " << trueFourMomentumNeutral.Px() << "	, " << trueFourMomentumNeutral.Py() << "	, " << trueFourMomentumNeutral.Pz() << "	, " << trueFourMomentumNeutral.E() << "	, " << "0.000000" << "		)" << std::endl;
	streamlog_out(DEBUG1) << "			Reco:(	" << "nnn" << "	, " << recoFourMomentumNeutral.M() << "	, " << recoFourMomentumNeutral.Px() << "	, " << recoFourMomentumNeutral.Py() << "	, " << recoFourMomentumNeutral.Pz() << "	, " << recoFourMomentumNeutral.E() << "	, " << "0.000000" << "		)" << std::endl;
	return 1;
}

int getLinkedNeutralMCPs( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentHadron , mcpVector &neutralMCPs )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
		streamlog_out(DEBUG0) << "	Daughter[ " << i_daughter <<" ]: genStatus = " << daughter->getGeneratorStatus() << " , PDG = " << daughter->getPDG() << " , Charge = " << daughter->getCharge() << std::endl;
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && daughter != SLDLepton && fabs( daughter->getCharge() ) <= 0.1 )
			{
				neutralMCPs.push_back( daughter );
			}
		}
		else
		{
			getLinkedNeutralMCPs( pLCEvent , SLDLepton , daughter , neutralMCPs );
		}
	}
	return 1;
}

int getLinkedNeutralPFOs( LCEvent *pLCEvent , MCParticle *SLDLepton , MCParticle *parentHadron , pfoVector &neutralPFOs , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
		ReconstructedParticle* linkedPFO = NULL;
		streamlog_out(DEBUG0) << "	Daughter[ " << i_daughter <<" ]: genStatus = " << daughter->getGeneratorStatus() << " , PDG = " << daughter->getPDG() << " , Charge = " << daughter->getCharge() << std::endl;
		if ( daughter->getGeneratorStatus() == 1 )
		{
			if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && daughter != SLDLepton && fabs( daughter->getCharge() ) <= 0.1 )
			{
				linkedPFO = getLinkedPFO( pLCEvent , daughter , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , false , true );
				if ( linkedPFO != NULL ) neutralPFOs.push_back( linkedPFO );
			}
		}
		else
		{
			getLinkedNeutralPFOs( pLCEvent , SLDLepton , daughter , neutralPFOs , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
		}
	}
	return 1;
}
*/
