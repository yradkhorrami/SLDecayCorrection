#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include "EVENT/Vertex.h"
#include "linkedPFO.h"
#include "flightDirection.h"

using namespace lcio;
using namespace marlin;

int getParentHadronFlightDirection( 	EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
					TVector3 &trueFlightDirection , TVector3 &recoFlightDirection ,
					std::string inputPrimaryVertex , std::string inputBuildUpVertex ,
					std::string inputJetCollection , int vertexingScenario ,
					std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection ,
					float &helicesDistance )
{
	int flightDirectionStatus = 0;
	std::vector<double> primaryVertex;
	std::vector<double> secondayVertex;
	std::vector<double> tritaryVertex;

	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	trueFlightDirection = TVector3( MotherHadron->getMomentum()[ 0 ] , MotherHadron->getMomentum()[ 1 ] , MotherHadron->getMomentum()[ 2 ] );
	trueFlightDirection.SetMag(1.0);

	int primaryVertexStatus = getPrimaryVertex( pLCEvent , SLDLepton , inputPrimaryVertex , false , primaryVertex );
	int secondayVertexStatus = getSecondaryVertex( pLCEvent , SLDLepton , inputPrimaryVertex , inputBuildUpVertex , false , secondayVertex , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , helicesDistance );
	flightDirectionStatus = secondayVertexStatus;
	streamlog_out(DEBUG1) << "	Primary Vertex status: " << primaryVertexStatus << std::endl;
	streamlog_out(DEBUG1) << "	Secondary Vertex status: " << secondayVertexStatus << std::endl;
	streamlog_out(DEBUG1) << "		Secondary Vertex" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << SLDLepton->getVertex()[ 0 ] << "	, " << SLDLepton->getVertex()[ 1 ] << "	, " << SLDLepton->getVertex()[ 2 ] << " 		)" << std::endl;
	if ( secondayVertexStatus == 2 || secondayVertexStatus == 4 ) streamlog_out(DEBUG1) << "			Reco:(	" << secondayVertex[ 0 ] << "	, " << secondayVertex[ 1 ] << "	, " << secondayVertex[ 2 ] << " 		)" << std::endl;

	if ( ( secondayVertexStatus == 2 || secondayVertexStatus == 4 ) && primaryVertexStatus == 2 )
	{
		if ( vertexingScenario == 1 )
		{
			recoFlightDirection = TVector3( secondayVertex[ 0 ] - primaryVertex[ 0 ] , secondayVertex[ 1 ] - primaryVertex[ 1 ] , secondayVertex[ 2 ] - primaryVertex[ 2 ] );
			recoFlightDirection.SetMag( 1.0 );
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: finding primary/secondary vertices and reconstruct flight direction" << std::endl;
		}
		else if ( vertexingScenario == 2 )
		{
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: Assigning jet axis to the flight direction of parent hadron" << std::endl;
			int jetAssigningStatus = getJetAxis( pLCEvent , SLDLepton , recoFlightDirection , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
			if ( jetAssigningStatus == 1 ) streamlog_out(DEBUG1) << "			Successfully assigned jet axis to the flight direction of parent hadron" << std::endl;
		}
		else if ( vertexingScenario == 3 )
		{
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: Assigning flight direction of leading particle in the jet to the flight direction of parent hadron" << std::endl;
			int leadingParticleStatus = getLeadingParticleFlightDirection( pLCEvent , SLDLepton , recoFlightDirection , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
			if ( leadingParticleStatus == 1 ) streamlog_out(DEBUG1) << "			Successfully assigned flight direction of leading particle in the jet to the flight direction of parent hadron" << std::endl;
		}
	}
	streamlog_out(DEBUG1) << "		Flight Direction" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << trueFlightDirection.X() << "	, " << trueFlightDirection.Y() << "	, " << trueFlightDirection.Z() << " 		)" << std::endl;
	streamlog_out(DEBUG1) << "			Reco:(	" << recoFlightDirection.X() << "	, " << recoFlightDirection.Y() << "	, " << recoFlightDirection.Z() << " 		)" << std::endl;
	recoFlightDirection.SetMag( 1.0 );
	return flightDirectionStatus;
}

int getPrimaryVertex( 	EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
			std::string inputPrimaryVertex , bool cheatVertices ,
			std::vector<double> &primaryVertex )
{
	int primaryVertexStatus = 0;

	try
	{
		LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
		streamlog_out(DEBUG0) << "	There is " << primaryVertexCollection->getNumberOfElements() << " Primary Vertex" << std::endl;
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the Primary Vertex collection" << std::endl;
		return -1;
	}

	if ( cheatVertices )
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		primaryVertex.push_back( MotherHadron->getVertex()[ 0 ] );
		primaryVertex.push_back( MotherHadron->getVertex()[ 1 ] );
		primaryVertex.push_back( MotherHadron->getVertex()[ 2 ] );
		primaryVertexStatus = 1;
		streamlog_out(DEBUG0) << "		true primary Vertex (x,y,z): 		" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	else
	{
		LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
		Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
		primaryVertex.push_back( primaryVtx->getPosition()[ 0 ] );
		primaryVertex.push_back( primaryVtx->getPosition()[ 1 ] );
		primaryVertex.push_back( primaryVtx->getPosition()[ 2 ] );
		primaryVertexStatus = 2;
		streamlog_out(DEBUG0) << "		reco primary Vertex (x,y,z): 	" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	return primaryVertexStatus;
}

int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
			std::string inputPrimaryVertex , std::string inputBuildUpVertex ,
			bool cheatVertices , std::vector<double> &secondayVertex , std::string inputJetCollection ,
			std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection ,
			float &helicesDistance )
{
	int secondayVertexStatus = -999;
	bool foundSecondaryVertex = false;
	try
	{
		LCCollection *buildUpVertexCollection = pLCEvent->getCollection( inputBuildUpVertex );
		LCCollection *jetCollection = pLCEvent->getCollection( inputJetCollection );
		streamlog_out(DEBUG0) << "	There are " << buildUpVertexCollection->getNumberOfElements() << " BuildUp Vertices and " << jetCollection->getNumberOfElements() << " jets" << std::endl;
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the BuildUp Vertex / Jet collection" << std::endl;
		secondayVertexStatus = -1;
	}
	if ( cheatVertices )
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		secondayVertex.push_back( MotherHadron->getEndpoint()[ 0 ] );
		secondayVertex.push_back( MotherHadron->getEndpoint()[ 1 ] );
		secondayVertex.push_back( MotherHadron->getEndpoint()[ 2 ] );
		foundSecondaryVertex = true;
		secondayVertexStatus = 0;
		streamlog_out(DEBUG0) << "		true secondary Vertex (x,y,z): 		" << secondayVertex[ 0 ] << "	, " << secondayVertex[ 1 ] << "	, " << secondayVertex[ 2 ] << std::endl;
	}
	else
	{
		ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
		LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
		Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
		double primaryVertex[ 3 ]{ primaryVtx->getPosition()[ 0 ] , primaryVtx->getPosition()[ 1 ] , primaryVtx->getPosition()[ 2 ] };
		LCCollection *jetCollection = pLCEvent->getCollection( inputJetCollection );
		int n_Jet = jetCollection->getNumberOfElements();
		LCCollection *buildUpVertexCollection = pLCEvent->getCollection( inputBuildUpVertex );
		int n_VTX = buildUpVertexCollection->getNumberOfElements();
		if ( linkedRecoLepton == NULL )
		{
			secondayVertexStatus = 1;
			foundSecondaryVertex = false;
			streamlog_out(DEBUG1) << "	There is no Reco Lepton in the event" << std::endl;
		}
		else
		{
			for ( int i_vtx = 0 ; i_vtx < n_VTX ; ++i_vtx )
			{
				Vertex* secondaryVtx = dynamic_cast<Vertex*>( buildUpVertexCollection->getElementAt( i_vtx ) );
				ReconstructedParticle* vertexRecoParticle = secondaryVtx->getAssociatedParticle();
				int nVertexParticles = ( vertexRecoParticle->getParticles() ).size();
				for ( int i_particle = 0 ; i_particle < nVertexParticles ; ++i_particle )
				{
					ReconstructedParticle* particle = vertexRecoParticle->getParticles()[ i_particle ];
					if ( linkedRecoLepton == particle )
					{
						secondayVertex.push_back( secondaryVtx->getPosition()[ 0 ] );
						secondayVertex.push_back( secondaryVtx->getPosition()[ 1 ] );
						secondayVertex.push_back( secondaryVtx->getPosition()[ 2 ] );
						secondayVertexStatus = 2;
						foundSecondaryVertex = true;
						streamlog_out(DEBUG1) << "	Found Reco Lepton in BuildUp Vertex" << std::endl;
					}
				}
			}
			if ( !foundSecondaryVertex )
			{
				std::vector<Vertex*> DownStreamVertices;
				for ( int i_vtx = 0 ; i_vtx < n_VTX ; ++i_vtx )
				{
					Vertex* secondaryVtx = dynamic_cast<Vertex*>( buildUpVertexCollection->getElementAt( i_vtx ) );
					ReconstructedParticle* vertexRecoParticle = secondaryVtx->getAssociatedParticle();
					int nVertexParticles = ( vertexRecoParticle->getParticles() ).size();
					float leadingEnergy = 0.0;
					ReconstructedParticle* leadingParticle = NULL;
					for ( int i_particle = 0 ; i_particle < nVertexParticles ; ++i_particle )
					{
						ReconstructedParticle* particle = vertexRecoParticle->getParticles()[ i_particle ];
						if ( particle->getEnergy() > leadingEnergy )
						{
							leadingEnergy = particle->getEnergy();
							leadingParticle = particle;
						}
					}
					ReconstructedParticle* DownStreamVertexJet = NULL;
					for ( int i_Jet = 0 ; i_Jet < n_Jet ; ++i_Jet )
					{
						ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt( i_Jet ) );
						int nParticles = ( jet->getParticles() ).size();
						for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
						{
							ReconstructedParticle* particle = jet->getParticles()[ i_particle ];
							if ( particle == leadingParticle )
							{
								for ( int i_lep = 0 ; i_lep < nParticles ; ++i_lep )
								{
									ReconstructedParticle* testLepton = jet->getParticles()[ i_lep ];
									if ( testLepton == linkedRecoLepton )
									{
										DownStreamVertices.push_back( secondaryVtx );
									}
								}
							}
						}
					}
				}
				streamlog_out(DEBUG1) << "	Found " << DownStreamVertices.size() << " DownStream Vertices in the same jet of semi-leptonic decay" << std::endl;
				if ( DownStreamVertices.size() == 0 )
				{
					secondayVertexStatus = 5;
				}
				else
				{
					float minFlightDistance = 1000000.0;
					Vertex* closetDownStreamVertex =NULL;
					bool foundDownStreamVertex = false;
					int intersectLeptonDSVertexstatus = -999;
					for ( int i_dsVTX = 0 ; i_dsVTX < DownStreamVertices.size() ; ++i_dsVTX )
					{
						Vertex* downStreamVertex = DownStreamVertices[ i_dsVTX ];
//						double primaryVertex[ 3 ]{ primaryVtx->getPosition()[ 0 ] , primaryVtx->getPosition()[ 1 ] , primaryVtx->getPosition()[ 2 ] };
						double dsVertex[ 3 ]{ downStreamVertex->getPosition()[ 0 ] , downStreamVertex->getPosition()[ 1 ] , downStreamVertex->getPosition()[ 2 ] };
						float dsDistance = std::sqrt( pow( dsVertex[ 0 ] - primaryVertex[ 0 ] , 2 ) + pow( dsVertex[ 1 ] - primaryVertex[ 1 ] , 2 ) + pow( dsVertex[ 2 ] - primaryVertex[ 2 ] , 2 ) );
						if ( dsDistance < minFlightDistance )
						{
							closetDownStreamVertex = downStreamVertex;
							minFlightDistance = dsDistance;
							foundDownStreamVertex = true;
						}
					}
					if ( foundDownStreamVertex )
					{
						intersectLeptonDSVertexstatus = intersectLeptonDSVertex( linkedRecoLepton , closetDownStreamVertex , secondayVertex , helicesDistance );
						secondayVertexStatus = 4;
					}
					else
					{
						secondayVertexStatus = 3;
					}
				}
			}
		}
	}
	return secondayVertexStatus;
}

int getJetAxis( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
		TVector3 &jetAxis , std::string inputJetCollection ,
		std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	int jetAssigningStatus = -999;
	ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
	LCCollection *jetCollection = pLCEvent->getCollection( inputJetCollection );
	int n_Jet = jetCollection->getNumberOfElements();
	for ( int i_Jet = 0 ; i_Jet < n_Jet ; ++i_Jet )
	{
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt( i_Jet ) );
		int nParticles = ( jet->getParticles() ).size();
		for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
		{
			ReconstructedParticle* particle = jet->getParticles()[ i_particle ];
			if ( particle == linkedRecoLepton )
			{
				jetAxis = TVector3( jet->getMomentum()[ 0 ] , jet->getMomentum()[ 1 ] , jet->getMomentum()[ 2 ] );
				jetAxis.SetMag( 1.0 );
				jetAssigningStatus = 1;
			}
		}
	}
	return jetAssigningStatus;
}

int getLeadingParticleFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
		TVector3 &leadingParticleFlightDirection , std::string inputJetCollection ,
		std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	int leadingParticleStatus = -999;
	ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
	LCCollection *jetCollection = pLCEvent->getCollection( inputJetCollection );
	int n_Jet = jetCollection->getNumberOfElements();
	for ( int i_Jet = 0 ; i_Jet < n_Jet ; ++i_Jet )
	{
		bool jetAssigned = false;
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt( i_Jet ) );
		int nParticles = ( jet->getParticles() ).size();
		float leadingEnergy = 0.0;
		ReconstructedParticle* leadingParticle = NULL;
		for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
		{
			ReconstructedParticle* particle = jet->getParticles()[ i_particle ];
			if ( particle->getEnergy() > leadingEnergy )
			{
				leadingParticle = particle;
				leadingEnergy = particle->getEnergy();
			}
			if ( particle == linkedRecoLepton )
			{
				jetAssigned = true;
			}
		}
		if ( jetAssigned )
		{
			leadingParticleFlightDirection = TVector3( leadingParticle->getMomentum()[ 0 ] , leadingParticle->getMomentum()[ 1 ] , leadingParticle->getMomentum()[ 2 ] );
			leadingParticleFlightDirection.SetMag( 1.0 );
			leadingParticleStatus = 1;
		}
	}
	return leadingParticleStatus;
}

int intersectLeptonDSVertex( 	EVENT::ReconstructedParticle *linkedRecoLepton ,
				EVENT::Vertex *downStreamVertex ,
				std::vector<double> &secondayVertex ,
			 	float &helicesDistance )
{
	float m_Bfield = 3.5;
	float scaleDSMomentum = -1000.0;
	float reverseScale = -1.0;
	float recoLeptonMomentum = std::sqrt( pow( linkedRecoLepton->getMomentum()[ 0 ] , 2 ) + pow( linkedRecoLepton->getMomentum()[ 1 ] , 2 ) + pow( linkedRecoLepton->getMomentum()[ 2 ] , 2 ) );
	Track *leptonTrack = linkedRecoLepton->getTracks()[ 0 ];
	HelixClass leptonHelix;
	leptonHelix.Initialize_Canonical( leptonTrack->getPhi() , leptonTrack->getD0() , leptonTrack->getZ0() , leptonTrack->getOmega() , leptonTrack->getTanLambda() , m_Bfield );
	float downStreamMomentum[ 3 ];
	downStreamMomentum[ 0 ] = scaleDSMomentum * downStreamVertex->getAssociatedParticle()->getMomentum()[ 0 ];
	downStreamMomentum[ 1 ] = scaleDSMomentum * downStreamVertex->getAssociatedParticle()->getMomentum()[ 1 ];
	downStreamMomentum[ 2 ] = scaleDSMomentum * downStreamVertex->getAssociatedParticle()->getMomentum()[ 2 ];
	float downStreamReverseMomentum[ 3 ];
	downStreamReverseMomentum[ 0 ] = reverseScale * downStreamVertex->getAssociatedParticle()->getMomentum()[ 0 ];
	downStreamReverseMomentum[ 1 ] = reverseScale * downStreamVertex->getAssociatedParticle()->getMomentum()[ 1 ];
	downStreamReverseMomentum[ 2 ] = reverseScale * downStreamVertex->getAssociatedParticle()->getMomentum()[ 2 ];
	float downStreamPosition[ 3 ];
	downStreamPosition[ 0 ] = downStreamVertex->getPosition()[ 0 ];
	downStreamPosition[ 1 ] = downStreamVertex->getPosition()[ 1 ];
	downStreamPosition[ 2 ] = downStreamVertex->getPosition()[ 2 ];
//	TVector3 dsMomentum = TVector3( downStreamMomentum[ 0 ] , downStreamMomentum[ 1 ] , downStreamMomentum[ 2 ] );
//	float dsPhi = dsMomentum.Phi();
//	float dsD0 = ( downStreamPosition[ 1 ] - downStreamMomentum[ 1 ] * downStreamPosition[ 0 ] / downStreamMomentum[ 1 ] ) * cos( dsPhi );
//	float dsOmega = 0.0;
	float downStreamCharge = -1.0 * linkedRecoLepton->getCharge(); //downStreamParameters[ 7 ];
	HelixClass downStreamHelix;
	downStreamHelix.Initialize_VP( downStreamPosition , downStreamMomentum , downStreamCharge , m_Bfield );
	float dsMomentum = std::sqrt( pow( downStreamMomentum[ 0 ] , 2 ) + pow( downStreamMomentum[ 1 ] , 2 ) + pow( downStreamMomentum[ 2 ] , 2 ) );
	float secVertexMomentum[ 3 ];
	float secVertexPosition[ 3 ];
//	float helicesDistance;
	if ( recoLeptonMomentum > dsMomentum )
	{
		helicesDistance = leptonHelix.getDistanceToHelix( &downStreamHelix , secVertexPosition , secVertexMomentum );
	}
	else
	{
		helicesDistance = downStreamHelix.getDistanceToHelix( &leptonHelix , secVertexPosition , secVertexMomentum );
	}
//	m_helicesDistance.push_back( helicesDistance );
//	h_helicesDistance->Fill( helicesDistance );
	secondayVertex.push_back( secVertexPosition[ 0 ] );
	secondayVertex.push_back( secVertexPosition[ 1 ] );
	secondayVertex.push_back( secVertexPosition[ 2 ] );
	streamlog_out(DEBUG1) << "	Down Stream Vertex at (	" << downStreamPosition[ 0 ] << "	,	" << downStreamPosition[ 1 ] << "	,	" << downStreamPosition[ 2 ] << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Down Stream Momentum: (	" << downStreamMomentum[ 0 ] << "	,	" << downStreamMomentum[ 1 ] << "	,	" << downStreamMomentum[ 2 ] << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Distance of Down Stream helix to SLDLepton helix = " << helicesDistance << " mm" << std::endl;
//	streamlog_out(DEBUG1) << "	TrueVertex (x,y,z) : ( " << SLDLepton->getVertex()[ 0 ] << " , " << SLDLepton->getVertex()[ 1 ] << " , " << SLDLepton->getVertex()[ 2 ] << " )" << std::endl;
//	streamlog_out(DEBUG1) << "	RecoVertex (x,y,z) : ( " << secVertexPosition[ 0 ] << " , " << secVertexPosition[ 1 ] << " , " << secVertexPosition[ 2 ] << " )" << std::endl;

	return 1;

}
