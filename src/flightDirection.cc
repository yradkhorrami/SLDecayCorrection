#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include "EVENT/Vertex.h"
#include "linkedPFO.h"
#include "flightDirection.h"
#include "DDMarlinCED.h"
#include <marlin/Processor.h>
#include "GeometryUtil.h"
using namespace lcio;
using namespace marlin;

int getParentHadronFlightDirection( 	EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
					TVector3 &trueFlightDirection , TVector3 &recoFlightDirection ,
					std::string inputPrimaryVertex , std::string inputBuildUpVertex ,
					std::string inputJetCollection , int vertexingScenario ,
					std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection ,
					float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent , SLDCorrection* thisProcessor )
{
	int flightDirectionStatus = 0;
	std::vector<double> primaryVertex;
	std::vector<double> secondayVertex;
	std::vector<double> tritaryVertex;



	if ( m_displayEvent )
	{
		DDMarlinCED::newEvent( thisProcessor ); // refresh
		DDMarlinCED::drawDD4hepDetector( thisProcessor->_theDetector , 0 , std::vector<std::string>{} ); // draw geometry
		DDCEDPickingHandler& pHandler = DDCEDPickingHandler::getInstance();
		pHandler.update(pLCEvent);
		ced_hit( 0.0 , 0.0 , 0.100 , 1 , 1 , 0xc765c4 );
		ced_hit( 0.0 , 0.0 , 1.000 , 1 , 1 , 0xc765c4 );
		ced_hit( 0.0 , 0.0 , 10.00 , 1 , 1 , 0xc765c4 );
	}

	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	trueFlightDirection = TVector3( MotherHadron->getMomentum()[ 0 ] , MotherHadron->getMomentum()[ 1 ] , MotherHadron->getMomentum()[ 2 ] );
	streamlog_out(DEBUG1) << "		Test 2, |flightDirection| = " << trueFlightDirection.Mag() << std::endl;
	trueFlightDirection.SetMag(1.0);

	if ( m_displayEvent )
	{
//		double bField = 3.5;
		double charge = MotherHadron->getCharge();
		double scale = 1.0;
		if ( std::fabs( charge ) <= 0.1 )
		{
			charge = 1.0;
			scale = 100.0;
		}
		streamlog_out( DEBUG1 ) << "	charge = " << charge << std::endl;
		double Px = scale * MotherHadron->getMomentumAtEndpoint()[ 0 ] ;
		double Py = scale * MotherHadron->getMomentumAtEndpoint()[ 1 ] ;
		double Pz = scale * MotherHadron->getMomentumAtEndpoint()[ 2 ] ;
		double Xs = MotherHadron->getVertex()[ 0 ];
		double Ys = MotherHadron->getVertex()[ 1 ];
		double Zs = MotherHadron->getVertex()[ 2 ];
		double Xe = MotherHadron->getEndpoint()[ 0 ];
		double Ye = MotherHadron->getEndpoint()[ 1 ];
		double Ze = MotherHadron->getEndpoint()[ 2 ];
		ced_hit( Xs , Ys , Zs , 2 , 1 , 0x279132 );
		ced_hit( Xs , Ys , Zs , 0 , 1 , 0x000000 );
		ced_hit( Xe , Ye , Ze , 2 , 1 , 0x279132 );
		ced_hit( Xe , Ye , Ze , 0 , 1 , 0x000000 );
		for ( double t = 0.0 ; t <= ( Xe - Xs ) / Px ; t += 0.005 * ( Xe - Xs ) / Px )
		{
			ced_hit( Xs + Px * t , Ys + Py * t , Zs + Pz * t , 0 , 1 , 0x000000 );
		}
	}


	int primaryVertexStatus = getPrimaryVertex( pLCEvent , SLDLepton , inputPrimaryVertex , false , primaryVertex , m_displayEvent );
	int secondayVertexStatus = getSecondaryVertex( pLCEvent , SLDLepton , inputPrimaryVertex , inputBuildUpVertex , false , secondayVertex , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , helicesDistance , SecondaryVertexPar , m_displayEvent );
	SecondaryVertexPar.push_back( MotherHadron->getEndpoint()[ 0 ] );
	SecondaryVertexPar.push_back( MotherHadron->getEndpoint()[ 1 ] );
	SecondaryVertexPar.push_back( MotherHadron->getEndpoint()[ 2 ] );
	flightDirectionStatus = secondayVertexStatus;
	streamlog_out(DEBUG1) << "	Parent Hadron Charge: " << MotherHadron->getCharge() << std::endl;
	streamlog_out(DEBUG1) << "	Primary Vertex status: " << primaryVertexStatus << std::endl;
	streamlog_out(DEBUG1) << "	Secondary Vertex status: " << secondayVertexStatus << std::endl;
	streamlog_out(DEBUG1) << "		Secondary Vertex" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << SLDLepton->getVertex()[ 0 ] << "	, " << SLDLepton->getVertex()[ 1 ] << "	, " << SLDLepton->getVertex()[ 2 ] << " 		)" << std::endl;

	if ( ( secondayVertexStatus == 2 || secondayVertexStatus == 4 || secondayVertexStatus == 5 ) && primaryVertexStatus == 2 )
	{
		streamlog_out(DEBUG1) << "			Reco:(	" << secondayVertex[ 0 ] << "	, " << secondayVertex[ 1 ] << "	, " << secondayVertex[ 2 ] << " 		)" << std::endl;
		recoFlightDirection = TVector3( secondayVertex[ 0 ] - primaryVertex[ 0 ] , secondayVertex[ 1 ] - primaryVertex[ 1 ] , secondayVertex[ 2 ] - primaryVertex[ 2 ] );
		streamlog_out(DEBUG1) << "		Test 3, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
		recoFlightDirection.SetMag( 1.0 );
		streamlog_out(DEBUG1) << "			Flight Direction Scenario: finding primary/secondary vertices and reconstruct flight direction" << std::endl;
/*
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
			recoFlightDirection.SetMag( 1.0 );
			if ( jetAssigningStatus == 1 ) streamlog_out(DEBUG1) << "			Successfully assigned jet axis to the flight direction of parent hadron" << std::endl;
		}
		else if ( vertexingScenario == 3 )
		{
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: Assigning flight direction of leading particle in the jet to the flight direction of parent hadron" << std::endl;
			int leadingParticleStatus = getLeadingParticleFlightDirection( pLCEvent , SLDLepton , recoFlightDirection , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
			recoFlightDirection.SetMag( 1.0 );
			if ( leadingParticleStatus == 1 ) streamlog_out(DEBUG1) << "			Successfully assigned flight direction of (Charged) leading particle in the jet to the flight direction of parent hadron" << std::endl;
			if ( leadingParticleStatus == -1 ) streamlog_out(DEBUG1) << "			Successfully assigned flight direction of (Neutral) leading particle in the jet to the flight direction of parent hadron" << std::endl;
		}
*/
	}
	else if ( secondayVertexStatus == 6 )
	{
		if ( vertexingScenario == 1 )
		{
			ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
			recoFlightDirection = TVector3( linkedRecoLepton->getMomentum()[ 0 ] , linkedRecoLepton->getMomentum()[ 1 ] , linkedRecoLepton->getMomentum()[ 2 ] );
			streamlog_out(DEBUG1) << "		Test 4, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
			recoFlightDirection.SetMag( 1.0 );
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: finding primary/secondary vertices and reconstruct flight direction" << std::endl;
		}
		else if ( vertexingScenario == 2 )
		{
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: Assigning jet axis to the flight direction of parent hadron" << std::endl;
			int jetAssigningStatus = getJetAxis( pLCEvent , SLDLepton , recoFlightDirection , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
			streamlog_out(DEBUG1) << "		Test 5, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
			recoFlightDirection.SetMag( 1.0 );
			if ( jetAssigningStatus == 1 ) streamlog_out(DEBUG1) << "			Successfully assigned jet axis to the flight direction of parent hadron" << std::endl;
		}
		else if ( vertexingScenario == 3 )
		{
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: Assigning flight direction of leading particle in the jet to the flight direction of parent hadron" << std::endl;
			int leadingParticleStatus = getLeadingParticleFlightDirection( pLCEvent , SLDLepton , recoFlightDirection , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
			streamlog_out(DEBUG1) << "		Test 6, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
			recoFlightDirection.SetMag( 1.0 );
			if ( leadingParticleStatus == 0 ) streamlog_out(DEBUG1) << "			Successfully assigned flight direction of photon leading particle in the jet to the flight direction of parent hadron" << std::endl;
			if ( leadingParticleStatus == 1 ) streamlog_out(DEBUG1) << "			Successfully assigned flight direction of (Charged) leading particle in the jet to the flight direction of parent hadron" << std::endl;
			if ( leadingParticleStatus == -1 ) streamlog_out(DEBUG1) << "			Successfully assigned flight direction of (Neutral) leading particle in the jet to the flight direction of parent hadron" << std::endl;
		}

	}
	if ( secondayVertexStatus == 1 || secondayVertexStatus == 6 )
	{
		SecondaryVertexPar.push_back( 0.0 );
		SecondaryVertexPar.push_back( 0.0 );
		SecondaryVertexPar.push_back( 0.0 );
	}
	else
	{
		SecondaryVertexPar.push_back( secondayVertex[ 0 ] );
		SecondaryVertexPar.push_back( secondayVertex[ 1 ] );
		SecondaryVertexPar.push_back( secondayVertex[ 2 ] );
	}
	streamlog_out(DEBUG1) << "		Flight Direction" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << trueFlightDirection.X() << "	, " << trueFlightDirection.Y() << "	, " << trueFlightDirection.Z() << " 		)" << std::endl;
	streamlog_out(DEBUG1) << "			Reco:(	" << recoFlightDirection.X() << "	, " << recoFlightDirection.Y() << "	, " << recoFlightDirection.Z() << " 		)" << std::endl;
	streamlog_out(DEBUG1) << "		Test 7, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
	recoFlightDirection.SetMag( 1.0 );
	if ( m_displayEvent )
	{
		DDMarlinCED::draw( thisProcessor , 1); // draw everything
	}
	return flightDirectionStatus;
}

int getPrimaryVertex( 	EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
			std::string inputPrimaryVertex , bool cheatVertices ,
			std::vector<double> &primaryVertex , bool m_displayEvent )
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
	if ( m_displayEvent )
	{
		ced_hit( primaryVertex[ 0 ] , primaryVertex[ 1 ] , primaryVertex[ 2 ] , 0 , 2 , 0xff7300 );
	}
	return primaryVertexStatus;
}

int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
			std::string inputPrimaryVertex , std::string inputBuildUpVertex ,
			bool cheatVertices , std::vector<double> &secondayVertex , std::string inputJetCollection ,
			std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection ,
			float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent )// , SLDCorrection* thisProcessor )
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
						if ( m_displayEvent )
						{
							ced_hit( secondayVertex[ 0 ] , secondayVertex[ 1 ] , secondayVertex[ 2 ] , 2 , 2 , 0x2b00ff );
							double bField = 3.5;//marlinutil::getBzAtOrigin();
							for ( int i_trk = 0 ; i_trk < nVertexParticles ; ++i_trk )
							{
								ReconstructedParticle* vtxRP = vertexRecoParticle->getParticles()[ i_trk ];
								int color;
								if ( linkedRecoLepton == vtxRP )
								{
									color = 0xff0000;
								}
								else
								{
									color = 0xb000e6;
								}
								Track *vtxTrack = vtxRP->getTracks()[ 0 ];
								double charge = ( vtxTrack->getOmega() > 0.0 ?  1.0 : -1.0 );
								double pT = bField * 3.0e-4 / std::fabs( vtxTrack->getOmega() );
								double Px = pT * std::cos( vtxTrack->getPhi() ) ;
								double Py = pT * std::sin( vtxTrack->getPhi() ) ;
								double Pz = pT * vtxTrack->getTanLambda() ;
								double Xs = vtxTrack->getReferencePoint()[ 0 ] - vtxTrack->getD0() * std::sin( vtxTrack->getPhi() );
								double Ys = vtxTrack->getReferencePoint()[ 1 ] + vtxTrack->getD0() * std::cos( vtxTrack->getPhi() );
								double Zs = vtxTrack->getReferencePoint()[ 2 ] + vtxTrack->getZ0();
								DDMarlinCED::drawHelix( bField , charge , Xs , Ys , Zs , Px , Py , Pz , 1 , 2 , color , 0.0 , 2100.0 , 3000.0 , 0 );
							}
						}
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
					secondayVertexStatus = 6;
					SecondaryVertexPar.push_back( 0.0 );
					SecondaryVertexPar.push_back( 0.0 );
					SecondaryVertexPar.push_back( 0.0 );
					SecondaryVertexPar.push_back( 0.0 );
					SecondaryVertexPar.push_back( 0.0 );
					SecondaryVertexPar.push_back( 0.0 );
				}
				else
				{
					float minFlightDistance = 1000000.0;
					Vertex* closetDownStreamVertex = NULL;
					bool foundDownStreamVertex = false;
					TVector3 trueDSVertex;
					TVector3 trueDSVertexMomentum;
					int trueDSVertexStatus = getTrueDownStreamVertex( SLDLepton , trueDSVertex , trueDSVertexMomentum );
					streamlog_out(DEBUG1) << "	Found " << trueDSVertexStatus << " particles at True Down Stream vertex" << std::endl;
					streamlog_out(DEBUG1) << "	True Down Stream vertex at (x,y,z) = 	( " << trueDSVertex.X() << "	,	" << trueDSVertex.Y() << "	,	" << trueDSVertex.Z() << "	)" << std::endl;
					streamlog_out(DEBUG1) << "	True Down Stream vertex Momentum (Px,Py,Pz) = 	( " << trueDSVertexMomentum.Px() << "	,	" << trueDSVertexMomentum.Py() << "	,	" << trueDSVertexMomentum.Pz() << "	)" << std::endl;
					for ( unsigned int i_dsVTX = 0 ; i_dsVTX < DownStreamVertices.size() ; ++i_dsVTX )
					{
						Vertex* downStreamVertex = DownStreamVertices[ i_dsVTX ];
						streamlog_out(DEBUG1) << "	Checking Down Stream vertex at (x,y,z) = 	( " << downStreamVertex->getPosition()[ 0 ] << "	,	" << downStreamVertex->getPosition()[ 1 ] << "	,	" << downStreamVertex->getPosition()[ 2 ] << "	)" << std::endl;
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
						streamlog_out(DEBUG1) << "	There are " << ( closetDownStreamVertex->getAssociatedParticle() )->getParticles().size() << " particle in DownStream Vertex" << std::endl;
						std::vector<double> PCAatLeptonTrack;
						std::vector<double> PCAatDownStreamLine;
						Track *leptonTrack	= linkedRecoLepton->getTracks()[ 0 ];
						TVector3 dsMomentum = TVector3( closetDownStreamVertex->getAssociatedParticle()->getMomentum()[ 0 ] , closetDownStreamVertex->getAssociatedParticle()->getMomentum()[ 1 ] , closetDownStreamVertex->getAssociatedParticle()->getMomentum()[ 2 ] );
						std::vector<double> downStreamPosition;
						downStreamPosition.push_back( closetDownStreamVertex->getPosition()[ 0 ] );
						downStreamPosition.push_back( closetDownStreamVertex->getPosition()[ 1 ] );
						downStreamPosition.push_back( closetDownStreamVertex->getPosition()[ 2 ] );


						int n_dsTracks = 0;
						for ( unsigned int i_dsPFOs = 0 ; i_dsPFOs < ( closetDownStreamVertex->getAssociatedParticle() )->getParticles().size() ; ++i_dsPFOs ) n_dsTracks += ( ( closetDownStreamVertex->getAssociatedParticle() )->getParticles()[ i_dsPFOs ] )->getTracks().size();
						if ( n_dsTracks == 3 || n_dsTracks == 5 || n_dsTracks == 7 || n_dsTracks == 9 )
						{
//							helicesDistance = intersectHelixHelix( linkedRecoLepton , dsMomentum , downStreamPosition , PCAatLeptonTrack , PCAatDownStreamLine );
							helicesDistance = intersectHelixLine( pLCEvent , leptonTrack , dsMomentum , downStreamPosition , PCAatLeptonTrack , PCAatDownStreamLine , inputPrimaryVertex );
							secondayVertexStatus = 4;
						}
						else
						{
							helicesDistance = intersectHelixLine( pLCEvent , leptonTrack , dsMomentum , downStreamPosition , PCAatLeptonTrack , PCAatDownStreamLine , inputPrimaryVertex );
							secondayVertexStatus = 5;
						}
						streamlog_out(DEBUG1) << "	Secondary Vertex on Lepton Track:(" << PCAatLeptonTrack[ 0 ] << "	,	" << PCAatLeptonTrack[ 1 ] << "	,	" << PCAatLeptonTrack[ 2 ] << "	)" << std::endl;
						streamlog_out(DEBUG1) << "	Secondary Vertex on DownStream Line:(" << PCAatDownStreamLine[ 0 ] << "	,	" << PCAatDownStreamLine[ 1 ] << "	,	" << PCAatDownStreamLine[ 2 ] << "	)" << std::endl;
						streamlog_out(DEBUG1) << "	Distance of Down Stream helix to SLDLepton helix = " << helicesDistance << " mm" << std::endl;
						secondayVertex = PCAatLeptonTrack;
						if ( m_displayEvent )
						{
							double bField = 3.5;//marlinutil::getBzAtOrigin();
							double lepcharge = ( leptonTrack->getOmega() > 0.0 ?  1.0 : -1.0 );
							double leppT = bField * 3.0e-4 / std::fabs( leptonTrack->getOmega() );
							double lepPx = leppT * std::cos( leptonTrack->getPhi() ) ;
							double lepPy = leppT * std::sin( leptonTrack->getPhi() ) ;
							double lepPz = leppT * leptonTrack->getTanLambda() ;
							double lepXs = leptonTrack->getReferencePoint()[ 0 ] - leptonTrack->getD0() * std::sin( leptonTrack->getPhi() );
							double lepYs = leptonTrack->getReferencePoint()[ 1 ] + leptonTrack->getD0() * std::cos( leptonTrack->getPhi() );
							double lepZs = leptonTrack->getReferencePoint()[ 2 ] + leptonTrack->getZ0();
							DDMarlinCED::drawHelix( bField , lepcharge , lepXs , lepYs , lepZs , lepPx , lepPy , lepPz , 1 , 2 , 0xff0000 , 0.0 , 2100.0 , 3000.0 , 0 );

							ced_hit( PCAatLeptonTrack[ 0 ] , PCAatLeptonTrack[ 1 ] , PCAatLeptonTrack[ 2 ] , 2 , 1 , 0xff0000 );
							ced_hit( PCAatLeptonTrack[ 0 ] , PCAatLeptonTrack[ 1 ] , PCAatLeptonTrack[ 2 ] , 0 , 1 , 0xff0000 );
							DDMarlinCED::drawHelix( bField , lepcharge , PCAatLeptonTrack[ 0 ] , PCAatLeptonTrack[ 1 ] , PCAatLeptonTrack[ 2 ] , lepPx , lepPy , lepPz , 1 , 2 , 0xff0000 , 0.0 , 2100.0 , 3000.0 , 0 );


							double dsCharge = -1.0 * lepcharge;
							double scale = ( secondayVertexStatus == 4 ? -1.0 : -100.0 );
							double dsPx = scale * dsMomentum[ 0 ];
							double dsPy = scale * dsMomentum[ 1 ];
							double dsPz = scale * dsMomentum[ 2 ];
							double dsXs = downStreamPosition[ 0 ];
							double dsYs = downStreamPosition[ 1 ];
							double dsZs = downStreamPosition[ 2 ];

							ced_hit( PCAatDownStreamLine[ 0 ] , PCAatDownStreamLine[ 1 ] , PCAatDownStreamLine[ 2 ] , 0 , 1 , 0x0000a8 );
							ced_hit( downStreamPosition[ 0 ] , downStreamPosition[ 1 ] , downStreamPosition[ 2 ] , 0 , 1 , 0x0000a8 );
							double vMin = -1.0 * ( dsMomentum.Px() * ( downStreamPosition[ 0 ] - PCAatDownStreamLine[ 0 ] ) + dsMomentum.Py() * ( downStreamPosition[ 1 ] - PCAatDownStreamLine[ 1 ] ) + dsMomentum.Pz() * ( downStreamPosition[ 2 ] - PCAatDownStreamLine[ 2 ] ) ) / dsMomentum.Mag2();
							double vMax = 0.0;
							double dsTestX	= downStreamPosition[ 0 ] + dsMomentum.Px() * vMin;
							double dsTestY	= downStreamPosition[ 1 ] + dsMomentum.Py() * vMin;
							double dsTestZ	= downStreamPosition[ 2 ] + dsMomentum.Pz() * vMin;
							ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
							dsTestX	= downStreamPosition[ 0 ] + dsMomentum.Px() * vMax;
							dsTestY	= downStreamPosition[ 1 ] + dsMomentum.Py() * vMax;
							dsTestZ	= downStreamPosition[ 2 ] + dsMomentum.Pz() * vMax;
							ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
							for ( double v = vMin ; v <= vMax ; v += ( vMax - vMin ) / 1000.0 )
							{
								dsTestX	= downStreamPosition[ 0 ] + dsMomentum.Px() * v;
								dsTestY	= downStreamPosition[ 1 ] + dsMomentum.Py() * v;
								dsTestZ	= downStreamPosition[ 2 ] + dsMomentum.Pz() * v;
								ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
							}

							ced_hit( PCAatDownStreamLine[ 0 ] , PCAatDownStreamLine[ 1 ] , PCAatDownStreamLine[ 2 ] , 2 , 1 , 0x0000ff );
//							ced_hit( downStreamPosition[ 0 ] , downStreamPosition[ 1 ] , downStreamPosition[ 2 ] , 2 , 1 , 0x0000ff );
							DDMarlinCED::drawHelix( bField , dsCharge , dsXs , dsYs , dsZs , dsPx , dsPy , dsPz , 1 , 2 , 0x0000ff , 0.0 , 2100.0 , 3000.0 , 0 );
							for ( unsigned int i_trk = 0 ; i_trk < ( closetDownStreamVertex->getAssociatedParticle() )->getParticles().size() ; ++i_trk )
							{
								ReconstructedParticle* vtxRP = ( closetDownStreamVertex->getAssociatedParticle() )->getParticles()[ i_trk ];
								Track *vtxTrack = vtxRP->getTracks()[ 0 ];
								double charge = ( vtxTrack->getOmega() > 0.0 ?  1.0 : -1.0 );
								double pT = bField * 3.0e-4 / std::fabs( vtxTrack->getOmega() );
								double Px = pT * std::cos( vtxTrack->getPhi() ) ;
								double Py = pT * std::sin( vtxTrack->getPhi() ) ;
								double Pz = pT * vtxTrack->getTanLambda() ;
								double Xs = downStreamPosition[ 0 ];
								double Ys = downStreamPosition[ 1 ];
								double Zs = downStreamPosition[ 2 ];
//								double Xs = vtxTrack->getReferencePoint()[ 0 ] - vtxTrack->getD0() * std::sin( vtxTrack->getPhi() );
//								double Ys = vtxTrack->getReferencePoint()[ 1 ] + vtxTrack->getD0() * std::cos( vtxTrack->getPhi() );
//								double Zs = vtxTrack->getReferencePoint()[ 2 ] + vtxTrack->getZ0();
								DDMarlinCED::drawHelix( bField , charge , Xs , Ys , Zs , Px , Py , Pz , 1 , 2 , 0xb000e6 , 0.0 , 2100.0 , 3000.0 , 0 );
							}
						}

//						secondayVertex = PCAatDownStreamLine;
					}
					else
					{
						secondayVertexStatus = 3;
					}
					SecondaryVertexPar.push_back( trueDSVertex.X() );
					SecondaryVertexPar.push_back( trueDSVertex.Y() );
					SecondaryVertexPar.push_back( trueDSVertex.Z() );
					SecondaryVertexPar.push_back( closetDownStreamVertex->getPosition()[ 0 ] );
					SecondaryVertexPar.push_back( closetDownStreamVertex->getPosition()[ 1 ] );
					SecondaryVertexPar.push_back( closetDownStreamVertex->getPosition()[ 2 ] );
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
	streamlog_out(DEBUG1) << "		looking for recoLepton " << linkedRecoLepton << " in " << n_Jet << " jets" << std::endl;
	for ( int i_Jet = 0 ; i_Jet < n_Jet ; ++i_Jet )
	{
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt( i_Jet ) );
		streamlog_out(DEBUG1) << "		Checking jet number " << i_Jet << " : " << jet << std::endl;
		int nParticles = ( jet->getParticles() ).size();
		for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
		{
			ReconstructedParticle* particle = jet->getParticles()[ i_particle ];
			streamlog_out(DEBUG1) << "		comparing particle " << particle << " in jet " << i_Jet << " with recoLespton" << std::endl;
			if ( particle == linkedRecoLepton )
			{
				streamlog_out(DEBUG1) << "		recoLepton matches eith particle " << i_particle << " : " << particle << std::endl;
				jetAxis = TVector3( jet->getMomentum()[ 0 ] , jet->getMomentum()[ 1 ] , jet->getMomentum()[ 2 ] );
				streamlog_out(DEBUG1) << "		Test 8, |flightDirection| = " << jetAxis.Mag() << std::endl;
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
			streamlog_out(DEBUG1) << "		Test 9, |flightDirection| = " << leadingParticleFlightDirection.Mag() << std::endl;
			leadingParticleFlightDirection.SetMag( 1.0 );
			if ( leadingParticle->getTracks().size() == 0 )
			{
				if ( leadingParticle->getType() == 22 )
				{
					leadingParticleStatus = 0;
				}
				else
				{
					leadingParticleStatus = -1;
				}
			}
			else
			{
				leadingParticleStatus = 1;
			}
		}
	}
	return leadingParticleStatus;
}

float intersectHelixHelix( 	EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 Momentum , std::vector<double> point ,
				std::vector<double> &PCAatLeptonTrack , std::vector<double> &PCAatDownStreamLine )
{
	float m_Bfield = 3.5;

	Track *leptonTrack = linkedRecoLepton->getTracks()[ 0 ];

	double leptonD0		= leptonTrack->getD0();
	double leptonZ0		= leptonTrack->getZ0();
	double leptonPhi0	= leptonTrack->getPhi();
	double leptonOmega	= leptonTrack->getOmega();
	double leptonTanLambda	= leptonTrack->getTanLambda();
	double xReferenceLepton	= leptonTrack->getReferencePoint()[ 0 ];
	double yReferenceLepton	= leptonTrack->getReferencePoint()[ 1 ];
	double zReferenceLepton	= leptonTrack->getReferencePoint()[ 2 ];
	while ( leptonPhi0 < 0.0 )
	{
		leptonPhi0 += 2.0 * 3.14159265359;
	}
	while ( leptonPhi0 >= 2.0 * 3.14159265359 )
	{
		leptonPhi0 -= 2.0 * 3.14159265359;
	}


	double lepcharge = ( leptonTrack->getOmega() > 0.0 ?  1.0 : -1.0 );
	double leppT = m_Bfield * 3.0e-4 / std::fabs( leptonTrack->getOmega() );
	double lepPx = leppT * std::cos( leptonTrack->getPhi() ) ;
	double lepPy = leppT * std::sin( leptonTrack->getPhi() ) ;
	double lepPz = leppT * leptonTrack->getTanLambda() ;
	double lepXs = leptonTrack->getReferencePoint()[ 0 ] - leptonTrack->getD0() * std::sin( leptonTrack->getPhi() );
	double lepYs = leptonTrack->getReferencePoint()[ 1 ] + leptonTrack->getD0() * std::cos( leptonTrack->getPhi() );
	double lepZs = leptonTrack->getReferencePoint()[ 2 ] + leptonTrack->getZ0();
	DDMarlinCED::drawHelix( m_Bfield , lepcharge , lepXs , lepYs , lepZs , lepPx , lepPy , lepPz , 2 , 1 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );

/*
	HelixClass leptonHelix;
	leptonHelix.Initialize_Canonical( leptonTrack->getPhi() , leptonTrack->getD0() , leptonTrack->getZ0() , leptonTrack->getOmega() , leptonTrack->getTanLambda() , m_Bfield );
	lepcharge = ( leptonHelix.getOmega() > 0.0 ?  1.0 : -1.0 );
	leppT = m_Bfield * 3.0e-4 / std::fabs( leptonHelix.getOmega() );
	lepPx = leppT * std::cos( leptonHelix.getPhi0() ) ;
	lepPy = leppT * std::sin( leptonHelix.getPhi0() ) ;
	lepPz = leppT * leptonHelix.getTanLambda() ;
	lepXs = leptonHelix.getReferencePoint()[ 0 ] - leptonHelix.getD0() * std::sin( leptonHelix.getPhi0() );
	lepYs = leptonHelix.getReferencePoint()[ 1 ] + leptonHelix.getD0() * std::cos( leptonHelix.getPhi0() );
	lepZs = leptonHelix.getReferencePoint()[ 2 ] + leptonHelix.getZ0();
*/

	DDMarlinCED::drawHelix( m_Bfield , lepcharge , lepXs , lepYs , lepZs , lepPx , lepPy , lepPz , 2 , 1 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );



	float downStreamMomentum[ 3 ];
	downStreamMomentum[ 0 ] = Momentum.Px();// downStreamVertex->getAssociatedParticle()->getMomentum()[ 0 ];
	downStreamMomentum[ 1 ] = Momentum.Py();// downStreamVertex->getAssociatedParticle()->getMomentum()[ 1 ];
	downStreamMomentum[ 2 ] = Momentum.Pz();// downStreamVertex->getAssociatedParticle()->getMomentum()[ 2 ];
	float downStreamPosition[ 3 ];
	downStreamPosition[ 0 ] = point[ 0 ];// downStreamVertex->getPosition()[ 0 ];
	downStreamPosition[ 1 ] = point[ 1 ];//downStreamVertex->getPosition()[ 1 ];
	downStreamPosition[ 2 ] = point[ 2 ];//downStreamVertex->getPosition()[ 2 ];
	float downStreamCharge = -1.0 * linkedRecoLepton->getCharge(); //downStreamParameters[ 7 ];
	HelixClass downStreamHelix;
	downStreamHelix.Initialize_VP( downStreamPosition , downStreamMomentum , downStreamCharge , m_Bfield );

	double downStreamD0		= downStreamHelix.getD0();
	double downStreamZ0		= downStreamHelix.getZ0();
	double downStreamPhi0		= downStreamHelix.getPhi0();
	double downStreamOmega		= downStreamHelix.getOmega();
	double downStreamTanLambda	= downStreamHelix.getTanLambda();
	double xReferenceDownStream	= 0.0;//downStreamHelix.getReferencePoint()[ 0 ];
	double yReferenceDownStream	= 0.0;//downStreamHelix.getReferencePoint()[ 1 ];
	double zReferenceDownStream	= 0.0;//downStreamHelix.getReferencePoint()[ 2 ];

	double dscharge = ( downStreamHelix.getOmega() > 0.0 ?  1.0 : -1.0 );
	double dspT = m_Bfield * 3.0e-4 / std::fabs( downStreamHelix.getOmega() );
	double dsPx = dspT * std::cos( downStreamHelix.getPhi0() ) ;
	double dsPy = dspT * std::sin( downStreamHelix.getPhi0() ) ;
	double dsPz = dspT * downStreamHelix.getTanLambda() ;
	double dsXs = xReferenceDownStream - downStreamHelix.getD0() * std::sin( downStreamHelix.getPhi0() );
	double dsYs = yReferenceDownStream + downStreamHelix.getD0() * std::cos( downStreamHelix.getPhi0() );
	double dsZs = zReferenceDownStream + downStreamHelix.getZ0();
//	double dsXs = downStreamHelix.getReferencePoint()[ 0 ] - downStreamHelix.getD0() * std::sin( downStreamHelix.getPhi0() );
//	double dsYs = downStreamHelix.getReferencePoint()[ 1 ] + downStreamHelix.getD0() * std::cos( downStreamHelix.getPhi0() );
//	double dsZs = downStreamHelix.getReferencePoint()[ 2 ] + downStreamHelix.getZ0();

	DDMarlinCED::drawHelix( m_Bfield , dscharge , point[ 0 ] , point[ 1 ] , point[ 2 ] , -dsPx , -dsPy , -dsPz , 2 , 3 , 0xff6a00 , 0.0 , 2100.0 , 3000.0 , 0 );

	double testPhi = downStreamPhi0;
	double testPointX	= xReferenceDownStream + ( 1.0 / downStreamOmega - downStreamD0 ) * std::sin( downStreamPhi0 ) - std::sin( testPhi ) / downStreamOmega;
	double testPointY	= yReferenceDownStream - ( 1.0 / downStreamOmega - downStreamD0 ) * std::cos( downStreamPhi0 ) + std::cos( testPhi ) / downStreamOmega;
	double testPointZ	= zReferenceDownStream + downStreamZ0 - downStreamPhi0 * downStreamTanLambda / std::fabs( downStreamOmega ) + testPhi * downStreamTanLambda / std::fabs( downStreamOmega );
	ced_hit( testPointX , testPointY , testPointZ , 0 , 1 , 0x000000 );
	testPhi = 4.983024;//downStreamPhi0;
	testPointX	= xReferenceDownStream + ( 1.0 / downStreamOmega - downStreamD0 ) * std::sin( downStreamPhi0 ) - std::sin( testPhi ) / downStreamOmega;
	testPointY	= yReferenceDownStream - ( 1.0 / downStreamOmega - downStreamD0 ) * std::cos( downStreamPhi0 ) + std::cos( testPhi ) / downStreamOmega;
	testPointZ	= zReferenceDownStream + downStreamZ0 - downStreamPhi0 * downStreamTanLambda / std::fabs( downStreamOmega ) + testPhi * downStreamTanLambda / std::fabs( downStreamOmega );
	ced_hit( testPointX , testPointY , testPointZ , 0 , 1 , 0x000000 );

//	2D-Function: x = leptonPhi , y = downStreamPhi

	double xMin = 0.0;
	double xMax = 2.0 * 3.14159265359;
	double yMin = 0.0;
	double yMax = 2.0 * 3.14159265359;
//	double xMin = leptonPhi0 - 3.14159265359 / 8.0;
//	double xMax = leptonPhi0 + 3.14159265359 / 8.0;
//	double yMin = downStreamPhi0 - 3.14159265359 / 8.0;;
//	double yMax = downStreamPhi0 + 3.14159265359 / 8.0;

	streamlog_out(DEBUG1) << "	Track Par. of Lepton: (	D0		,	Phi		,	Omega		,	Z0		,	TanLambda	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track Par. of Lepton: (	" << leptonD0 << "	,	" << leptonPhi0 << "	,	" << leptonOmega << "	,	" << leptonZ0 << "	,	" << leptonTanLambda << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	DownStreamVertex Position:  (	" << point[ 0 ] << "	,	" << point[ 1 ] << "	,	" << point[ 2 ] << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	DownStreamVertex Position:  (	" << dsXs << "	,	" << dsYs << "	,	" << dsZs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	DownStreamVertex Momentum:  (	" << Momentum.Px() << "	,	" << Momentum.Py() << "	,	" << Momentum.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Helix Par. of DSVertex:(D0		,	Phi		,	Omega		,	Z0		,	TanLambda	)" << std::endl;
	streamlog_out(DEBUG1) << "	Helix Par. of DSVertex:(" << downStreamD0 << "	,	" << downStreamPhi0 << "	,	" << downStreamOmega << "	,	" << downStreamZ0 << "	,	" << downStreamTanLambda << "	)" << std::endl;


	TF2 *distance = new TF2( "distance" , "std::sqrt( std::pow( ( [ 0 ] - std::sin( x ) / [ 1 ] ) - ( [ 5 ] - std::sin( y ) / [ 6 ] ) , 2 ) + std::pow( ( [ 2 ] + std::cos( x ) / [ 1 ] ) - ( [ 7 ] + std::cos( y ) / [ 6 ] ) , 2 ) + std::pow( ( [ 3 ] + x * [ 4 ] / std::fabs( [ 1 ] ) ) - ( [ 8 ] + y * [ 9 ] / std::fabs( [ 6 ] ) ) , 2 ) )" , xMin , xMax , yMin , yMax );

	distance->SetParameter( 0 , xReferenceLepton + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi0 ) );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << xReferenceLepton + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi0 ) << std::endl;
	distance->SetParameter( 1 , leptonOmega );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << leptonOmega << std::endl;
	distance->SetParameter( 2 , yReferenceLepton - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi0 ) );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << yReferenceLepton - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi0 ) << std::endl;
	distance->SetParameter( 3 , leptonZ0 + zReferenceLepton - leptonPhi0 * leptonTanLambda / std::fabs( leptonOmega ) );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << leptonZ0 + zReferenceLepton - leptonPhi0 * leptonTanLambda / std::fabs(leptonOmega) << std::endl;
	distance->SetParameter( 4 , leptonTanLambda );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << leptonTanLambda << std::endl;
	distance->SetParameter( 5 , xReferenceDownStream + ( 1.0 / downStreamOmega - downStreamD0 ) * std::sin( downStreamPhi0 ) );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] =  " << xReferenceDownStream + ( 1.0 / downStreamOmega - downStreamD0 ) * std::sin( downStreamPhi0 ) << std::endl;
	distance->SetParameter( 6 , downStreamOmega );
	streamlog_out(DEBUG1) << "	parameter [ 6 ] =  " << downStreamOmega << std::endl;
	distance->SetParameter( 7 , yReferenceDownStream - ( 1.0 / downStreamOmega - downStreamD0 ) * std::cos( downStreamPhi0 ) );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << yReferenceDownStream - ( 1.0 / downStreamOmega - downStreamD0 ) * std::cos( downStreamPhi0 ) << std::endl;
	distance->SetParameter( 8 , downStreamZ0 + zReferenceDownStream - downStreamPhi0 * downStreamTanLambda / std::fabs( downStreamOmega ) );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << downStreamZ0 + zReferenceDownStream - downStreamPhi0 * downStreamTanLambda / std::fabs( downStreamOmega ) << std::endl;
	distance->SetParameter( 9 , downStreamTanLambda );
	streamlog_out(DEBUG1) << "	parameter [ 9 ] =  " << downStreamTanLambda << std::endl;

	distance->SetRange( xMin , yMin , xMax , yMax );
	double minPhiLepton = 0.0;
	double minPhiDS = 0.0;
	double minDistance = distance->GetMinimumXY( minPhiLepton , minPhiDS );

	streamlog_out(DEBUG1) << "	Phi at min Distance on Lepton Helix = " << minPhiLepton << std::endl;
	streamlog_out(DEBUG1) << "	Phi at min Distance on Down Stream Helix = " << minPhiDS << std::endl;

	double xLeptonPCA	= xReferenceLepton + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi0 ) - std::sin( minPhiLepton ) / leptonOmega;
	double yLeptonPCA	= yReferenceLepton - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi0 ) + std::cos( minPhiLepton ) / leptonOmega;
	double zLeptonPCA	= zReferenceLepton + leptonZ0 - leptonPhi0 * leptonTanLambda / std::fabs( leptonOmega ) + minPhiLepton * leptonTanLambda / std::fabs( leptonOmega );
	streamlog_out(DEBUG1) << "	Lepton PCA (x,y,z) = 	( " << xLeptonPCA << "	,	" << yLeptonPCA << "	,	" << zLeptonPCA << "	)" << std::endl;
	PCAatLeptonTrack.push_back( xLeptonPCA );
	PCAatLeptonTrack.push_back( yLeptonPCA );
	PCAatLeptonTrack.push_back( zLeptonPCA );

	ced_hit( xLeptonPCA , yLeptonPCA , zLeptonPCA , 0 , 1 , 0x000000 );

	double xDownStreamPCA	= xReferenceDownStream + ( 1.0 / downStreamOmega - downStreamD0 ) * std::sin( downStreamPhi0 ) - std::sin( minPhiDS ) / downStreamOmega;
	double yDownStreamPCA	= yReferenceDownStream - ( 1.0 / downStreamOmega - downStreamD0 ) * std::cos( downStreamPhi0 ) + std::cos( minPhiDS ) / downStreamOmega;
	double zDownStreamPCA	= zReferenceDownStream + downStreamZ0 - downStreamPhi0 * downStreamTanLambda / std::fabs( downStreamOmega ) + minPhiDS * downStreamTanLambda / std::fabs( downStreamOmega );
	streamlog_out(DEBUG1) << "	DS PCA (x,y,z) = 	( " << xDownStreamPCA << "	,	" << yDownStreamPCA << "	,	" << zDownStreamPCA << "	)" << std::endl;
	PCAatDownStreamLine.push_back( xDownStreamPCA );
	PCAatDownStreamLine.push_back( yDownStreamPCA );
	PCAatDownStreamLine.push_back( zDownStreamPCA );

	ced_hit( xDownStreamPCA , yDownStreamPCA , zDownStreamPCA , 0 , 1 , 0xff6a00 );

	return minDistance;

}

float intersectHelixLine( EVENT::LCEvent *pLCEvent , EVENT::Track *track ,
				TVector3 Momentum , std::vector<double> point ,
				std::vector<double> &PCAatLeptonTrack ,
				std::vector<double> &PCAatDownStreamLine ,
				std::string inputPrimaryVertex )
{

	LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
	Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );

	double leptonD0		= track->getD0();
	double leptonZ0		= track->getZ0();
	double leptonPhi	= track->getPhi();
	double leptonOmega	= track->getOmega();
	double leptonTanLambda	= track->getTanLambda();
	double leptoncharge	= ( track->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference	= track->getReferencePoint()[ 0 ];
	double yReference	= track->getReferencePoint()[ 1 ];
	double zReference	= track->getReferencePoint()[ 2 ];
	double xCenter		= xReference + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi );
	double yCenter		= yReference - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi );
	double zCenter		= zReference + leptonZ0;

	while ( leptonPhi < 0.0 )
	{
		leptonPhi += 2.0 * 3.14159265359;
	}
	while ( leptonPhi >= 2.0 * 3.14159265359 )
	{
		leptonPhi -= 2.0 * 3.14159265359;
	}

	float m_Bfield = 3.5;
	double leppT = m_Bfield * 3.0e-4 / std::fabs( track->getOmega() );
	double lepPx = leppT * std::cos( track->getPhi() ) ;
	double lepPy = leppT * std::sin( track->getPhi() ) ;
	double lepPz = leppT * track->getTanLambda() ;
	double lepXs = track->getReferencePoint()[ 0 ] - track->getD0() * std::sin( track->getPhi() );
	double lepYs = track->getReferencePoint()[ 1 ] + track->getD0() * std::cos( track->getPhi() );
	double lepZs = track->getReferencePoint()[ 2 ] + track->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Lepton :	" << leptoncharge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Lepton Track ( Xc , Yc , Zc ) : (	" << xCenter << "	,	" << yCenter << "	,	" << zCenter << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Lepton Track ( Xc , Yc , Zc ) : (	" << lepXs << "	,	" << lepYs << "	,	" << lepZs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Lepton Track ( Px , Py , Pz ) : (	" << lepPx << "	,	" << lepPy << "	,	" << lepPz << "	)" << std::endl;

////////////////////////////////////////////////////////////////////////////////
//	double X = helixFreeParameter;
//
//	double xLepton		= xReference + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi ) - std::sin( X ) / leptonOmega;
//	double yLepton		= yReference - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi ) + std::cos( X ) / leptonOmega;
//	double zLepton		= zReference + leptonZ0 - ( X - leptonPhi ) * leptonTanLambda / std::fabs( leptonOmega );
////////////////////////////////////////////////////////////////////////////////

	double vMin = -1.0 * ( Momentum.Px() * ( point[ 0 ] - primaryVtx->getPosition()[ 0 ] ) + Momentum.Py() * ( point[ 1 ] - primaryVtx->getPosition()[ 1 ] ) + Momentum.Pz() * ( point[ 2 ] - primaryVtx->getPosition()[ 2 ] ) ) / Momentum.Mag2();
	double vMax = 0.0;

	streamlog_out(DEBUG1) << "	Track Par. of Lepton: (	D0		,	Phi		,	Omega		,	Z0		,	TanLambda	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track Par. of Lepton: (	" << leptonD0 << "	,	" << leptonPhi << "	,	" << leptonOmega << "	,	" << leptonZ0 << "	,	" << leptonTanLambda << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Lepton Track Ref. (x,y,z) : (	" << xReference << "	,	" << yReference << "	,	" << zReference << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	DownStreamVertex Position:  (	" << point[ 0 ] << "	,	" << point[ 1 ] << "	,	" << point[ 2 ] << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	DownStreamVertex Momentum:  (	" << Momentum.Px() << "	,	" << Momentum.Py() << "	,	" << Momentum.Pz() << "	)" << std::endl;


////////////////////////////////////////////////////////////////////////////////
//	double y = lineFreeParameter;
//
//	double xDownStream	= point[ 0 ] + Momentum.Px() * y;
//	double yDownStream	= point[ 1 ] + Momentum.Py() * y;
//	double zDownStream	= point[ 2 ] + Momentum.Pz() * y;
////////////////////////////////////////////////////////////////////////////////

	double xMin = leptonPhi - 3.14159265359 / 50.0;
	double xMax = leptonPhi + 3.14159265359 / 50.0;
	double yMin = vMin;
	double yMax = vMax;

	TF2 *distance = new TF2( "distance" , "std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - ( [ 7 ] + [ 8 ] * y ) , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - ( [ 9 ] + [ 10 ] * y ) , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * ( x - [ 6 ] ) ) - ( [ 11 ] + [ 12 ] * y ) , 2 ) )" , xMin , xMax , yMin , yMax );

	distance->SetParameter( 0 , xReference + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << xReference + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi ) << std::endl;
	distance->SetParameter( 1 , -1.0 / leptonOmega );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << -1.0 / leptonOmega << std::endl;
	distance->SetParameter( 2 , yReference - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << yReference - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi ) << std::endl;
	distance->SetParameter( 3 , 1.0 / leptonOmega );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << 1.0 / leptonOmega << std::endl;
	distance->SetParameter( 4 , zReference + leptonZ0 );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << zReference + leptonZ0 << std::endl;
	distance->SetParameter( 5 , -1.0 * leptonTanLambda / leptonOmega );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] =  " << -1.0 * leptonTanLambda / ( leptonOmega ) << std::endl;
	distance->SetParameter( 6 , leptonPhi );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << leptonPhi << std::endl;
	distance->SetParameter( 7 , point[ 0 ] );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << point[ 0 ] << std::endl;
	distance->SetParameter( 8 , Momentum.Px() );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << Momentum.Px() << std::endl;
	distance->SetParameter( 9 , point[ 1 ] );
	streamlog_out(DEBUG1) << "	parameter [ 9 ] =  " << point[ 1 ] << std::endl;
	distance->SetParameter( 10 , Momentum.Py() );
	streamlog_out(DEBUG1) << "	parameter [ 10 ] =  " << Momentum.Py() << std::endl;
	distance->SetParameter( 11 , point[ 2 ] );
	streamlog_out(DEBUG1) << "	parameter [ 11 ] =  " << point[ 2 ] << std::endl;
	distance->SetParameter( 12 , Momentum.Pz() );
	streamlog_out(DEBUG1) << "	parameter [ 12 ] = " << Momentum.Pz() << std::endl;

	distance->SetRange( xMin , yMin , xMax , yMax );
	double minPhi = 0.0;
	double minV = 0.0;
	double minDistance = distance->GetMinimumXY( minPhi , minV );

	streamlog_out(DEBUG1) << "	Phi at min Distance = " << minPhi << std::endl;
	streamlog_out(DEBUG1) << "	V at min Distance = " << minV << std::endl;

	double xLeptonPCA	= xReference + ( 1.0 / leptonOmega - leptonD0 ) * std::sin( leptonPhi ) - std::sin( minPhi ) / leptonOmega;
	double yLeptonPCA	= yReference - ( 1.0 / leptonOmega - leptonD0 ) * std::cos( leptonPhi ) + std::cos( minPhi ) / leptonOmega;
	double zLeptonPCA	= zReference + leptonZ0 - ( minPhi - leptonPhi ) * leptonTanLambda / leptonOmega;
	streamlog_out(DEBUG1) << "	Lepton PCA (x,y,z) = 	( " << xLeptonPCA << "	,	" << yLeptonPCA << "	,	" << zLeptonPCA << "	)" << std::endl;
	PCAatLeptonTrack.push_back( xLeptonPCA );
	PCAatLeptonTrack.push_back( yLeptonPCA );
	PCAatLeptonTrack.push_back( zLeptonPCA );

	double xDownStream	= point[ 0 ] + Momentum.Px() * minV;
	double yDownStream	= point[ 1 ] + Momentum.Py() * minV;
	double zDownStream	= point[ 2 ] + Momentum.Pz() * minV;
	streamlog_out(DEBUG1) << "	DS PCA (x,y,z) = 	( " << xDownStream << "	,	" << yDownStream << "	,	" << zDownStream << "	)" << std::endl;
	PCAatDownStreamLine.push_back( xDownStream );
	PCAatDownStreamLine.push_back( yDownStream );
	PCAatDownStreamLine.push_back( zDownStream );
	return minDistance;
}

int getTrueDownStreamVertex( EVENT::MCParticle *SLDLepton , TVector3 &trueDSVertex , TVector3 &trueDSVertexMomentum )
{
	EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	int nChargedParticles = getTrueVertex( MotherHadron , trueDSVertex ,trueDSVertexMomentum );
	return nChargedParticles;
}

int getTrueVertex( EVENT::MCParticle *MotherParticle , TVector3 &trueDSVertex ,TVector3 &trueDSVertexMomentum )
{
	bool foundDSVertex = false;
	int nChargedParticles = 0;
	for ( unsigned int i_daughter1 = 0 ; i_daughter1 < MotherParticle->getDaughters().size() ; ++i_daughter1 )
	{
		const EVENT::MCParticle *Daughter1 = MotherParticle->getDaughters()[ i_daughter1 ];
		if ( Daughter1->getGeneratorStatus() == 1 && Daughter1->getCharge() >= 0.9 )
		{
			for ( unsigned int i_daughter2 = 0 ; i_daughter2 < MotherParticle->getDaughters().size() ; ++i_daughter2 )
			{
				const EVENT::MCParticle *Daughter2 = MotherParticle->getDaughters()[ i_daughter2 ];
				if ( Daughter2->getGeneratorStatus() == 1 && Daughter2->getCharge() <= -0.9 )
				{
					foundDSVertex = true;
				}
			}
		}
	}
	if ( foundDSVertex )
	{
		trueDSVertex = TVector3( MotherParticle->getEndpoint()[ 0 ] , MotherParticle->getEndpoint()[ 1 ] , MotherParticle->getEndpoint()[ 2 ] );
		trueDSVertexMomentum = TVector3( MotherParticle->getMomentum()[ 0 ] , MotherParticle->getMomentum()[ 1 ] , MotherParticle->getMomentum()[ 2 ] );
		nChargedParticles = MotherParticle->getDaughters().size();
		for ( unsigned int i_daughter = 0 ; i_daughter < MotherParticle->getDaughters().size() ; ++i_daughter )
		{
			EVENT::MCParticle *Daughter = MotherParticle->getDaughters()[ i_daughter ];
			if ( Daughter->getGeneratorStatus() == 1 && fabs( Daughter->getCharge() ) >= 0.5 ) nChargedParticles++;
		}
	}
	else
	{
		for ( unsigned int i_daughter = 0 ; i_daughter < MotherParticle->getDaughters().size() ; ++i_daughter )
		{
			EVENT::MCParticle *Daughter = MotherParticle->getDaughters()[ i_daughter ];
			nChargedParticles = getTrueVertex( Daughter , trueDSVertex ,trueDSVertexMomentum );
		}
	}
	return nChargedParticles;
}
