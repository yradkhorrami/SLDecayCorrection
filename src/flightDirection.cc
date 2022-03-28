#include "flightDirection.h"
using namespace lcio;
using namespace marlin;

void getTrueFlightDirection( EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , std::vector<double> &trueStartVertex , std::vector<double> &trueSLDVertex )
{
	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	trueFlightDirection = TVector3( MotherHadron->getMomentumAtEndpoint()[ 0 ] , MotherHadron->getMomentumAtEndpoint()[ 1 ] , MotherHadron->getMomentumAtEndpoint()[ 2 ] );
	trueFlightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG1) << "	True Flight Direction ( nx , ny , nz ):	(	" << trueFlightDirection.X() << "	,	" << trueFlightDirection.Y() << "	,	" << trueFlightDirection.Z() << "	)" << std::endl;
	trueStartVertex.clear();
	trueStartVertex.push_back( MotherHadron->getVertex()[ 0 ] );
	trueStartVertex.push_back( MotherHadron->getVertex()[ 1 ] );
	trueStartVertex.push_back( MotherHadron->getVertex()[ 2 ] );
	streamlog_out(DEBUG1) << "	True Start Vertex ( x , y , z ):	(	" << trueStartVertex[ 0 ] << "	,	" << trueStartVertex[ 1 ] << "	,	" << trueStartVertex[ 2 ] << "	)" << std::endl;
	trueSLDVertex.clear();
	trueSLDVertex.push_back( MotherHadron->getEndpoint()[ 0 ] );
	trueSLDVertex.push_back( MotherHadron->getEndpoint()[ 1 ] );
	trueSLDVertex.push_back( MotherHadron->getEndpoint()[ 2 ] );
	streamlog_out(DEBUG1) << "	True End Vertex ( x , y , z ):		(	" << trueSLDVertex[ 0 ] << "	,	" << trueSLDVertex[ 1 ] << "	,	" << trueSLDVertex[ 2 ] << "	)" << std::endl;
	return;
}

int getRecoFlightDirection( 	EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 &recoFlightDirection ,
				double &hadronFlightLength , EVENT::Vertex *primaryVertex , EVENT::Vertex *startVertex ,
				vtxVector &SLDVertices , pfoVector &SLDVerticesRP , EVENT::ReconstructedParticle *assignedJet ,
				std::vector<EVENT::Vertex*> verticesInJet , std::vector<EVENT::ReconstructedParticle*> PFOswithAloneTracks ,
				float &helicesDistance , int vertexingScenario , TVector3 &daughterHadronFlightDirection ,
				double &daughterHadronFlightDistance , floatVector& sldVertexPosition )
{
	int SLDStatus = -999;
//	drawReconstructedParticle( linkedRecoLepton , primaryVertex , 0xf00000 , 0xf00000 );

	bool recoLeptonIsInVertex = false;
	TVector3 startVertexPosition( startVertex->getPosition()[ 0 ] , startVertex->getPosition()[ 1 ] , startVertex->getPosition()[ 2 ] );
	EVENT::Vertex* recoLeptonVertex = getParticleVertex( linkedRecoLepton , verticesInJet , recoLeptonIsInVertex );
	recoFlightDirection = TVector3( 0.0 , 0.0 , 0.0 );
	hadronFlightLength = 0.0;
	helicesDistance = 0.0;
	daughterHadronFlightDistance = 0.0;
	daughterHadronFlightDirection = TVector3( 0.0 , 0.0 , 0.0 );
	SLDVertices.clear();
	SLDVerticesRP.clear();
	sldVertexPosition.clear();
	if ( recoLeptonIsInVertex )
	{
		SLDStatus = 4;
		streamlog_out(DEBUG1) << "	(" << SLDStatus << ") Lepton from semi-leptonic decay found in a BuildUp Vertex, BuildUp Vertex is used as vertex of semi-leptonic decay!" << std::endl;
		SLDVertices.push_back( recoLeptonVertex );
		SLDVerticesRP.push_back( recoLeptonVertex->getAssociatedParticle() );
		recoFlightDirection = TVector3( recoLeptonVertex->getPosition()[ 0 ] - startVertex->getPosition()[ 0 ] , recoLeptonVertex->getPosition()[ 1 ] - startVertex->getPosition()[ 1 ] , recoLeptonVertex->getPosition()[ 2 ] - startVertex->getPosition()[ 2 ] );
		hadronFlightLength = recoFlightDirection.Mag();
		recoFlightDirection.SetMag( 1.0 );
		ReconstructedParticle* recoLeptonVertexRP = recoLeptonVertex->getAssociatedParticle();
		for ( unsigned int i_par = 0 ; i_par < recoLeptonVertexRP->getParticles().size() ; ++i_par )
		{
			if ( recoLeptonVertexRP->getParticles()[ i_par ] != linkedRecoLepton ) daughterHadronFlightDirection += TVector3( ( recoLeptonVertexRP->getParticles()[ i_par ] )->getMomentum() );
		}
		sldVertexPosition.push_back( recoLeptonVertex->getPosition()[ 0 ] );
		sldVertexPosition.push_back( recoLeptonVertex->getPosition()[ 1 ] );
		sldVertexPosition.push_back( recoLeptonVertex->getPosition()[ 2 ] );
		daughterHadronFlightDistance = daughterHadronFlightDirection.Mag();
		daughterHadronFlightDirection.SetMag( 1.0 );
/*
		for ( unsigned int i_vtx = 0 ; i_vtx < verticesInJet.size() ; ++i_vtx )
		{
			Vertex* testVertex = verticesInJet[ i_vtx ];
			if ( testVertex != recoLeptonVertex )
			{
				SLDVertices.push_back( testVertex );
				SLDVerticesRP.push_back( testVertex->getAssociatedParticle() );
			}
		}
*/
	}
	else if ( verticesInJet.size() != 0 )
	{
		streamlog_out(DEBUG1) << "	Lepton from semi-leptonic decay not found in a BuildUp Vertex, Investigating NuildUp vertices in jet" << std::endl;
		EVENT::Vertex* thirdVertex = NULL;
		float minDistanceToPrimaryVertex = 1000000.0;
		for ( unsigned int i_vtx = 0 ; i_vtx < verticesInJet.size() ; ++i_vtx )
		{
			EVENT::Vertex* testVertex = verticesInJet[ i_vtx ];
			float distanceToPrimaryVertex = std::sqrt( pow( testVertex->getPosition()[ 0 ] - startVertex->getPosition()[ 0 ] , 2 ) + pow( testVertex->getPosition()[ 1 ] - startVertex->getPosition()[ 1 ] , 2 ) + pow( testVertex->getPosition()[ 2 ] - startVertex->getPosition()[ 2 ] , 2 ) );
			if ( distanceToPrimaryVertex < minDistanceToPrimaryVertex )
			{
				minDistanceToPrimaryVertex = distanceToPrimaryVertex;
				thirdVertex = testVertex;
			}
		}
		if ( thirdVertex != NULL )
		{
			SLDStatus = 5;
			streamlog_out(DEBUG1) << "	(" << SLDStatus << ") There is One BuildUp Vertex in jet" << std::endl;
			streamlog_out(DEBUG1) << "		Intersection point of Lepton and other BuildUp Vertices in jet is used as vertex of semi-leptonic decay!" << std::endl;
			TVector3 PCAatTrack;
			TVector3 PCAatLine;
			Track* leptonTrack = linkedRecoLepton->getTracks()[ 0 ];
			TVector3 momentumOfLine( ( thirdVertex->getAssociatedParticle() )->getMomentum() );
			std::vector<double> pointOnLine{};
			pointOnLine.push_back( thirdVertex->getPosition()[ 0 ] );
			pointOnLine.push_back( thirdVertex->getPosition()[ 1 ] );
			pointOnLine.push_back( thirdVertex->getPosition()[ 2 ] );
			helicesDistance = intersectTrackLine( leptonTrack , primaryVertex , momentumOfLine , pointOnLine , PCAatTrack , PCAatLine );
			recoFlightDirection = PCAatTrack - startVertexPosition;
			hadronFlightLength = recoFlightDirection.Mag();
			daughterHadronFlightDirection = TVector3( thirdVertex->getPosition()[ 0 ] - PCAatTrack[ 0 ] , thirdVertex->getPosition()[ 1 ] - PCAatTrack[ 1 ] , thirdVertex->getPosition()[ 2 ] - PCAatTrack[ 2 ] );
			daughterHadronFlightDistance = daughterHadronFlightDirection.Mag();
			daughterHadronFlightDirection.SetMag( 1.0 );
//			recoFlightDirection = TVector3( PCAatTrack[ 0 ] - startVertex->getPosition()[ 0 ] , PCAatTrack[ 1 ] - startVertex->getPosition()[ 1 ] , PCAatTrack[ 2 ] - startVertex->getPosition()[ 2 ] );
			recoFlightDirection.SetMag( 1.0 );
			SLDVertices.push_back( thirdVertex );
			SLDVerticesRP.push_back( thirdVertex->getAssociatedParticle() );
/*
			for ( unsigned int i_vtx = 0 ; i_vtx < verticesInJet.size() ; ++i_vtx )
			{
				Vertex* testVertex = verticesInJet[ i_vtx ];
				if ( testVertex != thirdVertex )
				{
					SLDVertices.push_back( testVertex );
					SLDVerticesRP.push_back( testVertex->getAssociatedParticle() );
				}
			}
*/
			sldVertexPosition.push_back( PCAatTrack.X() );
			sldVertexPosition.push_back( PCAatTrack.Y() );
			sldVertexPosition.push_back( PCAatTrack.Z() );
		}
	}
	else
	{
		streamlog_out(DEBUG1) << "	There is NO BuildUp Vertex in jet" << std::endl;
		streamlog_out(DEBUG1) << "	Other scenarios are investigated, but Nu-correction will not be performed" << std::endl;
		TVector3 PCAatLepton;
		TVector3 PCAatOtherParticle;
		std::vector<double> chargedParticlePosition( 3 , 0.0 );
		TVector3 chargedParticleMomentum( 0.0 , 0.0 , 0.0 );
		ReconstructedParticle* leadingChargedParticle = getLeadingChargedParticle( assignedJet );
		if ( leadingChargedParticle != NULL )
		{
			chargedParticleMomentum = TVector3( leadingChargedParticle->getMomentum() );
			chargedParticlePosition[ 0 ] = leadingChargedParticle->getReferencePoint()[ 0 ];
			chargedParticlePosition[ 1 ] = leadingChargedParticle->getReferencePoint()[ 1 ];
			chargedParticlePosition[ 2 ] = leadingChargedParticle->getReferencePoint()[ 2 ];
		}
		std::vector<double> neutralParticlePosition( 3 , 0.0 );
		TVector3 neutralParticleMomentum( 0.0 , 0.0 , 0.0 );
		ReconstructedParticle* leadingNeutralParticle = getLeadingNeutralParticle( assignedJet );
		if ( leadingNeutralParticle != NULL )
		{
			neutralParticleMomentum = TVector3( leadingNeutralParticle->getMomentum() );
			neutralParticlePosition[ 0 ] = leadingNeutralParticle->getReferencePoint()[ 0 ];
			neutralParticlePosition[ 1 ] = leadingNeutralParticle->getReferencePoint()[ 1 ];
			neutralParticlePosition[ 1 ] = leadingNeutralParticle->getReferencePoint()[ 2 ];
		}
		if ( vertexingScenario == 2 )
		{
			if ( leadingChargedParticle != NULL )
			{
				SLDStatus = 7;
				if ( leadingChargedParticle->getTracks().size() == 1 )
				{
					helicesDistance = intersectTrackTrack( linkedRecoLepton->getTracks()[ 0 ] , leadingChargedParticle->getTracks()[ 0 ] , PCAatLepton , PCAatOtherParticle );
				}
				else
				{
					helicesDistance = intersectTrackLine( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , chargedParticleMomentum , chargedParticlePosition , PCAatLepton , PCAatOtherParticle );
				}
				recoFlightDirection = PCAatLepton - startVertexPosition;
				hadronFlightLength = recoFlightDirection.Mag();
				sldVertexPosition.push_back( PCAatLepton.X() );
				sldVertexPosition.push_back( PCAatLepton.Y() );
				sldVertexPosition.push_back( PCAatLepton.Z() );
//				recoFlightDirection = TVector3( PCAatLepton[ 0 ] - startVertex->getPosition()[ 0 ] , PCAatLepton[ 1 ] - startVertex->getPosition()[ 1 ] , PCAatLepton[ 2 ] - startVertex->getPosition()[ 2 ] );
			}
			else
			{
				vertexingScenario = 1;
			}
		}
		else if ( vertexingScenario == 3 )
		{
			if ( leadingNeutralParticle != NULL )
			{
				SLDStatus = 7;
				helicesDistance = intersectTrackLine( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , neutralParticleMomentum , neutralParticlePosition , PCAatLepton , PCAatOtherParticle );
				recoFlightDirection = PCAatLepton - startVertexPosition;
				hadronFlightLength = recoFlightDirection.Mag();
				sldVertexPosition.push_back( PCAatLepton.X() );
				sldVertexPosition.push_back( PCAatLepton.Y() );
				sldVertexPosition.push_back( PCAatLepton.Z() );
//					recoFlightDirection = TVector3( PCAatLepton[ 0 ] - startVertex->getPosition()[ 0 ] , PCAatLepton[ 1 ] - startVertex->getPosition()[ 1 ] , PCAatLepton[ 2 ] - startVertex->getPosition()[ 2 ] );
			}
			else
			{
				vertexingScenario = 1;
			}
		}
		else if ( vertexingScenario == 4  )
		{
			SLDStatus = 6;
			helicesDistance = 10000.0;
			double temphelicesDistance = 0.0;
			for ( unsigned int i_par = 0 ; i_par < PFOswithAloneTracks.size() ; ++i_par )
			{
				TVector3 tempPCAatLepton;
				TVector3 tempPCAatOtherParticle;
				ReconstructedParticle* alonePFO = PFOswithAloneTracks[ i_par ];
				TVector3 alonePFOMomentum( alonePFO->getMomentum() );
				std::vector<double> alonePFOPosition{};
				alonePFOPosition.push_back( alonePFO->getReferencePoint()[ 0 ] );
				alonePFOPosition.push_back( alonePFO->getReferencePoint()[ 1 ] );
				alonePFOPosition.push_back( alonePFO->getReferencePoint()[ 2 ] );
				if ( alonePFO->getTracks().size() == 1 )
				{
					temphelicesDistance = intersectTrackTrack( linkedRecoLepton->getTracks()[ 0 ] , alonePFO->getTracks()[ 0 ] , tempPCAatLepton , tempPCAatOtherParticle );
				}
				else
				{
					temphelicesDistance = intersectTrackLine( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , alonePFOMomentum , alonePFOPosition , tempPCAatLepton , tempPCAatOtherParticle );
				}
				if ( temphelicesDistance < helicesDistance )
				{
					helicesDistance = temphelicesDistance;
					PCAatLepton = tempPCAatLepton;
					PCAatOtherParticle = tempPCAatOtherParticle;
				}
			}
			recoFlightDirection = PCAatLepton - startVertexPosition;
			hadronFlightLength = recoFlightDirection.Mag();
			sldVertexPosition.push_back( PCAatLepton.X() );
			sldVertexPosition.push_back( PCAatLepton.Y() );
			sldVertexPosition.push_back( PCAatLepton.Z() );
//			recoFlightDirection = TVector3( PCAatLepton[ 0 ] - startVertex->getPosition()[ 0 ] , PCAatLepton[ 1 ] - startVertex->getPosition()[ 1 ] , PCAatLepton[ 2 ] - startVertex->getPosition()[ 2 ] );
		}
		if ( vertexingScenario == 1 )
		{
			recoFlightDirection = TVector3( assignedJet->getMomentum() );
			hadronFlightLength = -1.0;
			SLDStatus = 7;
		}
		daughterHadronFlightDirection = recoFlightDirection;
		daughterHadronFlightDistance = daughterHadronFlightDirection.Mag();
		recoFlightDirection.SetMag( 1.0 );
		daughterHadronFlightDirection.SetMag( 1.0 );
	}
	if ( SLDVertices.size() != 0 )
	{
		streamlog_out(DEBUG2) << "Lepton is in this vertex:" << std::endl;
		streamlog_out(DEBUG2) << *SLDVertices[ 0 ] << std::endl;
	}
	recoFlightDirection.SetMag( 1.0 );
	return SLDStatus;
}

double intersectTrackLine( 	EVENT::Track *track , EVENT::Vertex* primaryVertex ,
				TVector3 momentumOfLine , std::vector<double> pointOnLine ,
				TVector3 &PCAatTrack ,
				TVector3 &PCAatLine )
{
	streamlog_out(DEBUG1) << "<<<<<<<<<<<<<<<<<<<<   Intersecting a Track and a Line   >>>>>>>>>>>>>>>>>>>>" << std::endl;
	streamlog_out(DEBUG1) << "	TRACK:" << std::endl;
	streamlog_out(DEBUG1) << *track << std::endl;
	double trackD0		= track->getD0();
	double trackZ0		= track->getZ0();
	double trackPhi		= track->getPhi();
	double trackOmega	= track->getOmega();
	double trackTanLambda	= track->getTanLambda();
	double trackcharge	= ( track->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference	= track->getReferencePoint()[ 0 ];
	double yReference	= track->getReferencePoint()[ 1 ];
	double zReference	= track->getReferencePoint()[ 2 ];
	double xCenter		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi );
	double yCenter		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi );
	double zCenter		= zReference + trackZ0;

	while ( trackPhi < 0.0 )
	{
		trackPhi += 2.0 * 3.14159265359;
	}
	while ( trackPhi >= 2.0 * 3.14159265359 )
	{
		trackPhi -= 2.0 * 3.14159265359;
	}

	double m_Bfield = MarlinUtil::getBzAtOrigin();
	double trackPt = m_Bfield * 3.0e-4 / std::fabs( track->getOmega() );
	double trackPx = trackPt * std::cos( track->getPhi() ) ;
	double trackPy = trackPt * std::sin( track->getPhi() ) ;
	double trackPz = trackPt * track->getTanLambda() ;
	double trackXs = track->getReferencePoint()[ 0 ] - track->getD0() * std::sin( track->getPhi() );
	double trackYs = track->getReferencePoint()[ 1 ] + track->getD0() * std::cos( track->getPhi() );
	double trackZs = track->getReferencePoint()[ 2 ] + track->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Track :	" << trackcharge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Track ( Xc , Yc , Zc ) :	(	" << xCenter << "	,	" << yCenter << "	,	" << zCenter << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Track ( Xc , Yc , Zc ) :	(	" << trackXs << "	,	" << trackYs << "	,	" << trackZs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Track ( Px , Py , Pz ) :	(	" << trackPx << "	,	" << trackPy << "	,	" << trackPz << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track Reference Point (x,y,z) :		(	" << xReference << "	,	" << yReference << "	,	" << zReference << "	)" << std::endl;

////////////////////////////////////////////////////////////////////////////////
//	double X = helixFreeParameter;
//
//	double xTrack		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( X ) / trackOmega;
//	double yTrack		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( X ) / trackOmega;
//	double zTrack		= zReference + trackZ0 - ( X - trackPhi ) * trackTanLambda / std::fabs( trackOmega );
////////////////////////////////////////////////////////////////////////////////

	double vMin = -1.0 * ( momentumOfLine.Px() * ( pointOnLine[ 0 ] - primaryVertex->getPosition()[ 0 ] ) + momentumOfLine.Py() * ( pointOnLine[ 1 ] - primaryVertex->getPosition()[ 1 ] ) + momentumOfLine.Pz() * ( pointOnLine[ 2 ] - primaryVertex->getPosition()[ 2 ] ) ) / momentumOfLine.Mag2();
	double vMax = 0.0;

	streamlog_out(DEBUG1) << "	DownStreamVertex Position:  (	" << pointOnLine[ 0 ] << "	,	" << pointOnLine[ 1 ] << "	,	" << pointOnLine[ 2 ] << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	DownStreamVertex Momentum:  (	" << momentumOfLine.Px() << "	,	" << momentumOfLine.Py() << "	,	" << momentumOfLine.Pz() << "	)" << std::endl;


////////////////////////////////////////////////////////////////////////////////
//	double y = lineFreeParameter;
//
//	double xDownStream	= point[ 0 ] + Momentum.Px() * y;
//	double yDownStream	= point[ 1 ] + Momentum.Py() * y;
//	double zDownStream	= point[ 2 ] + Momentum.Pz() * y;
////////////////////////////////////////////////////////////////////////////////

	double xMin = trackPhi - 3.14159265359 / 40.0;
	double xMax = trackPhi + 3.14159265359 / 40.0;
	double yMin = vMin;
	double yMax = vMax;

	DDMarlinCED::drawHelix( m_Bfield , trackcharge , trackXs , trackYs , trackZs , trackPx , trackPy , trackPz , 2 , 1 , 0xff0000 , 0.0 , 1500.0 , 2000.0 , 0 );
	double t = 20.0;
	double endPointX = pointOnLine[ 0 ] + t * momentumOfLine.Px();
	double endPointY = pointOnLine[ 1 ] + t * momentumOfLine.Py();
	double endPointZ = pointOnLine[ 2 ] + t * momentumOfLine.Pz();
	ced_line( endPointX , endPointY , endPointZ , pointOnLine[ 0 ] , pointOnLine[ 1 ] , pointOnLine[ 2 ] , 2 , 1 , 0x0061ff );

	double dphi = ( xMax - xMin ) / 100.0;
	double dv = ( yMax - yMin ) / 10.0;
	for ( double phi = xMin ; phi <= xMax ; phi += dphi )
	{
		double trackx = xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( phi ) / trackOmega;
		double tracky = yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( phi ) / trackOmega;
		double trackz = zReference + trackZ0 - ( phi - trackPhi ) * trackTanLambda / trackOmega;
		ced_hit( trackx , tracky , trackz , 2 , 3 , 0x000000 );
	}
	for ( double v = yMin ; v <= yMax ; v += dv )
	{
		double linex = pointOnLine[ 0 ] + v * momentumOfLine.Px();
		double liney = pointOnLine[ 1 ] + v * momentumOfLine.Py();
		double linez = pointOnLine[ 2 ] + v * momentumOfLine.Pz();
		ced_hit( linex , liney , linez , 2 , 3 , 0x000000 );

	}



	TF2 *distance = new TF2( "distance" , "std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - ( [ 7 ] + [ 8 ] * y ) , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - ( [ 9 ] + [ 10 ] * y ) , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * ( x - [ 6 ] ) ) - ( [ 11 ] + [ 12 ] * y ) , 2 ) )" , xMin , xMax , yMin , yMax );

	distance->SetParameter( 0 , xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) << std::endl;
	distance->SetParameter( 1 , -1.0 / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << -1.0 / trackOmega << std::endl;
	distance->SetParameter( 2 , yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) << std::endl;
	distance->SetParameter( 3 , 1.0 / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << 1.0 / trackOmega << std::endl;
	distance->SetParameter( 4 , zReference + trackZ0 );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << zReference + trackZ0 << std::endl;
	distance->SetParameter( 5 , -1.0 * trackTanLambda / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] =  " << -1.0 * trackTanLambda / ( trackOmega ) << std::endl;
	distance->SetParameter( 6 , trackPhi );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << trackPhi << std::endl;
	distance->SetParameter( 7 , pointOnLine[ 0 ] );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << pointOnLine[ 0 ] << std::endl;
	distance->SetParameter( 8 , momentumOfLine.Px() );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << momentumOfLine.Px() << std::endl;
	distance->SetParameter( 9 , pointOnLine[ 1 ] );
	streamlog_out(DEBUG1) << "	parameter [ 9 ] =  " << pointOnLine[ 1 ] << std::endl;
	distance->SetParameter( 10 , momentumOfLine.Py() );
	streamlog_out(DEBUG1) << "	parameter [ 10 ] =  " << momentumOfLine.Py() << std::endl;
	distance->SetParameter( 11 , pointOnLine[ 2 ] );
	streamlog_out(DEBUG1) << "	parameter [ 11 ] =  " << pointOnLine[ 2 ] << std::endl;
	distance->SetParameter( 12 , momentumOfLine.Pz() );
	streamlog_out(DEBUG1) << "	parameter [ 12 ] = " << momentumOfLine.Pz() << std::endl;

	distance->SetRange( xMin , yMin , xMax , yMax );
//	distance->SetNpx( 1000 );
//	distance->SetNpy( 1000 );
	double minPhi = trackPhi;
	double minV = 0.0;
	double minDistance = distance->GetMinimumXY( minPhi , minV );

	streamlog_out(DEBUG1) << "	Phi at min Distance = " << minPhi << std::endl;
	streamlog_out(DEBUG1) << "	V at min Distance = " << minV << std::endl;

	double xTrackPCA	= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( minPhi ) / trackOmega;
	double yTrackPCA	= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( minPhi ) / trackOmega;
	double zTrackPCA	= zReference + trackZ0 - ( minPhi - trackPhi ) * trackTanLambda / trackOmega;
	streamlog_out(DEBUG1) << "	Lepton PCA (x,y,z) = 	( " << xTrackPCA << "	,	" << yTrackPCA << "	,	" << zTrackPCA << "	)" << std::endl;
	ced_hit( xTrackPCA , yTrackPCA , zTrackPCA , 5 , 5 , 0xff0000 );
	PCAatTrack = TVector3( xTrackPCA , yTrackPCA , zTrackPCA );

	double xLine	= pointOnLine[ 0 ] + momentumOfLine.Px() * minV;
	double yLine	= pointOnLine[ 1 ] + momentumOfLine.Py() * minV;
	double zLine	= pointOnLine[ 2 ] + momentumOfLine.Pz() * minV;
	streamlog_out(DEBUG1) << "	DS PCA (x,y,z) = 	( " << xLine << "	,	" << yLine << "	,	" << zLine << "	)" << std::endl;
	ced_hit( xLine , yLine , zLine , 5 , 5 , 0x0061ff );
	PCAatLine = TVector3( xLine , yLine , zLine );
	streamlog_out(DEBUG1) << "	Distance of Track and Line = " << minDistance << std::endl;
	return minDistance;
}

double intersectTrackTrack(	EVENT::Track *track1 , EVENT::Track *track2 ,
				TVector3 &PCAatTrack1 ,
				TVector3 &PCAatTrack2 )
{
	double m_Bfield = MarlinUtil::getBzAtOrigin();
	streamlog_out(DEBUG1) << "<<<<<<<<<<<<<<<<<<<<<<<<   Intersecting Two Tracks   >>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	streamlog_out(DEBUG1) << "	TRACK1:" << std::endl;
	streamlog_out(DEBUG1) << *track1 << std::endl;
	double track1D0		= track1->getD0();
	double track1Z0		= track1->getZ0();
	double track1Phi	= track1->getPhi();
	double track1Omega	= track1->getOmega();
	double track1TanLambda	= track1->getTanLambda();
	double track1charge	= ( track1->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference1	= track1->getReferencePoint()[ 0 ];
	double yReference1	= track1->getReferencePoint()[ 1 ];
	double zReference1	= track1->getReferencePoint()[ 2 ];
	double xCenter1		= xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi );
	double yCenter1		= yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi );
	double zCenter1		= zReference1 + track1Z0;

	while ( track1Phi < 0.0 )
	{
		track1Phi += 2.0 * 3.14159265359;
	}
	while ( track1Phi >= 2.0 * 3.14159265359 )
	{
		track1Phi -= 2.0 * 3.14159265359;
	}

	double track1Pt = m_Bfield * 3.0e-4 / std::fabs( track1->getOmega() );
	double track1Px = track1Pt * std::cos( track1->getPhi() ) ;
	double track1Py = track1Pt * std::sin( track1->getPhi() ) ;
	double track1Pz = track1Pt * track1->getTanLambda() ;
	double track1Xs = track1->getReferencePoint()[ 0 ] - track1->getD0() * std::sin( track1->getPhi() );
	double track1Ys = track1->getReferencePoint()[ 1 ] + track1->getD0() * std::cos( track1->getPhi() );
	double track1Zs = track1->getReferencePoint()[ 2 ] + track1->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Track1 :	" << track1charge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Track1 ( Xc , Yc , Zc ) :	(	" << xCenter1 << "	,	" << yCenter1 << "	,	" << zCenter1 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Track1 ( Xc , Yc , Zc ) :	(	" << track1Xs << "	,	" << track1Ys << "	,	" << track1Zs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Track1 ( Px , Py , Pz ) :	(	" << track1Px << "	,	" << track1Py << "	,	" << track1Pz << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track1 Reference Point (x,y,z) :		(	" << xReference1 << "	,	" << yReference1 << "	,	" << zReference1 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "-----------------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "	TRACK2:" << std::endl;
	streamlog_out(DEBUG1) << *track2 << std::endl;
	double track2D0		= track2->getD0();
	double track2Z0		= track2->getZ0();
	double track2Phi	= track2->getPhi();
	double track2Omega	= track2->getOmega();
	double track2TanLambda	= track2->getTanLambda();
	double track2charge	= ( track2->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference2	= track2->getReferencePoint()[ 0 ];
	double yReference2	= track2->getReferencePoint()[ 1 ];
	double zReference2	= track2->getReferencePoint()[ 2 ];
	double xCenter2		= xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi );
	double yCenter2		= yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi );
	double zCenter2		= zReference2 + track2Z0;

	while ( track2Phi < 0.0 )
	{
		track2Phi += 2.0 * 3.14159265359;
	}
	while ( track2Phi >= 2.0 * 3.14159265359 )
	{
		track2Phi -= 2.0 * 3.14159265359;
	}

	double track2Pt = m_Bfield * 3.0e-4 / std::fabs( track2->getOmega() );
	double track2Px = track2Pt * std::cos( track2->getPhi() ) ;
	double track2Py = track2Pt * std::sin( track2->getPhi() ) ;
	double track2Pz = track2Pt * track2->getTanLambda() ;
	double track2Xs = track2->getReferencePoint()[ 0 ] - track2->getD0() * std::sin( track2->getPhi() );
	double track2Ys = track2->getReferencePoint()[ 1 ] + track2->getD0() * std::cos( track2->getPhi() );
	double track2Zs = track2->getReferencePoint()[ 2 ] + track2->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Track2 :	" << track2charge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Track2 ( Xc , Yc , Zc ) :	(	" << xCenter2 << "	,	" << yCenter2 << "	,	" << zCenter2 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Track2 ( Xc , Yc , Zc ) :	(	" << track2Xs << "	,	" << track2Ys << "	,	" << track2Zs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Track2 ( Px , Py , Pz ) :	(	" << track2Px << "	,	" << track2Py << "	,	" << track2Pz << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track2 Reference Point (x,y,z) :		(	" << xReference1 << "	,	" << yReference1 << "	,	" << zReference1 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "-----------------------------------------------------------------------------" << std::endl;

	////////////////////////////////////////////////////////////////////////////////
	//	double X = helixFreeParameter;
	//
	//	double xTrack		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( X ) / trackOmega;
	//	double yTrack		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( X ) / trackOmega;
	//	double zTrack		= zReference + trackZ0 - ( X - trackPhi ) * trackTanLambda / trackOmega;
	////////////////////////////////////////////////////////////////////////////////

//	double trackLengthToScan2D = 1000.0;
//	double deltaPhi1 = trackLengthToScan2D * track1Omega;
//	double deltaPhi2 = trackLengthToScan2D * track2Omega;
	double xMin = track1Phi - 3.14159265359 / 100.0;
	double xMax = track1Phi + 3.14159265359 / 100.0;
	double yMin = track2Phi - 3.14159265359 / 100.0;
	double yMax = track2Phi + 3.14159265359 / 100.0;


	DDMarlinCED::drawHelix( m_Bfield , track1charge , track1Xs , track1Ys , track1Zs , track1Px , track1Py , track1Pz , 2 , 1 , 0xff0000 , 0.0 , 1500.0 , 2000.0 , 0 );
	DDMarlinCED::drawHelix( m_Bfield , track2charge , track2Xs , track2Ys , track2Zs , track2Px , track2Py , track2Pz , 2 , 1 , 0x0061ff , 0.0 , 1500.0 , 2000.0 , 0 );

	double dphi1 = ( xMax - xMin ) / 100.0;
	double dphi2 = ( yMax - yMin ) / 100.0;
	for ( double phi1 = xMin ; phi1 <= xMax ; phi1 += dphi1 )
	{
		double track1x = xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) - std::sin( phi1 ) / track1Omega;
		double track1y = yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) + std::cos( phi1 ) / track1Omega;
		double track1z = zReference1 + track1Z0 - ( phi1 - track1Phi ) * track1TanLambda / track1Omega;
		ced_hit( track1x , track1y , track1z , 2 , 5 , 0x000000 );
	}
	for ( double phi2 = yMin ; phi2 <= yMax ; phi2 += dphi2 )
	{
		double track2x = xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) - std::sin( phi2 ) / track2Omega;
		double track2y = yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) + std::cos( phi2 ) / track2Omega;
		double track2z = zReference2 + track2Z0 - ( phi2 - track2Phi ) * track2TanLambda / track2Omega;
		ced_hit( track2x , track2y , track2z , 2 , 5 , 0x000000 );

	}

//	dist = std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - ( [ 7 ] + [ 8 ] * std::sin( x ) ) , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - ( [ 9 ] + [ 10 ] * std::cos( x ) ) , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * x ) - ( [ 11 ] + [ 12 ] * y ) , 2 ) )

	TF2 *distance = new TF2( "distance" , "std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - ( [ 6 ] + [ 7 ] * std::sin( y ) ) , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - ( [ 8 ] + [ 9 ] * std::cos( y ) ) , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * x ) - ( [ 10 ] + [ 11 ] * y ) , 2 ) )" , xMin , xMax , yMin , yMax );

	distance->SetParameter( 0 , xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) << std::endl;
	distance->SetParameter( 1 , -1.0 / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << -1.0 / track1Omega << std::endl;
	distance->SetParameter( 2 , yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) << std::endl;
	distance->SetParameter( 3 , 1.0 / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << 1.0 / track1Omega << std::endl;
	distance->SetParameter( 4 , zReference1 + track1Z0 + track1Phi * track1TanLambda / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << zReference1 + track1Z0 + track1Phi * track1TanLambda / track1Omega;
	distance->SetParameter( 5 , -1.0 * track1TanLambda / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] =  " << -1.0 * track1TanLambda / ( track1Omega ) << std::endl;
	distance->SetParameter( 6 , xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 6 ] =  " << xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) << std::endl;
	distance->SetParameter( 7 , -1.0 / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << -1.0 / track2Omega << std::endl;
	distance->SetParameter( 8 , yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << yReference1 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) << std::endl;
	distance->SetParameter( 9 , 1.0 / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 9 ] =  " << 1.0 / track2Omega << std::endl;
	distance->SetParameter( 10 , zReference2 + track2Z0 + track2Phi * track2TanLambda / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 10 ] =  " << zReference2 + track2Z0 + track2Phi * track2TanLambda / track2Omega;
	distance->SetParameter( 11 , -1.0 * track2TanLambda / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 11 ] =  " << -1.0 * track2TanLambda / ( track2Omega ) << std::endl;

	distance->SetRange( xMin , yMin , xMax , yMax );
//	distance->SetNpx( 1000 );
//	distance->SetNpy( 1000 );
	double minPhi1 = track1Phi;
	double minPhi2 = track2Phi;
	double minDistance = distance->GetMinimumXY( minPhi1 , minPhi2 );

	streamlog_out(DEBUG1) << "	Phi at Track 1 = " << minPhi1 << std::endl;
	streamlog_out(DEBUG1) << "	Phi at Track 2 = " << minPhi2 << std::endl;

	double xTrack1PCA	= xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) - std::sin( minPhi1 ) / track1Omega;
	double yTrack1PCA	= yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) + std::cos( minPhi1 ) / track1Omega;
	double zTrack1PCA	= zReference1 + track1Z0 - ( minPhi1 - track1Phi ) * track1TanLambda / track1Omega;
	streamlog_out(DEBUG1) << "	Track 1 PCA (x,y,z) = 	( " << xTrack1PCA << "	,	" << yTrack1PCA << "	,	" << zTrack1PCA << "	)" << std::endl;
	ced_hit( xTrack1PCA , yTrack1PCA , zTrack1PCA , 5 , 5 , 0xff0000 );
	PCAatTrack1 = TVector3( xTrack1PCA , yTrack1PCA , zTrack1PCA );

	double xTrack2PCA	= xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) - std::sin( minPhi2 ) / track2Omega;
	double yTrack2PCA	= yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) + std::cos( minPhi2 ) / track2Omega;
	double zTrack2PCA	= zReference2 + track2Z0 - ( minPhi2 - track2Phi ) * track2TanLambda / track2Omega;
	streamlog_out(DEBUG1) << "	Track 2 PCA (x,y,z) = 	( " << xTrack2PCA << "	,	" << yTrack2PCA << "	,	" << zTrack2PCA << "	)" << std::endl;
	ced_hit( xTrack2PCA , yTrack2PCA , zTrack2PCA , 5 , 5 , 0x0061ff );
	PCAatTrack2 = TVector3( xTrack2PCA , yTrack2PCA , zTrack2PCA );


	streamlog_out(DEBUG1) << "	Distance of Two Tracks = " << minDistance << std::endl;
	return minDistance;

}

/*

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
	double trackPt = m_Bfield * 3.0e-4 / std::fabs( track->getOmega() );
	double trackPx = trackPt * std::cos( track->getPhi() ) ;
	double trackPy = trackPt * std::sin( track->getPhi() ) ;
	double trackPz = trackPt * track->getTanLambda() ;
	double trackXs = track->getReferencePoint()[ 0 ] - track->getD0() * std::sin( track->getPhi() );
	double trackYs = track->getReferencePoint()[ 1 ] + track->getD0() * std::cos( track->getPhi() );
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


*/







/*

























void cheatTrueFlightDirection( EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , bool m_displayEvent )
{
	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	trueFlightDirection = TVector3( MotherHadron->getMomentum()[ 0 ] , MotherHadron->getMomentum()[ 1 ] , MotherHadron->getMomentum()[ 2 ] );
	trueFlightDirection.SetMag(1.0);
	if ( m_displayEvent ) ced_hit( MotherHadron->getVertex()[ 0 ] , MotherHadron->getVertex()[ 1 ] , MotherHadron->getVertex()[ 2 ] , 2 , 2 , 0x000000 );
}

int getParentHadronFlightDirection( 	EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
					TVector3 &trueFlightDirection , TVector3 &recoFlightDirection ,
					std::string inputPrimaryVertex , std::string inputBuildUpVertex ,
					std::string inputJetCollection , int vertexingScenario ,
					std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection ,
					float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent )
{
	int flightDirectionStatus = 0;
	std::vector<double> startVertexPosition;
	std::vector<double> endVertexPosition;
	std::vector<double> secondayVertex;
	std::vector<double> tritaryVertex;

	if ( m_displayEvent )
	{
		ced_hit( 0.0 , 0.0 , 0.100 , 1 , 1 , 0xc765c4 );
		ced_hit( 0.0 , 0.0 , 1.000 , 1 , 1 , 0xc765c4 );
		ced_hit( 0.0 , 0.0 , 10.00 , 1 , 1 , 0xc765c4 );
	}

	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	trueFlightDirection = TVector3( MotherHadron->getMomentum()[ 0 ] , MotherHadron->getMomentum()[ 1 ] , MotherHadron->getMomentum()[ 2 ] );
	streamlog_out(DEBUG1) << "		Test 2, |flightDirection| = " << trueFlightDirection.Mag() << std::endl;
	trueFlightDirection.SetMag(1.0);

	ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
	LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
	Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );

	std::vector<EVENT::Vertex *> buildUpVertices;
	int startVertexStatus = getStartVertexPosition( primaryVtx , startVertexPosition , m_displayEvent );
	int endVertexStatus = getEndVertexPosition( primaryVtx , buildUpVertices , linkedRecoLepton , endVertexPosition , helicesDistance , m_displayEvent );


	int secondayVertexStatus = getSecondaryVertex( pLCEvent , SLDLepton , inputPrimaryVertex , inputBuildUpVertex , false , secondayVertex , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , helicesDistance , SecondaryVertexPar , m_displayEvent );
	SecondaryVertexPar.push_back( MotherHadron->getEndpoint()[ 0 ] );
	SecondaryVertexPar.push_back( MotherHadron->getEndpoint()[ 1 ] );
	SecondaryVertexPar.push_back( MotherHadron->getEndpoint()[ 2 ] );
	streamlog_out(DEBUG1) << "	Parent Hadron Charge: " << MotherHadron->getCharge() << std::endl;
	streamlog_out(DEBUG1) << "	Primary Vertex status: " << startVertexStatus << std::endl;
	streamlog_out(DEBUG1) << "	Secondary Vertex status: " << secondayVertexStatus << std::endl;
	streamlog_out(DEBUG1) << "		Secondary Vertex" << std::endl;
	streamlog_out(DEBUG1) << "			True:(	" << SLDLepton->getVertex()[ 0 ] << "	, " << SLDLepton->getVertex()[ 1 ] << "	, " << SLDLepton->getVertex()[ 2 ] << " 		)" << std::endl;

	if ( ( secondayVertexStatus == 3 || secondayVertexStatus == 4 || secondayVertexStatus == 5 ) && startVertexStatus == 2 )
	{
		flightDirectionStatus = secondayVertexStatus;
		streamlog_out(DEBUG1) << "			Reco:(	" << secondayVertex[ 0 ] << "	, " << secondayVertex[ 1 ] << "	, " << secondayVertex[ 2 ] << " 		)" << std::endl;
		recoFlightDirection = TVector3( secondayVertex[ 0 ] - startVertexPosition[ 0 ] , secondayVertex[ 1 ] - startVertexPosition[ 1 ] , secondayVertex[ 2 ] - startVertexPosition[ 2 ] );
		streamlog_out(DEBUG1) << "		Test 3, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
		recoFlightDirection.SetMag( 1.0 );
		streamlog_out(DEBUG1) << "			Flight Direction Scenario: finding primary/secondary vertices and reconstruct flight direction" << std::endl;
	}
	else if ( secondayVertexStatus == 2 || secondayVertexStatus == 10 )
	{
		if ( vertexingScenario == 1 )
		{
//			ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
			recoFlightDirection = TVector3( linkedRecoLepton->getMomentum()[ 0 ] , linkedRecoLepton->getMomentum()[ 1 ] , linkedRecoLepton->getMomentum()[ 2 ] );
			streamlog_out(DEBUG1) << "		Test 4, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
			recoFlightDirection.SetMag( 1.0 );
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: finding primary/secondary vertices and reconstruct flight direction" << std::endl;
			if ( secondayVertexStatus == 2 ) flightDirectionStatus = 10;
			if ( secondayVertexStatus == 10 ) flightDirectionStatus = 11;
		}
		else if ( vertexingScenario == 2 )
		{
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: Assigning jet axis to the flight direction of parent hadron" << std::endl;
			int jetAssigningStatus = getJetAxis( pLCEvent , SLDLepton , recoFlightDirection , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
			streamlog_out(DEBUG1) << "		Test 5, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
			recoFlightDirection.SetMag( 1.0 );
			if ( jetAssigningStatus == 1 )
			{
				streamlog_out(DEBUG1) << "			Successfully assigned jet axis to the flight direction of parent hadron" << std::endl;
				if ( secondayVertexStatus == 2 ) flightDirectionStatus = 6;
				if ( secondayVertexStatus == 10 ) flightDirectionStatus = 11;
			}
			else
			{
				streamlog_out(DEBUG1) << "			!!! recoLepton does not belong to any jet !!!" << std::endl;
				flightDirectionStatus = 2;
			}
		}
		else if ( vertexingScenario == 3 )
		{
			streamlog_out(DEBUG1) << "			Flight Direction Scenario: Assigning flight direction of leading particle in the jet to the flight direction of parent hadron" << std::endl;
			int leadingParticleStatus = getLeadingParticleFlightDirection( pLCEvent , SLDLepton , recoFlightDirection , inputJetCollection , recoMCTruthLinkCollection , mcTruthRecoLinkCollection );
			streamlog_out(DEBUG1) << "		Test 6, |flightDirection| = " << recoFlightDirection.Mag() << std::endl;
			recoFlightDirection.SetMag( 1.0 );
			if ( leadingParticleStatus == 7 ) streamlog_out(DEBUG1) << "			Leading particle in the jet is a charged particle " << std::endl;
			if ( leadingParticleStatus == 8 ) streamlog_out(DEBUG1) << "			Leading particle in the jet is a neutral hadron " << std::endl;
			if ( leadingParticleStatus == 9 ) streamlog_out(DEBUG1) << "			Leading particle in the jet is a photon " << std::endl;
			if ( leadingParticleStatus == 2 ) streamlog_out(DEBUG1) << "			!!! recoLepton does not belong to any jet !!!" << std::endl;
			if ( secondayVertexStatus == 2 ) flightDirectionStatus = leadingParticleStatus;
			if ( secondayVertexStatus == 10 ) flightDirectionStatus = 11;
		}

	}
	else if ( secondayVertexStatus == 1 )
	{
		flightDirectionStatus = 1;
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
	streamlog_out(DEBUG1) << "			Flight Direction Status = " << flightDirectionStatus << std::endl;
	streamlog_out(DEBUG1) << "			Flight Direction Error: " << std::endl;
	streamlog_out(DEBUG1) << "					Alpha = " << acos( trueFlightDirection.Dot( recoFlightDirection ) ) * 180.0 / 3.14159265 << " deg" << std::endl;
	streamlog_out(DEBUG1) << "				   Cos(Alpha) = " << trueFlightDirection.Dot( recoFlightDirection ) << std::endl;
	return flightDirectionStatus;
}

int getStartVertexPosition( 	EVENT::Vertex *startVertex , std::vector<double> &startVertexPosition , bool m_displayEvent )
{
	int startVertexPositionStatus = 0;
	if( startVertex != NULL )
	{
		streamlog_out(DEBUG0) << "	There is a start vertex for parent hadron" << std::endl;
		startVertexPosition.push_back( startVertex->getPosition()[ 0 ] );
		startVertexPosition.push_back( startVertex->getPosition()[ 1 ] );
		startVertexPosition.push_back( startVertex->getPosition()[ 2 ] );
		startVertexPositionStatus = 1;
		streamlog_out(DEBUG0) << "		reco primary Vertex (x,y,z): 	" << startVertexPosition[ 0 ] << "	, " << startVertexPosition[ 1 ] << "	, " << startVertexPosition[ 2 ] << std::endl;
		if ( m_displayEvent ) ced_hit( startVertexPosition[ 0 ] , startVertexPosition[ 1 ] , startVertexPosition[ 2 ] , 2 , 2 , 0xff7300 );
	}
	else
	{
		streamlog_out(DEBUG0) << "	There is no start vertex for parent hadron" << std::endl;
		startVertexPositionStatus = -1;
	}
	return startVertexPositionStatus;
}

int getEndVertexPosition(	EVENT::Vertex *primaryVtx , std::vector<EVENT::Vertex *> buildUpVertices ,
				EVENT::ReconstructedParticle *linkedRecoLepton , std::vector<double> &endVertexPosition ,
				double helicesDistance , bool m_displayEvent )
{
	int secondayVertexStatus = -999;
	bool foundSecondaryVertex = false;

	return secondayVertexStatus;
}

int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton ,
			std::string inputPrimaryVertex , std::string inputBuildUpVertex ,
			bool cheatVertices , std::vector<double> &secondayVertex , std::string inputJetCollection ,
			std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection ,
			float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent )
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
		bool leptonInPrimaryVertex = false;
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
			ReconstructedParticle* primaryParticle = primaryVtx->getAssociatedParticle();
			int nPrimPar = ( primaryParticle->getParticles() ).size();
			for ( int i_par = 0 ; i_par < nPrimPar ; ++i_par )
			{
				ReconstructedParticle* particle = primaryParticle->getParticles()[ i_par ];
				if ( linkedRecoLepton == particle )
				{
					leptonInPrimaryVertex = true;
					secondayVertexStatus = 10;
				}
			}
			if ( !leptonInPrimaryVertex )
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
							secondayVertexStatus = 3;
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
						secondayVertexStatus = 2;
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
						streamlog_out(DEBUG1) << *linkedRecoLepton << std::endl;
						Track *leptonTrack	= linkedRecoLepton->getTracks()[ 0 ];

						int trueDSVertexStatus = getTrueDownStreamVertex( SLDLepton , trueDSVertex , trueDSVertexMomentum );
						streamlog_out(DEBUG1) << "	Found " << trueDSVertexStatus << " particles at True Down Stream vertex" << std::endl;
						streamlog_out(DEBUG1) << "	True Down Stream vertex at (x,y,z) = 	( " << trueDSVertex.X() << "	,	" << trueDSVertex.Y() << "	,	" << trueDSVertex.Z() << "	)" << std::endl;
						streamlog_out(DEBUG1) << "	True Down Stream vertex Momentum (Px,Py,Pz) = 	( " << trueDSVertexMomentum.Px() << "	,	" << trueDSVertexMomentum.Py() << "	,	" << trueDSVertexMomentum.Pz() << "	)" << std::endl;
						for ( unsigned int i_dsVTX = 0 ; i_dsVTX < DownStreamVertices.size() ; ++i_dsVTX )
						{
							Vertex* downStreamVertex = DownStreamVertices[ i_dsVTX ];
							std::vector<double> PCAatLeptonTrack;
							std::vector<double> PCAatDownStreamLine;
							streamlog_out(DEBUG1) << "	Checking Down Stream vertex at (x,y,z) = 	( " << downStreamVertex->getPosition()[ 0 ] << "	,	" << downStreamVertex->getPosition()[ 1 ] << "	,	" << downStreamVertex->getPosition()[ 2 ] << "	)" << std::endl;
							std::vector<double> dsVertex;
							dsVertex.push_back( downStreamVertex->getPosition()[ 0 ] );
							dsVertex.push_back( downStreamVertex->getPosition()[ 1 ] );
							dsVertex.push_back( downStreamVertex->getPosition()[ 2 ] );
							TVector3 dsMomentum = TVector3( downStreamVertex->getAssociatedParticle()->getMomentum()[ 0 ] , downStreamVertex->getAssociatedParticle()->getMomentum()[ 1 ] , downStreamVertex->getAssociatedParticle()->getMomentum()[ 2 ] );
	//						float dsDistance = std::sqrt( pow( dsVertex[ 0 ] - primaryVertex[ 0 ] , 2 ) + pow( dsVertex[ 1 ] - primaryVertex[ 1 ] , 2 ) + pow( dsVertex[ 2 ] - primaryVertex[ 2 ] , 2 ) );
							float dsDistance = intersectHelixLine( pLCEvent , leptonTrack , dsMomentum , dsVertex , PCAatLeptonTrack , PCAatDownStreamLine , inputPrimaryVertex );
							if ( dsDistance < minFlightDistance )
							{
								closetDownStreamVertex = downStreamVertex;
								minFlightDistance = dsDistance;
								foundDownStreamVertex = true;
							}
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
								// DDMarlinCED::drawHelix( bField , lepcharge , lepXs , lepYs , lepZs , lepPx , lepPy , lepPz , 1 , 2 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );

	//							DDMarlinCED::drawHelix( bField , lepcharge , PCAatLeptonTrack[ 0 ] , PCAatLeptonTrack[ 1 ] , PCAatLeptonTrack[ 2 ] , lepPx , lepPy , lepPz , 1 , 2 , 0xff0000 , 0.0 , 2100.0 , 3000.0 , 0 );


								double dsCharge = -1.0 * lepcharge;
								double scale = ( secondayVertexStatus == 4 ? -1.0 : -100.0 );
								double dsPx = scale * dsMomentum[ 0 ];
								double dsPy = scale * dsMomentum[ 1 ];
								double dsPz = scale * dsMomentum[ 2 ];
								double dsXs = dsVertex[ 0 ];
								double dsYs = dsVertex[ 1 ];
								double dsZs = dsVertex[ 2 ];

								ced_hit( PCAatDownStreamLine[ 0 ] , PCAatDownStreamLine[ 1 ] , PCAatDownStreamLine[ 2 ] , 2 , 1 , 0x0000ff );
								ced_hit( PCAatDownStreamLine[ 0 ] , PCAatDownStreamLine[ 1 ] , PCAatDownStreamLine[ 2 ] , 0 , 1 , 0x0000a8 );
								ced_hit( dsVertex[ 0 ] , dsVertex[ 1 ] , dsVertex[ 2 ] , 0 , 1 , 0x0000a8 );
								ced_line( dsVertex[ 0 ] , dsVertex[ 1 ] , dsVertex[ 2 ] , PCAatDownStreamLine[ 0 ] , PCAatDownStreamLine[ 1 ] , PCAatDownStreamLine[ 2 ] , 2 , 2 , 0x0000a8 );
//								double vMin = -1.0 * ( dsMomentum.Px() * ( dsVertex[ 0 ] - PCAatDownStreamLine[ 0 ] ) + dsMomentum.Py() * ( dsVertex[ 1 ] - PCAatDownStreamLine[ 1 ] ) + dsMomentum.Pz() * ( dsVertex[ 2 ] - PCAatDownStreamLine[ 2 ] ) ) / dsMomentum.Mag2();
//								double vMax = 0.0;
//								double dsTestX	= dsVertex[ 0 ] + dsMomentum.Px() * vMin;
//								double dsTestY	= dsVertex[ 1 ] + dsMomentum.Py() * vMin;
//								double dsTestZ	= dsVertex[ 2 ] + dsMomentum.Pz() * vMin;
//								ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
//								dsTestX	= dsVertex[ 0 ] + dsMomentum.Px() * vMax;
//								dsTestY	= dsVertex[ 1 ] + dsMomentum.Py() * vMax;
//								dsTestZ	= dsVertex[ 2 ] + dsMomentum.Pz() * vMax;
//								ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
//								for ( double v = vMin ; v <= vMax ; v += ( vMax - vMin ) / 1000.0 )
//								{
//									dsTestX	= dsVertex[ 0 ] + dsMomentum.Px() * v;
//									dsTestY	= dsVertex[ 1 ] + dsMomentum.Py() * v;
//									dsTestZ	= dsVertex[ 2 ] + dsMomentum.Pz() * v;
//									ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
//								}
								DDMarlinCED::drawHelix( bField , dsCharge , dsXs , dsYs , dsZs , dsPx , dsPy , dsPz , 1 , 2 , 0x0000ff , 0.0 , 2100.0 , 3000.0 , 0 );
								for ( unsigned int i_trk = 0 ; i_trk < ( downStreamVertex->getAssociatedParticle() )->getParticles().size() ; ++i_trk )
								{
									ReconstructedParticle* vtxRP = ( downStreamVertex->getAssociatedParticle() )->getParticles()[ i_trk ];
									Track *vtxTrack = vtxRP->getTracks()[ 0 ];
									double charge = ( vtxTrack->getOmega() > 0.0 ?  1.0 : -1.0 );
									double pT = bField * 3.0e-4 / std::fabs( vtxTrack->getOmega() );
									double Px = pT * std::cos( vtxTrack->getPhi() ) ;
									double Py = pT * std::sin( vtxTrack->getPhi() ) ;
									double Pz = pT * vtxTrack->getTanLambda() ;
									double Xs = dsVertex[ 0 ];
									double Ys = dsVertex[ 1 ];
									double Zs = dsVertex[ 2 ];
									DDMarlinCED::drawHelix( bField , charge , Xs , Ys , Zs , Px , Py , Pz , 1 , 2 , 0xb000e6 , 0.0 , 2100.0 , 3000.0 , 0 );
								}
							}
						}
						if ( foundDownStreamVertex )
						{
							streamlog_out(DEBUG1) << "	There are " << ( closetDownStreamVertex->getAssociatedParticle() )->getParticles().size() << " particle in DownStream Vertex" << std::endl;
							std::vector<double> PCAatLeptonTrack;
							std::vector<double> PCAatDownStreamLine;
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
							double lepXs = leptonTrack->getReferencePoint()[ 0 ] - leptonTrack->getD0() * std::sin( leptonTrack->getPhi() );
							double lepYs = leptonTrack->getReferencePoint()[ 1 ] + leptonTrack->getD0() * std::cos( leptonTrack->getPhi() );
							double lepZs = leptonTrack->getReferencePoint()[ 2 ] + leptonTrack->getZ0();
							double lepRs = std::sqrt( std::pow( lepXs - primaryVertex[ 0 ] , 2 ) + std::pow( lepYs - primaryVertex[ 1 ] , 2 ) + std::pow( lepZs - primaryVertex[ 2 ] , 2 ) );
							double dsRs  = std::sqrt( std::pow( downStreamPosition[ 0 ] - primaryVertex[ 0 ] , 2 ) + std::pow( downStreamPosition[ 1 ] - primaryVertex[ 1 ] , 2 ) + std::pow( downStreamPosition[ 2 ] - primaryVertex[ 2 ] , 2 ) );
							if ( lepRs <= dsRs )
							{
								secondayVertex = PCAatLeptonTrack;
							}
							else
							{
								secondayVertex = PCAatDownStreamLine;
							}
//							if ( m_displayEvent )
//							{
//								double bField = 3.5;//marlinutil::getBzAtOrigin();
//								double lepcharge = ( leptonTrack->getOmega() > 0.0 ?  1.0 : -1.0 );
//								double leppT = bField * 3.0e-4 / std::fabs( leptonTrack->getOmega() );
//								double lepPx = leppT * std::cos( leptonTrack->getPhi() ) ;
//								double lepPy = leppT * std::sin( leptonTrack->getPhi() ) ;
//								double lepPz = leppT * leptonTrack->getTanLambda() ;
//								lepXs = leptonTrack->getReferencePoint()[ 0 ] - leptonTrack->getD0() * std::sin( leptonTrack->getPhi() );
//								lepYs = leptonTrack->getReferencePoint()[ 1 ] + leptonTrack->getD0() * std::cos( leptonTrack->getPhi() );
//								lepZs = leptonTrack->getReferencePoint()[ 2 ] + leptonTrack->getZ0();
//								DDMarlinCED::drawHelix( bField , lepcharge , lepXs , lepYs , lepZs , lepPx , lepPy , lepPz , 1 , 2 , 0xff0000 , 0.0 , 2100.0 , 3000.0 , 0 );
//								DDMarlinCED::drawHelix( -bField , -lepcharge , lepXs , lepYs , lepZs , lepPx , lepPy , lepPz , 1 , 2 , 0xff0000 , 0.0 , 2100.0 , 3000.0 , 0 );
//								DDMarlinCED::drawHelix( bField , lepcharge , 0 , 0 , 0 , lepPx , lepPy , lepPz , 1 , 2 , 0x00ff40 , 0.0 , 2100.0 , 3000.0 , 0 );
//								DDMarlinCED::drawHelix( -bField , -lepcharge , 0 , 0 , 0 , lepPx , lepPy , lepPz , 1 , 2 , 0x00ff40 , 0.0 , 2100.0 , 3000.0 , 0 );
//								DDMarlinCED::drawHelix( bField , lepcharge , -3.34e+00 , 1.65e-01 , 5.36e+00 , -1.706949 , 0.402711 , 4.861839 , 1 , 2 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );
//								DDMarlinCED::drawHelix( -bField , -lepcharge , -3.34e+00 , 1.65e-01 , 5.36e+00 , -1.706949 , 0.402711 , 4.861839 , 1 , 2 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );

//								ced_hit( PCAatLeptonTrack[ 0 ] , PCAatLeptonTrack[ 1 ] , PCAatLeptonTrack[ 2 ] , 2 , 1 , 0xff0000 );
//								ced_hit( PCAatLeptonTrack[ 0 ] , PCAatLeptonTrack[ 1 ] , PCAatLeptonTrack[ 2 ] , 0 , 1 , 0xff0000 );
	//							DDMarlinCED::drawHelix( bField , lepcharge , PCAatLeptonTrack[ 0 ] , PCAatLeptonTrack[ 1 ] , PCAatLeptonTrack[ 2 ] , lepPx , lepPy , lepPz , 1 , 2 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );


//								double dsCharge = -1.0 * lepcharge;
//								double scale = ( secondayVertexStatus == 4 ? -1.0 : -100.0 );
//								double dsPx = scale * dsMomentum[ 0 ];
//								double dsPy = scale * dsMomentum[ 1 ];
//								double dsPz = scale * dsMomentum[ 2 ];
//								double dsXs = downStreamPosition[ 0 ];
//								double dsYs = downStreamPosition[ 1 ];
//								double dsZs = downStreamPosition[ 2 ];

//								ced_hit( PCAatDownStreamLine[ 0 ] , PCAatDownStreamLine[ 1 ] , PCAatDownStreamLine[ 2 ] , 0 , 1 , 0x0000a8 );
//								ced_hit( downStreamPosition[ 0 ] , downStreamPosition[ 1 ] , downStreamPosition[ 2 ] , 0 , 1 , 0x0000a8 );
//								double vMin = -1.0 * ( dsMomentum.Px() * ( downStreamPosition[ 0 ] - PCAatDownStreamLine[ 0 ] ) + dsMomentum.Py() * ( downStreamPosition[ 1 ] - PCAatDownStreamLine[ 1 ] ) + dsMomentum.Pz() * ( downStreamPosition[ 2 ] - PCAatDownStreamLine[ 2 ] ) ) / dsMomentum.Mag2();
//								double vMax = 0.0;
//								double dsTestX	= downStreamPosition[ 0 ] + dsMomentum.Px() * vMin;
//								double dsTestY	= downStreamPosition[ 1 ] + dsMomentum.Py() * vMin;
//								double dsTestZ	= downStreamPosition[ 2 ] + dsMomentum.Pz() * vMin;
//								ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
//								dsTestX	= downStreamPosition[ 0 ] + dsMomentum.Px() * vMax;
//								dsTestY	= downStreamPosition[ 1 ] + dsMomentum.Py() * vMax;
//								dsTestZ	= downStreamPosition[ 2 ] + dsMomentum.Pz() * vMax;
//								ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
//								for ( double v = vMin ; v <= vMax ; v += ( vMax - vMin ) / 1000.0 )
//								{
//									dsTestX	= downStreamPosition[ 0 ] + dsMomentum.Px() * v;
//									dsTestY	= downStreamPosition[ 1 ] + dsMomentum.Py() * v;
//									dsTestZ	= downStreamPosition[ 2 ] + dsMomentum.Pz() * v;
//									ced_hit( dsTestX , dsTestY , dsTestZ , 0 , 1 , 0x0000a8 );
//								}

//								ced_hit( PCAatDownStreamLine[ 0 ] , PCAatDownStreamLine[ 1 ] , PCAatDownStreamLine[ 2 ] , 2 , 1 , 0x0000ff );
	//							ced_hit( downStreamPosition[ 0 ] , downStreamPosition[ 1 ] , downStreamPosition[ 2 ] , 2 , 1 , 0x0000ff );
	//							DDMarlinCED::drawHelix( bField , dsCharge , dsXs , dsYs , dsZs , dsPx , dsPy , dsPz , 1 , 2 , 0x0000ff , 0.0 , 2100.0 , 3000.0 , 0 );
//								for ( unsigned int i_trk = 0 ; i_trk < ( closetDownStreamVertex->getAssociatedParticle() )->getParticles().size() ; ++i_trk )
//								{
//									ReconstructedParticle* vtxRP = ( closetDownStreamVertex->getAssociatedParticle() )->getParticles()[ i_trk ];
//									Track *vtxTrack = vtxRP->getTracks()[ 0 ];
//									double charge = ( vtxTrack->getOmega() > 0.0 ?  1.0 : -1.0 );
//									double pT = bField * 3.0e-4 / std::fabs( vtxTrack->getOmega() );
//									double Px = pT * std::cos( vtxTrack->getPhi() ) ;
//									double Py = pT * std::sin( vtxTrack->getPhi() ) ;
//									double Pz = pT * vtxTrack->getTanLambda() ;
//									double Xs = downStreamPosition[ 0 ];
//									double Ys = downStreamPosition[ 1 ];
//									double Zs = downStreamPosition[ 2 ];
	//								double Xs = vtxTrack->getReferencePoint()[ 0 ] - vtxTrack->getD0() * std::sin( vtxTrack->getPhi() );
	//								double Ys = vtxTrack->getReferencePoint()[ 1 ] + vtxTrack->getD0() * std::cos( vtxTrack->getPhi() );
	//								double Zs = vtxTrack->getReferencePoint()[ 2 ] + vtxTrack->getZ0();
	//								DDMarlinCED::drawHelix( bField , charge , Xs , Ys , Zs , Px , Py , Pz , 1 , 2 , 0xb000e6 , 0.0 , 2100.0 , 3000.0 , 0 );
//								}
//							}
	//						secondayVertex = PCAatDownStreamLine;
						}
						else
						{
							secondayVertexStatus = 2;
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
	bool jetAssigned = false;
	ReconstructedParticle* leadingParticle = NULL;
	ReconstructedParticle* linkedJet = NULL;
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
				jetAssigned = true;
				linkedJet = jet;
			}
		}
	}
	if ( jetAssigned )
	{
		int nParticles = ( linkedJet->getParticles() ).size();
		float leadingEnergy = 0.0;
		for ( int i_par = 0 ; i_par < nParticles ; ++i_par )
		{
			ReconstructedParticle* par = linkedJet->getParticles()[ i_par ];
			if ( par->getEnergy() > leadingEnergy )
			{
				leadingParticle = par;
				leadingEnergy = par->getEnergy();
			}
		}
		leadingParticleFlightDirection = TVector3( leadingParticle->getMomentum()[ 0 ] , leadingParticle->getMomentum()[ 1 ] , leadingParticle->getMomentum()[ 2 ] );
		streamlog_out(DEBUG1) << "		Test 9, |flightDirection| = " << leadingParticleFlightDirection.Mag() << std::endl;
		leadingParticleFlightDirection.SetMag( 1.0 );
		if ( leadingParticle->getCharge() == 0 )
		{
			if ( leadingParticle->getType() == 22 )
			{
				leadingParticleStatus = 9;
				streamlog_out(DEBUG1) << "		Leading particle in the jet is a photon " << std::endl;
			}
			else
			{
				leadingParticleStatus = 8;
				streamlog_out(DEBUG1) << "		Leading particle in the jet is a neutral hadron " << std::endl;
			}
		}
		else
		{
			leadingParticleStatus = 7;
			streamlog_out(DEBUG1) << "		Leading particle in the jet is a charged particle " << std::endl;
		}
	}
	else
	{
		leadingParticleStatus = 2;
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


//	HelixClass leptonHelix;
//	leptonHelix.Initialize_Canonical( leptonTrack->getPhi() , leptonTrack->getD0() , leptonTrack->getZ0() , leptonTrack->getOmega() , leptonTrack->getTanLambda() , m_Bfield );
//	lepcharge = ( leptonHelix.getOmega() > 0.0 ?  1.0 : -1.0 );
//	leppT = m_Bfield * 3.0e-4 / std::fabs( leptonHelix.getOmega() );
//	lepPx = leppT * std::cos( leptonHelix.getPhi0() ) ;
//	lepPy = leppT * std::sin( leptonHelix.getPhi0() ) ;
//	lepPz = leppT * leptonHelix.getTanLambda() ;
//	lepXs = leptonHelix.getReferencePoint()[ 0 ] - leptonHelix.getD0() * std::sin( leptonHelix.getPhi0() );
//	lepYs = leptonHelix.getReferencePoint()[ 1 ] + leptonHelix.getD0() * std::cos( leptonHelix.getPhi0() );
//	lepZs = leptonHelix.getReferencePoint()[ 2 ] + leptonHelix.getZ0();


//	DDMarlinCED::drawHelix( m_Bfield , lepcharge , lepXs , lepYs , lepZs , lepPx , lepPy , lepPz , 2 , 1 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );



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

*/


void drawMCParticles( EVENT::MCParticle *MotherHadron )
{
	double m_Bfield = MarlinUtil::getBzAtOrigin();
	ced_line( MotherHadron->getEndpoint()[ 0 ] , MotherHadron->getEndpoint()[ 1 ] , MotherHadron->getEndpoint()[ 2 ] , MotherHadron->getVertex()[ 0 ] , MotherHadron->getVertex()[ 1 ] , MotherHadron->getVertex()[ 2 ] , 1 , 1 , 0x00E0FF ); //DRAW UNSTABLE PARTICLES IN CYAN
	streamlog_out(DEBUG0) << " An unstable decay product: " << std::endl;
	streamlog_out(DEBUG0) << *MotherHadron << std::endl;
	for ( unsigned int i_daughter = 0 ; i_daughter < MotherHadron->getDaughters().size() ; ++i_daughter )
	{
		MCParticle *duaughter = MotherHadron->getDaughters()[ i_daughter ];
		if ( duaughter->getGeneratorStatus() == 1 )
		{
			if ( std::fabs( duaughter->getPDG() ) == 12 || std::fabs( duaughter->getPDG() ) == 14 || std::fabs( duaughter->getPDG() ) == 16 ) //DRAW Neutrinos in GRAY
			{
				ced_line( duaughter->getEndpoint()[ 0 ] , duaughter->getEndpoint()[ 1 ] , duaughter->getEndpoint()[ 2 ] , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 2 , 2 , 0x949494 );
				streamlog_out(DEBUG0) << " True Neutrino: " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
			else if ( std::fabs( duaughter->getPDG() ) == 11 || std::fabs( duaughter->getPDG() ) == 13 || std::fabs( duaughter->getPDG() ) == 15 ) // DRAW Leptons in RED
			{
				DDMarlinCED::drawHelix( m_Bfield , +1.0 , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 100.0 * duaughter->getMomentum()[ 0 ] , 100.0 * duaughter->getMomentum()[ 1 ] , 100.0 * duaughter->getMomentum()[ 2 ] , 1 , 1 , 0xfe1100 , 0.0 , 2100.0 , 3000.0 , 0 );
				streamlog_out(DEBUG0) << " True Lepton: " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
			else if ( std::fabs( duaughter->getCharge() ) <=0.1 ) //DRAW Neutral Particles in GREEN
			{
				ced_line( duaughter->getEndpoint()[ 0 ] , duaughter->getEndpoint()[ 1 ] , duaughter->getEndpoint()[ 2 ] , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 1 , 1 , 0x01be4b );
				streamlog_out(DEBUG0) << " True Neutral Particle: " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
			else //DRAW Charged Particles in BLUE
			{
				DDMarlinCED::drawHelix( m_Bfield , duaughter->getCharge() , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , duaughter->getMomentum()[ 0 ] , duaughter->getMomentum()[ 1 ] , duaughter->getMomentum()[ 2 ] , 1 , 1 , 0x5a6ffa , 0.0 , 2100.0 , 3000.0 , 0 );
				streamlog_out(DEBUG0) << " True Charged Particle: " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
		}
		else if ( duaughter->getGeneratorStatus() == 2 )
		{
			drawMCParticles( duaughter );
		}
		else //DRAW Unstable Particles in BLACK
		{
			if ( std::fabs( duaughter->getCharge() ) <=0.1 )
			{
				ced_line( duaughter->getEndpoint()[ 0 ] , duaughter->getEndpoint()[ 1 ] , duaughter->getEndpoint()[ 2 ] , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 1 , 1 , 0xbcbcbc );
				streamlog_out(DEBUG0) << " Other Unstable Neutral Particles (genStatus != 1 && genStatus != 2 ): " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
			else
			{
				DDMarlinCED::drawHelix( m_Bfield , duaughter->getCharge() , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , duaughter->getMomentum()[ 0 ] , duaughter->getMomentum()[ 1 ] , duaughter->getMomentum()[ 2 ] , 1 , 1 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );
				streamlog_out(DEBUG0) << " Other Unstable Charged Particles (genStatus != 1 && genStatus != 2 ): " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
		}
	}
}

void drawReconstructedParticle( EVENT::ReconstructedParticle *reconstructedParticle , EVENT::Vertex *primaryVertex , int colorCharged , int colorNeutral )
{
	double m_Bfield = MarlinUtil::getBzAtOrigin();
	int nTracks = ( reconstructedParticle->getTracks() ).size();
	if ( nTracks > 0 )
	{
		for ( int i_trk = 0 ; i_trk < nTracks ; ++i_trk )
		{
			Track *track = reconstructedParticle->getTracks()[ i_trk ];
			double trackCharge = ( track->getOmega() > 0.0 ?  1.0 : -1.0 );
			double trackPt = m_Bfield * 3.0e-4 / std::fabs( track->getOmega() );
			double trackPx = trackPt * std::cos( track->getPhi() ) ;
			double trackPy = trackPt * std::sin( track->getPhi() ) ;
			double trackPz = trackPt * track->getTanLambda() ;
			double trackXs = track->getReferencePoint()[ 0 ] - track->getD0() * std::sin( track->getPhi() );
			double trackYs = track->getReferencePoint()[ 1 ] + track->getD0() * std::cos( track->getPhi() );
			double trackZs = track->getReferencePoint()[ 2 ] + track->getZ0();
			DDMarlinCED::drawHelix( m_Bfield , trackCharge , trackXs , trackYs , trackZs , trackPx , trackPy , trackPz , 2 , 2 , colorCharged , 0.0 , 1500.0 , 2000.0 , 0 );
		}
	}
	else
	{
		double startPointX = primaryVertex->getPosition()[ 0 ];
		double startPointY = primaryVertex->getPosition()[ 1 ];
		double startPointZ = primaryVertex->getPosition()[ 2 ];
		double endPointX = ( reconstructedParticle->getClusters()[ 0 ]->getPosition()[ 0 ] - primaryVertex->getPosition()[ 0 ] ) * 2.0 / 3.0;
		double endPointY = ( reconstructedParticle->getClusters()[ 0 ]->getPosition()[ 1 ] - primaryVertex->getPosition()[ 1 ] ) * 2.0 / 3.0;
		double endPointZ = ( reconstructedParticle->getClusters()[ 0 ]->getPosition()[ 2 ] - primaryVertex->getPosition()[ 2 ] ) * 2.0 / 3.0;
		ced_line( endPointX , endPointY , endPointZ , startPointX , startPointY , startPointZ , 2 , 2 , colorNeutral );
	}
}
