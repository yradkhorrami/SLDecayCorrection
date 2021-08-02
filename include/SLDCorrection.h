#ifndef SLDCorrection_h
#define SLDCorrection_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include "UTIL/LCRelationNavigator.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;
class TF1;

using namespace lcio ;
using namespace marlin ;

class SLDCorrection : public Processor
{
public:
	virtual Processor *newProcessor()
	{
		return new SLDCorrection;
	}
	SLDCorrection();
	virtual ~SLDCorrection() = default;
	SLDCorrection( const SLDCorrection& ) = delete;
	SLDCorrection &operator = ( const SLDCorrection& ) = delete;
	virtual void init();
	virtual void Clear();
	virtual void processRunHeader();
	virtual void processEvent( EVENT::LCEvent *pLCEvent );
	virtual void check( EVENT::LCEvent *pLCEvent );
	virtual void end();

private:
	typedef std::vector<int>			IntVector;
	typedef std::vector<double>			DoubleVector;
	typedef std::vector<float>			FloatVector;

	std::string				m_mcParticleCollection{};
	std::string				m_inputPfoCollection{};
	std::string				m_inputJetCollection{};
	std::string				m_jetFlavour{};
	std::string				m_inputPrimaryVertex{};
	std::string				m_inputBuildUpVertex{};
	std::string				m_RecoMCTruthLinkCollection{};
	std::string				m_MCTruthRecoLinkCollection{};
	std::string				m_TrackMCTruthLinkCollection{};
	std::string				m_MCTruthTrackLinkCollection{};
	std::string				m_ClusterMCTruthLinkCollection{};
	std::string				m_MCTruthClusterLinkCollection{};
	std::string				m_SLDNuCollection{};
	std::string				m_rootFile{};



	TFile					*m_pTFile{};
	TTree					*m_pTTree{};

};
#endif
