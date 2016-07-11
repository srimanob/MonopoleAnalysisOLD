//C++
#include <memory>
#include <vector>
#include <string>
#include <set>
#include <stdio.h>
#include <math.h>

//ROOT
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "TMath.h"

//CMSSW
//#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"



//Maybe use later
//#include "DataFormats/JetReco/interface/GenJetCollection.h"

//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/Candidate/interface/Candidate.h"

//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"

//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//#include "DataFormats/Math/interface/deltaR.h"

// Ecal includes
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/Records/interface/CaloTopologyRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
//#include "DataFormats/EcalDetId/interface/EBDetId.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
//#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/HLTReco/interface/TriggerObject.h"
//#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
//#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"    
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/MuonReco/interface/MuonCocktails.h"
//#include "DataFormats/Common/interface/RefCore.h"
//#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
//#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
//#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
//#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
//#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"




//#include "DataFormats/Common/interface/Handle.h"

//#include "DataFormats/Common/interface/ValueMap.h"


//#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
//#include "FWCore/ParameterSet/interface/FileInPath.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLHcal.h"
