//======================================================
// Make the root tree for mono-pole analysis
// CMSSW_7_6_3: Sherif Elgammal 17/03/2016
//  + work with AOD format
// CMSSW_8_0_X: Phat Srimanobhas
//  + Adapt to 8_0_X
//  + Analysis framework
//  + GIT
//======================================================

#include "NtupleAnalyzer.h"

using namespace edm;
using namespace reco;
using namespace std;

class NtupleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit NtupleAnalyzer(const edm::ParameterSet&);
  ~NtupleAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  void GenParticleTree(const edm::Event& evt); 
  void TrackTree(const Event& event, const EventSetup& setup);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override; 

  // ----------member data ---------------------------  
  TFile*  rootFile_;
  TTree* mytree;
  std::string outputFile_; // output file

  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesColl_;  
  edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;
  
  //branch: GEN particles
  std::vector<int> iGen;
  std::vector<int> idGen;
  std::vector<int> statusGen; 
  std::vector<float> ptGen;
  std::vector<float> etaGen;
  std::vector<float> phiGen;
  std::vector<float> massGen;
  std::vector<int> chargeGen;
  std::vector<float> EnergyGen;
  std::vector<float> pxGen;
  std::vector<float> pyGen;
  std::vector<float> pzGen;  

  //branch: tracks
  std::vector<float> Track_P;
  std::vector<float> Track_Pt;
  std::vector<float> Track_Eta;
  std::vector<float> Track_D0;
  std::vector<float> Track_Z0;
  std::vector<float> Track_Phi0;
  std::vector<float> Track_DeDx;
  std::vector<float> Track_DeDxError;
  std::vector<int> Track_SatStrips;
  std::vector<int> Track_TotStrips;
  std::vector<int> Track_NTracks;
  std::vector<int> Track_HitSatStrips;
  std::vector<int> Track_HitTotStrips;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NtupleAnalyzer::NtupleAnalyzer(const edm::ParameterSet& iConfig):
  outputFile_(iConfig.getParameter<std::string>("outputFile"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   rootFile_   = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms  
   genParticlesColl_ = consumes<std::vector<GenParticle> >(iConfig.getParameter<edm::InputTag>("genInfo"));
   generalTracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTracks"));
}


NtupleAnalyzer::~NtupleAnalyzer()
{
  delete rootFile_;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
NtupleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  
  //#ifdef THIS_IS_AN_EVENT_EXAMPLE
  //Handle<ExampleData> pIn;
  //iEvent.getByLabel("example",pIn);
  //#endif
  
  //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //ESHandle<SetupData> pSetup;
  //iSetup.get<SetupRecord>().get(pSetup);
  //#endif
    
  GenParticleTree(iEvent); 
  TrackTree(iEvent,iSetup);
  mytree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
NtupleAnalyzer::beginJob()
{  
  // go to *OUR* rootfile and book histograms
  rootFile_->cd();
  // Declare histograms
  mytree  = new TTree("tree","tr");
  
  //branch: GEN particles
  mytree->Branch("iGen",&iGen);
  mytree->Branch("idGen",&idGen);
  mytree->Branch("statusGen",&statusGen);
  mytree->Branch("ptGen",&ptGen);
  mytree->Branch("etaGen",&etaGen);
  mytree->Branch("phiGen",&phiGen);
  mytree->Branch("massGen",&massGen);
  mytree->Branch("chargeGen",&chargeGen);
  mytree->Branch("EnergyGen",&EnergyGen);
  mytree->Branch("pxGen",&pxGen);
  mytree->Branch("pyGen",&pyGen);
  mytree->Branch("pzGen",&pzGen);

  //branch: tracks
  mytree->Branch("Track_NTracks", &Track_NTracks);
  mytree->Branch("Track_Pt", &Track_Pt);
  mytree->Branch("Track_P", &Track_P);
  mytree->Branch("Track_Eta", &Track_Eta);
  mytree->Branch("Track_Phi0", &Track_Phi0);
  mytree->Branch("Track_D0", &Track_D0);
  mytree->Branch("Track_Z0", &Track_Z0);
  mytree->Branch("Track_DeDx", &Track_DeDx);
  mytree->Branch("Track_DeDxError", &Track_DeDxError);
  mytree->Branch("Track_HitSatStrips", &Track_HitSatStrips);
  mytree->Branch("Track_HitTotStrips", &Track_HitTotStrips);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtupleAnalyzer::endJob() 
{
  rootFile_->cd();
  mytree->Write();
  rootFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleAnalyzer);

//=============================================================
//
//            Method for Genrated Particles Tree
//
//=============================================================
void NtupleAnalyzer::GenParticleTree(const edm::Event& evt){
  int NbGen = 0;
  edm::Handle<GenParticleCollection> genParticles;
  evt.getByToken(genParticlesColl_, genParticles);
  if (!(genParticles.isValid())) return;  
  iGen.clear();
  idGen.clear();
  statusGen.clear();
  ptGen.clear();
  etaGen.clear();
  phiGen.clear();
  massGen.clear();
  chargeGen.clear();
  EnergyGen.clear();
  pxGen.clear();
  pyGen.clear();
  pzGen.clear();
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    iGen.push_back(NbGen);
    math::XYZTLorentzVector iParticle;
    iParticle = p.p4();
    idGen.push_back(p.pdgId());
    statusGen.push_back(p.status());
    ptGen.push_back(p.pt());
    etaGen.push_back(p.eta());
    phiGen.push_back(p.phi());
    massGen.push_back(iParticle.M());
    chargeGen.push_back(p.charge());
    EnergyGen.push_back(p.energy());
    pxGen.push_back(p.px());
    pyGen.push_back(p.py());
    pzGen.push_back(p.pz());
    NbGen++;
  }
}

void NtupleAnalyzer::TrackTree(const edm::Event& event, const edm::EventSetup& setup){
  int NbTracks = 0;
  Track_NTracks.clear(); 
  Track_Pt.clear(); 
  Track_P.clear(); 
  Track_Eta.clear();  
  Track_D0.clear(); 
  Track_Z0.clear(); 
  Track_Phi0.clear(); 
  Track_DeDx.clear();
  Track_DeDxError.clear();
  Track_HitTotStrips.clear(); 
  Track_HitSatStrips.clear(); 
  Handle<reco::TrackCollection> _hTracks;
  //Handle<ValueMap<reco::DeDxData> > _hDeDx;
  //Handle<reco::DeDxHitInfo> _hDeDx;
  event.getByToken(generalTracksToken_, _hTracks);
  //event.getByToken(dedxHarmonic2_, _hDeDx);
  // Loop over the tracks
  int _NTracks = _hTracks->size();
  for(int i = 0; i!=_NTracks; i++){
    edm::Ref<std::vector<reco::Track> > TrackRef(_hTracks, i);
    NbTracks++;
    Track_NTracks.push_back(NbTracks);
    Track_Pt.push_back(TrackRef->pt());
    Track_P.push_back(TrackRef->p());
    Track_Eta.push_back(TrackRef->eta());
    Track_Phi0.push_back(TrackRef->phi());
    Track_D0.push_back(TrackRef->d0());
    Track_Z0.push_back(TrackRef->dsz());
    //Track_DeDx.push_back((*_hDeDx.product())[TrackRef].dEdx());
    //Track_DeDxError.push_back((*_hDeDx.product())[TrackRef].dEdxError());
    //Track_HitTotStrips.push_back((*_hDeDx.product())[TrackRef].numberOfMeasurements());
    //Track_HitSatStrips.push_back((*_hDeDx.product())[TrackRef].numberOfSaturatedMeasurements());
  }
}
