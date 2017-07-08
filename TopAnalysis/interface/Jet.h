#ifndef _Jet_h_
#define _Jet_h_

#include <TLorentzVector.h>
#include <vector>

class pfTrack {

 public:
  pfTrack(TLorentzVector p4, float dxy, float dxyE, float dz, float dzE, int pfid, int quality, bool highPurity);
  pfTrack(TLorentzVector p4, int pfid);
  ~pfTrack();
  float Pt();
  float Eta();
  float Phi();
  float M();
  float DeltaR(pfTrack &rhs);
  int getPfid();
  inline int getPdgId() { return getPfid(); }
  int charge();
  int getQuality();
  bool highPurity();
  float getDxy();
  float getDxyE();
  float getDz();
  float getDzE();
  TLorentzVector &getVec();
  //inline TLorentzVector operator+(pfTrack &rhs) { return vec_+rhs.getVec() ; }
  void setPfid(int pfid);
  void setMass(float mass);

 private:
  TLorentzVector vec_;
  float dxy_;
  float dxyE_;
  float dz_;
  float dzE_;
  int pfid_;
  int quality_;
  bool highPurity_;

};

typedef std::pair<pfTrack,int> IdTrack;

class Jet {

 public:
  Jet(TLorentzVector p4, float csv, int idx);
  Jet(TLorentzVector p4, float csv, int idx, float chargedPt, float PFPt);
  Jet(TLorentzVector p4, float csv, int idx, float chargedPt, float PFPt, int genJet);
  ~Jet();
  //void addTrack(TLorentzVector p4, int pfid);
  //void addTrack(pfTrack pf, int pfid);
  void addTrack(pfTrack pf);
  void addTrack(int idx);
  void addDxy(float dxy, float dxyE);
  void addDz(float dz, float dzE);
  void addDz(int idx);
  TLorentzVector &getVec();
  float &getCSV();
  float getPt();
  float &getChargedPt();
  float &getPFPt();
  int &getGenJet();
  int &getJetIndex();
  int &getIndex(int idx);
  float &getDxy(int idx);
  float &getDz(int idx);
  float &getDxyE(int idx);
  float &getDzE(int idx);
  std::vector<pfTrack> &getTracks();
  pfTrack getTrack(int idx,float mass);
  void sortTracksByPt();
 
 private:
  TLorentzVector p4_;
  float csv_;
  int jetindex_;
  int idx_;
  float chargedPt_;
  float PFPt_;
  int genJet_;
  std::vector<pfTrack> trks_;
  std::vector<int> index_;
  std::vector<float> dxy_;
  std::vector<float> dxyE_;
  std::vector<float> dz_;
  std::vector<float> dzE_;

};


inline bool sortJetsByPt(Jet i, Jet j)  { return i.getVec().Pt() > j.getVec().Pt(); }
inline bool sortJetsByCSV(Jet i, Jet j) { return i.getCSV() > j.getCSV(); }
inline bool sortIdTracksByPt(IdTrack i, IdTrack j)  { return i.first.Pt() > j.first.Pt(); }
inline bool sortpfTracksByPt(pfTrack i, pfTrack j)  { return i.Pt() > j.Pt(); }

#endif
