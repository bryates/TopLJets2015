{
  //gROOT->ProcessLine(".L /afs/cern.ch/user/b/byates/TopAnalysis/src/CharmEvent.cc");
  gROOT->ProcessLine(".L splot_d0_mu_tag.C");
  gROOT->ProcessLine(".L convert.C");
  bool all=0;
  int epoch=2;
  TString sample("d0_mu_tag");
submit("172v5",0,"up",1);
submit("172v5",0,"725",1);
submit("isr-up",0,"uup",1);
submit("isr-down",0,"scentral",1);
submit("isr-down",0,"1025",1);
submit("isr-down",0,"675",1);
submit("fsr-up",0,"uup",1);
submit("fsr-up",0,"central",1);
submit("fsr-up",0,"ccentral",1);
submit("fsr-up",0,"cccentral",1);
submit("fsr-up",0,"scentral",1);
submit("fsr-up",0,"down",1);
submit("fsr-up",0,"dddown",1);
submit("fsr-down",0,"925",1);
submit("QCD_erdON",0,"675",1,true);
  if(all) splot("172.5", false, epoch);
  if(all) splot("isr-up", false, epoch);
  if(all) splot("isr-down", false, epoch);
  if(all) splot("fsr-up", false, epoch);
  if(all) splot("fsr-down", false, epoch);
  if(all) splot("ueup", false, epoch);
  if(all) splot("uedown", false, epoch);
  if(all) splot("erdON", false, epoch);
  if(all) splot("QCD_erdON", false, epoch);
  if(all) splot("GluonMove_erdON", false, epoch);
  std::vector<TString> syst = {"LEP", "TRIGGER", "TRK", "PU", "PI", "JER" };
  for(auto & it : syst) {
    if(all) splot("up_"+it, false, epoch);
    if(all) splot("down_"+it, false, epoch);
  }
  if(all) splot("hdampup", false, epoch);
  if(all) splot("hdampdown", false, epoch);
  if(all) splot("Data", true, epoch);
  /*
  splot_d0(name,false,"d0");
  splot_d0(name,false,"d0up");
  splot_d0(name,false,"d0down");
  */
}

TString current("");
bool bNew(false);
void submit(TString name="", bool isData=false, TString frag="", int ep=0, bool sub=false) {
  /*
  TString cmd = TString::Format("bsub -q 8nm /afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh root -q -b '/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/splot_d0.C(\"%s\",%d,\"%s\")'",name.Data(),isData,frag.Data());
  std::cout << cmd << std::endl;
  gSystem->Exec(cmd);
  */
  TString subfile(TString::Format("sub_d0_mu_tag%s",name.Data()));
  if(current != name) {
    subfile = TString::Format("sub_d0_mu_tag%s",current.Data());
    bNew = current == "" ? true : false;
    std::cout << current << "-> " << name << std::endl;
    if(!bNew) gSystem->Exec(TString::Format("condor_submit %s -batch-name %s-%s%d",subfile.Data(),current.Data(),sample.Data(),epoch));
    if(!bNew) gSystem->Exec(TString::Format("rm %s",subfile.Data()));
    current = name;
  }
  gSystem->Exec(TString::Format("echo \"universe              = vanilla\" >> %s",subfile.Data()));
  gSystem->Exec(TString::Format("echo \"executable            = cond_submit_d0_mu_tag.sh\" >> %s",subfile.Data()));
  gSystem->Exec(TString::Format("echo \"arguments             = \\$(ClusterID) \\$(ProcId) %s %d %d %s\" >> %s",name.Data(),isData,ep,frag.Data(),subfile.Data()));
  gSystem->Exec(TString::Format("echo \"output                = log_%s_%s_mu_tag_%d.out\" >> %s",name.Data(),frag.Data(),ep,subfile.Data()));
  gSystem->Exec(TString::Format("echo \"error                 = log_%s_%s_mu_tag_%d.err\" >> %s",name.Data(),frag.Data(),ep,subfile.Data()));
  gSystem->Exec(TString::Format("echo \"log                   = log_%s_%s_mu_tag_%d.log\" >> %s",name.Data(),frag.Data(),ep,subfile.Data()));
  gSystem->Exec(TString::Format("echo \"+JobFlavour           = longlunch\" >> %s",subfile.Data()));
  gSystem->Exec("echo \"request_disk=2G\" >> sub");
  gSystem->Exec("echo \"request_memory=2G\" >> sub");
  gSystem->Exec("echo \"request_cpus=4\" >> sub");
  gSystem->Exec(TString::Format("echo \"Should_Transfer_Files = NO\" >> %s",subfile.Data()));
  gSystem->Exec(TString::Format("echo \"queue\" >> %s",subfile.Data()));
  //if(bNew) gSystem->Exec(TString::Format("condor_submit sub_d0_mu -batch-name %s-%s%d",name.Data(),sample.Data(),ep));
  //if(bNew) gSystem->Exec("rm sub_d0_mu");
  /*
  if(sub) gSystem->Exec(TString::Format("condor_submit sub_d0_mu -batch-name %s_%s_mu_tag",name.Data(),frag.Data()));
  if(sub) gSystem->Exec("rm sub_d0_mu");
  */
  bNew = false;
  if(sub) gSystem->Exec(TString::Format("condor_submit %s -batch-name %s-%s%d",subfile.Data(),current.Data(),sample.Data(),epoch));
  if(sub) gSystem->Exec(TString::Format("rm %s",subfile.Data()));
}

void splot(TString name="fsr-up", bool isData=false, int epoch=0) {
  //splot_d0(name);
  if(!isData) {
  submit(name,false,"",epoch,false);
  submit(name,false,"up",epoch,false);
  submit(name,false,"uup",epoch,false);
  submit(name,false,"uuup",epoch,false);
  submit(name,false,"central",epoch,false);
  submit(name,false,"ccentral",epoch,false);
  submit(name,false,"cccentral",epoch,false);
  submit(name,false,"scentral",epoch,false);
  submit(name,false,"down",epoch,false);
  submit(name,false,"ddown",epoch,false);
  submit(name,false,"dddown",epoch,false);
  submit(name,false,"sdown",epoch,false);
  submit(name,false,"1025",epoch,false);
  submit(name,false,"1075",epoch,false);
  submit(name,false,"925",epoch,false);
  submit(name,false,"725",epoch,false);
  submit(name,false,"700",epoch,false);
  submit(name,false,"675",epoch,false);
  submit(name,false,"625",epoch,false);
  submit(name,false,"600",epoch,false);
  submit(name,false,"555",epoch,false);
  submit(name,false,"fit",epoch,true);
  }
  else submit("Data",true,"",epoch,true);
}
void frag() {
  frag_ratio();
}
