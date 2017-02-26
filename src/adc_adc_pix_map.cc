#include<iostream>
#include<stdint.h>
#include<boost/container/flat_map.hpp>
#include<TChain.h>
#include<TFile.h>
#include<TTree.h>
#include<THnSparse.h>
#include<TString.h>
#include<TRegexp.h>
#include"feio.h"

//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 if(argc<3){
	std::cerr<<"USAGE: "<<argv[0]<<" channel# root_filename[s]"<<std::endl;
	exit(111);
 }

 int nchannel = TString(argv[1]).Atoi();
 int nasic = TString(argv[2]).Contains("2ASIC") ? 2 : 3;
 int iasic = nchannel/64;
 if(nchannel<0 || nchannel>191 || TString(argv[1]).Contains(TRegexp("[a-zA-Z]")) || (nasic==2 && iasic==1)){
	std::cerr<<"you chose empty channel"<<std::endl;
	exit(222);
 }

 UShort_t fadc[RICHfrontend::NCHANNELS];

 TChain* h22 = new TChain("h22");
 for(int iarg=2;iarg<argc;iarg++)
	if(TString(argv[iarg]).Contains(TRegexp(".root$")))
		h22->AddFile(argv[iarg]);
 h22->SetBranchAddress("fadc", fadc);

 Int_t vbin[] = {300,300,64};
 Double_t vmin[] = {0.5,0.5,0.5};
 Double_t vmax[] = {3000.5,3000.5,64.5};
 THnSparseI *hh = new THnSparseI("hs","hs",3,vbin,vmin,vmax);

 boost::container::flat_map<uint32_t, uint32_t> hmap;
 long unsigned int nen = h22->GetEntries();
 for(int ien=0; ien<nen; ien++){
	h22->GetEntry(ien);

//	for(int ich=iasic*64;ich<(iasic+1)*64;ich++){
//		double vv[] = {fadc[nchannel], fadc[ich], RICHfrontend::chan2pix[ich%64]};
//		hh->Fill(vv);
//		uint32_t key = (RICHfrontend::chan2pix[ich%64]<<24) | (fadc[nchannel]<<12) | (fadc[ich]);
//		hmap[key]++;
//	}
 }

 delete h22;

/*
 TFile* ff = new TFile(Form("adcPixMap_%s", basename(argv[2])), "recreate");
 TTree* tt = new TTree("tmap", "adc pixel mapping");

 UInt_t key, val;
 tt->Branch("key", &key, "key/i");
 tt->Branch("val", &val, "val/i");

 for (boost::container::flat_map<uint32_t, uint32_t>::iterator it=hmap.begin(); it!=hmap.end(); ++it){
	key = it->first;
	val = it->second;
	tt->Fill();
 }

 ff->Write();

 delete tt, ff;
*/

 return 0;
}

