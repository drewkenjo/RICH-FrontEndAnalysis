#include<iostream>
#include<stdint.h>
#include<TChain.h>
#include<TFile.h>
#include<TTree.h>
#include<TString.h>
#include<TRegexp.h>
#include<boost/container/flat_map.hpp>
#include"feio.h"

int main(int argc, char** argv)
{
 if(argc<3 || TString(argv[1]).Contains(TRegexp("[a-zA-Z]"))){
	std::cerr<<"USAGE: "<<argv[0]<<" ch# root_filename[s]"<<std::endl;
	exit(111);
 }
 int channel = TString(argv[1]).Atoi();
 int nasic = TString(argv[2]).Contains("2ASIC") ? 2 : 3;
 if(channel<0 || channel>191 || (nasic==2 && channel/64==1)){
	std::cerr<<"you chose wrong channel"<<std::endl;
	exit(111);
 }

 TChain* tt = new TChain("tmap");
 for(int iarg=2;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);

 UInt_t key, val;
 UChar_t chan;
 tt->SetBranchAddress("chan", &chan);
 tt->SetBranchAddress("key", &key);
 tt->SetBranchAddress("val", &val);

 boost::container::flat_map<uint32_t, uint32_t> hmap;

 long unsigned int nen = tt->GetEntries();
 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);
	if(chan==channel)
		hmap[key] += val;
 }
 delete tt;

 TFile* fmap = new TFile(Form("full%03d_%s", channel, basename(argv[2])), "recreate");
 TTree* tmap = new TTree("tmap", Form("adc map for channel %03d, pixel %02d", channel, RICHfrontend::chan2pix[channel%64]));
 tmap->Branch("key", &key, "key/i");
 tmap->Branch("val", &val, "val/i");

 for (boost::container::flat_map<uint32_t, uint32_t>::iterator it=hmap.begin(); it!=hmap.end(); ++it){                                             
	key = it->first;
	val = it->second;
	tmap->Fill();
 }

 fmap->Write();
 delete tmap, fmap;

 return 0;
}

