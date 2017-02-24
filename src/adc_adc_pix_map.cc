#include<iostream>
#include<stdint.h>
#include<boost/container/flat_map.hpp>
#include<list>
#include<TFile.h>
#include<TTree.h>
#include<TString.h>
#include<TRegexp.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(TString);
	~goodRICHEvent();
	void Fill(rawEvent&);

  private:
	int nasic;
	TString fname;
	boost::container::flat_map<uint32_t, uint32_t> hh[NCHANNELS];
	uint8_t ipix[64];
};


goodRICHEvent::goodRICHEvent(TString _fname):fname(_fname)
{
	nasic = fname.Contains("2ASIC") ? 2 : 3;
	uint8_t ipixels[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};
	std::copy(ipixels, ipixels+64, ipix);
}


void goodRICHEvent::Fill(rawEvent &rev){
	std::copy(rev.fadc, rev.fadc+NCHANNELS, fadc);

	for(int ichan=0; ichan<NCHANNELS; ichan++)
	{
		int iasic = ichan/64;
		for(int ich=iasic*64;ich<(iasic+1)*64;ich++){
			uint32_t key = (ipix[ich%64]<<24) | (fadc[ichan]<<12) | (fadc[ich]);
			hh[ichan][key]++;
		}
		if(nasic==2 && ichan==63) ichan+=64;
	}
}


goodRICHEvent::~goodRICHEvent(){
	TFile* ff = new TFile(Form("adcPixMap_%s", basename(fname.Data())), "recreate");
	TTree* tt = new TTree("tmap", "adc pixel mapping");

	UInt_t key, val;
	UChar_t chan;
	tt->Branch("chan", &chan, "chan/b");
	tt->Branch("key", &key, "key/i");
	tt->Branch("val", &val, "val/i");

	for(int chan=0;chan<NCHANNELS;chan++)
	for (boost::container::flat_map<uint32_t, uint32_t>::iterator it=hh[chan].begin(); it!=hh[chan].end(); ++it){
		key = it->first;
		val = it->second;
		tt->Fill();
	}

	ff->Write();

	delete tt, ff;
}


//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev(argv[1]);

 for(int iarg=1;iarg<argc;iarg++){
	if(!TString(argv[iarg]).Contains(TRegexp(".root$"))) continue;
	if(iarg%1==0) std::cout<<iarg<<" / "<<argc-1<<std::endl;

	TFile* ff = new TFile(argv[iarg]);
	TTree* tt = (TTree*) ff->Get("h22");

	tt->SetBranchAddress("fadc", rawEv.fadc);

	int nen = tt->GetEntries();
	for(int ien=0; ien<nen; ien++){
		tt->GetEntry(ien);
		ev.Fill(rawEv);
	}

	delete tt, ff;
 }

 return 0;
}

