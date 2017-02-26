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
};


goodRICHEvent::goodRICHEvent(TString _fname):fname(_fname)
{
	nasic = fname.Contains("2ASIC") ? 2 : 3;
}


void goodRICHEvent::Fill(rawEvent &rev){
	std::copy(rev.fadc, rev.fadc+NCHANNELS, fadc);

	for(int ichan=0; ichan<NCHANNELS; ichan++)
	{
		for(int ich=0;ich<NCHANNELS;ich++){
			uint32_t key = (ich<<24) | (fadc[ichan]<<12) | (fadc[ich]);
			hh[ichan][key]++;

			if(nasic==2 && ich==63) ich += 64;
		}
		if(nasic==2 && ichan==63) ichan += 64;
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

	for(chan=0;chan<NCHANNELS;chan++)
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

