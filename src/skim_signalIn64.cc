#include<iostream>
#include<TChain.h>
#include<TFile.h>
#include<TTree.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TCanvas.h>
#include<TH1I.h>
#include<TF1.h>
#include"feio.h"


using namespace RICHfrontend;

/////////////////////////////////////////////////////
class skimmer{
  public:
	skimmer(TString, UShort_t*);
	~skimmer();
	void Fill(rawEvent&);

  private:
	TFile* ff[3];
	TTree* tt[3];

	UShort_t m4sig[NCHANNELS];

	int nasic;

	UShort_t fadc[NPIXELS];

	UShort_t fnedge;
	UShort_t ftdc[MAXEDGES];
	UShort_t fchan[MAXEDGES];
	Bool_t fpolar[MAXEDGES];
};


////////////////////////////////////////////////////
skimmer::skimmer(TString fname, UShort_t *_m4sig){
	nasic = fname.Contains("2ASIC") ? 2 : 3;
	std::copy(_m4sig, _m4sig+NCHANNELS, m4sig);

	for(int iasic=0;iasic<3;iasic++){
		if(iasic==1 && nasic==2) continue;

		ff[iasic] = new TFile(Form("skim64_%d_%s", iasic, basename(fname.Data())), "recreate");
		tt[iasic] = new TTree("h22", Form("skimmed data for 64 pixels in asic %d/%d", iasic+1, nasic));

		tt[iasic]->Branch("fadc", fadc, "fadc[64]/s");

		tt[iasic]->Branch("fnedge", &fnedge, "fnedge/s");
		tt[iasic]->Branch("fchan", fchan, "fchan[fnedge]/s");
		tt[iasic]->Branch("ftdc", ftdc, "ftdc[fnedge]/s");
		tt[iasic]->Branch("fpolar", fpolar, "fpolar[fnedge]/O");
	}
}

skimmer::~skimmer(){
	for(int iasic=0;iasic<3;iasic++){
		if(iasic==1 && nasic==2) continue;

		ff[iasic]->Write();
	}
}

void skimmer::Fill(rawEvent &rev)
{
	for(int iasic=0;iasic<3;iasic++){
		if(iasic==1 && nasic==2) continue;

		bool fillTree = false;
		for(int ichPMT=0; ichPMT<NPIXELS; ichPMT++){
			int ich = ichPMT+iasic*64;
			if(rev.fadc[ich] > m4sig[ich]){
				fillTree = true;
				break;
			}
		}

		if(fillTree){
			for(int ich=0;ich<NPIXELS;ich++)
				fadc[ich] = rev.fadc[iasic*64+ich];

			fnedge=0;
			for(int iedge=0; iedge<rev.fnedge; iedge++){
				if(rev.fchan[iedge]/64==iasic){
					fchan[fnedge] = rev.fchan[iedge];
					ftdc[fnedge] = rev.ftdc[iedge];
					fpolar[fnedge] = rev.fpolar[iedge];
					fnedge++;
				}
			}

			tt[iasic]->Fill();
		}
	}
}


/////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 TString fname(argv[1]);

 int nasic = fname.Contains("2ASIC") ? 2 : 3;
 rawEvent rawEv;

 TChain* h22 = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
 if(TString(argv[iarg]).Contains(TRegexp(".root$")))
	h22->AddFile(argv[iarg]);

 h22->SetBranchAddress("fadc", rawEv.fadc);

 h22->SetBranchAddress("fnedge", &rawEv.fnedge);
 h22->SetBranchAddress("ftdc", rawEv.ftdc);
 h22->SetBranchAddress("fchan", rawEv.fchan);
 h22->SetBranchAddress("fpolar", rawEv.fpolar);


////////////////////////////////////////////////////
 TH1I* h1[NCHANNELS];
 for(int ich=0;ich<NCHANNELS;ich++)
	h1[ich] = new TH1I(Form("h1_%03d",ich), Form("channel %03d",ich), 400,300.5,700.5);
 TF1* f1 = new TF1("f1","gaus",0,1000);
 UShort_t m4sig[NCHANNELS]={};


////////////////////////////////////////////////////
 unsigned long int nen = h22->GetEntries();
 for(long unsigned int ien=0;ien<std::min(1000000.0, (double) nen);ien++){
	h22->GetEntry(ien);
	for(int ich=0;ich<NCHANNELS;ich++){
		h1[ich]->Fill(rawEv.fadc[ich]);
		if(nasic==2 && ich==63) ich+=64;
	}
 }


////////////////////////////////////////////////////
 for(int ich=0;ich<NCHANNELS;ich++){
	double ymax = h1[ich]->GetMaximum();
	h1[ich]->GetXaxis()->SetRange(h1[ich]->FindFirstBinAbove(0.001*ymax)-5, h1[ich]->FindLastBinAbove(0.001*ymax)+5);
	h1[ich]->Fit(f1,"Q");

	m4sig[ich] = (UShort_t) (f1->GetParameter(1) + 4*fabs(f1->GetParameter(2)));

	if(nasic==2 && ich==63) ich+=64;
 }

 for(int ich=0;ich<NCHANNELS;ich++)
	delete h1[ich];


////////////////////////////////////////////////////
 skimmer skimEv(fname, m4sig);
 for(long unsigned int ien=0;ien<nen;ien++){
	h22->GetEntry(ien);
	skimEv.Fill(rawEv);
 }

 delete h22;
 return 0;
}


