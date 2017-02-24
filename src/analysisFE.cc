#include<iostream>
#include<TChain.h>
#include<TH2F.h>
#include<TH1F.h>
#include<TMultiGraph.h>
#include<TGraph.h>
#include<TCanvas.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent();
	~goodRICHEvent();
	void Fill(rawEvent&);
  private:
	TH2I* h2[NCHANNELS];
	TH1I* h1[NCHANNELS];
};

goodRICHEvent::goodRICHEvent(){
	for(int ich=0;ich<NCHANNELS;ich++){
		h1[ich]  = new TH1I(Form("h1_%03d",ich), Form("Channel %d",ich), 4100, -0.5, 4099.5 );
		h2[ich]  = new TH2I(Form("h2_%03d",ich), Form("Channel %d; duration; leading time",ich), 500, -0.5, 499.5, 500, -0.5, 499.5);
	}
}

void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);

	for(int ichan=0; ichan<NCHANNELS; ichan++)
	for(int iedge=0; iedge<ftdc[ichan].size(); iedge++)
	if(fpolar[ichan][iedge]==fLeadingEdge
		&& iedge<ftdc[ichan].size()-1
		&& fpolar[ichan][iedge+1]==fTrailingEdge
		){
			ftime[ichan].push_back(ftdc[ichan][iedge]);
			fdur[ichan].push_back(ftdc[ichan][iedge+1] - ftdc[ichan][iedge]);
			iedge++;
	}

	for(int ich=0;ich<NCHANNELS;ich++){
		h1[ich]->Fill(fadc[ich]);
		for(int isig=0;isig<ftime[ich].size();isig++)
			h2[ich]->Fill(fdur[ich][isig], ftime[ich][isig]);
	}

}

goodRICHEvent::~goodRICHEvent(){
	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	c1->Print("fe_adc.pdf[");

	int activeCh[] = {22, 85, 170};
	for(int iactive=0;iactive<3;iactive++){
		int ich = activeCh[iactive];
	
	     h1[ich]->GetXaxis()->SetRange(h1[ich]->FindFirstBinAbove(1), h1[ich]->FindLastBinAbove(1));
	     h2[ich]->GetXaxis()->SetRange(1, h2[ich]->FindLastBinAbove(1, 1)+5);
	     h2[ich]->GetYaxis()->SetRange(h2[ich]->FindFirstBinAbove(1, 2)-5, h2[ich]->FindLastBinAbove(1,2)+5);
		h2[ich]->Draw("colz");
		c1->Print("fe_adc.pdf");
	}
	c1->Print("fe_adc.pdf]");
	
	for(int ich=0;ich<NCHANNELS;ich++)
		delete h1[ich], h2[ich];
}



//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;

 TChain *tt = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);
 tt->SetBranchAddress("trigID", &rawEv.trigID);
 tt->SetBranchAddress("timeStamp", &rawEv.timeStamp);

 tt->SetBranchAddress("fadc", rawEv.fadc);

 tt->SetBranchAddress("fnedge", &rawEv.fnedge);
 tt->SetBranchAddress("ftdc", rawEv.ftdc);
 tt->SetBranchAddress("fchan", rawEv.fchan);
 tt->SetBranchAddress("fpolar", rawEv.fpolar);

 goodRICHEvent ev;

 int nen = tt->GetEntries();
 for(int ien=0; ien<nen; ien++){
	tt->GetEntry(ien);
	ev.Fill(rawEv);
 }

 return 0;
}

