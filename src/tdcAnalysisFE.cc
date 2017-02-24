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
	int ipix[NCHANNELS];
};


goodRICHEvent::goodRICHEvent(){
	int ipixel[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};
	for(int ich=0;ich<NCHANNELS;ich++)
		ipix[ich] = ich/64*64 + ipixel[ich%64] - 1;

	for(int ich=0;ich<NCHANNELS;ich++){
		h2[ipix[ich]]  = new TH2I(Form("h2_%03d",ich), Form("Channel %d, pixel %d; duration; leading time",ich,ipix[ich]%64+1), 500, -0.5, 499.5, 500, -0.5, 499.5);
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

	int ich=170;
	int isdelayed = false;
	for(int isig=0;isig<ftime[ich].size();isig++)
	if(ftime[ich][isig] > -0.293*fdur[ich][isig]+90.25)
		isdelayed = true;

	if(isdelayed)
	for(int ich=0;ich<NCHANNELS;ich++)
	for(int isig=0;isig<ftime[ich].size();isig++){
		h2[ipix[ich]]->Fill(fdur[ich][isig], ftime[ich][isig]);
	}
}


goodRICHEvent::~goodRICHEvent(){
	TCanvas* c1 = new TCanvas("c1","c1",900,900);
	c1->Divide(3,3);
	c1->Print("fe_tdc.pdf[");

	int activeCh[] = {22, 85, 170};
	int activePx[] = {19, 18, 22};
	for(int iactive=0; iactive<3; iactive++){
		for(int irow=-1;irow<2;irow++)
		for(int icol=-1;icol<2;icol++){
			c1->cd(5+irow*3+icol);

			int ipx = activePx[iactive]-1 + 64*iactive + irow*8 + icol;
	
		     h2[ipx]->GetXaxis()->SetRange(1, h2[ipx]->FindLastBinAbove(1, 1)+5);
		     h2[ipx]->GetYaxis()->SetRange(h2[ipx]->FindFirstBinAbove(1, 2)-5, h2[ipx]->FindLastBinAbove(1,2)+5);
			h2[ipx]->Draw("colz");
		}
		c1->Print("fe_tdc.pdf");
	}
	c1->Print("fe_tdc.pdf]");
	
	for(int ipx=0;ipx<NCHANNELS;ipx++)
		delete h2[ipx];
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

