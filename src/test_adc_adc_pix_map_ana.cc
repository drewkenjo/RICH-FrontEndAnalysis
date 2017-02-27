#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TChain.h>
#include<TTree.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TH2I.h>
#include"feio.h"

int main(int argc, char** argv)
{
 UShort_t fadc[RICHfrontend::NCHANNELS];
 TChain* tt = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);
 tt->SetBranchAddress("fadc", fadc);
 long unsigned int nen = tt->GetEntries();

//////////////////////////////////////////////////////////////
 TH1I *hadc[64];
 TH2I *hadc2adc[64], *zadc2adc[64];
 for(int ipix=0;ipix<64;ipix++){
	hadc[ipix] = new TH1I(Form("h1_%02d",ipix), Form("adc for pix %02d",ipix), 1000,0.5,3000.5);
	hadc2adc[ipix] = new TH2I(Form("hadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix), 300,0.5,3000.5, 300,0.5,3000.5);
	zadc2adc[ipix] = new TH2I(Form("zadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix), 750,0.5,1500.5, 750,0.5,1500.5);
 }

 int channel = 42;
 int iasic = channel/64;

 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);

	for(int ich = iasic*64; ich<(iasic+1)*64; ich++){
		hadc[ich]->Fill(fadc[channel]);
		hadc2adc[ich]->Fill(fadc[channel], fadc[ich]);
		zadc2adc[ich]->Fill(fadc[channel], fadc[ich]);
	}
 }

 delete tt;


//////////////////////////////////////////////////////////////
 gStyle->SetLineScalePS(.1);
 gStyle->SetOptStat(0);
 gStyle->SetPadBottomMargin(0.0001);
 gStyle->SetPadTopMargin(0.0001);
 gStyle->SetPadLeftMargin(0.0001);
 gStyle->SetPadRightMargin(0.0001);
 gStyle->SetOptLogz();

 TCanvas* c1 = new TCanvas("c1","c1", 1300,1300);
 c1->Divide(8,8,0.00001,0.00001);
 c1->Print("1.pdf[");

 for(int ipix=0;ipix<64;ipix++){
	c1->cd(ipix+1);
	hadc2adc[ipix]->Draw("colz");
 }
 c1->Print("1.pdf");

 double xmaxbin = zadc2adc[channel]->ProjectionX()->GetMaximumBin();
 double ymaxbin[64];
 for(int ipix=0;ipix<64;ipix++){
	c1->cd(ipix+1);
	ymaxbin[ipix] = zadc2adc[ipix]->ProjectionY()->GetMaximumBin();
	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-50, xmaxbin+500);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-50, ymaxbin[ipix]+500);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->Print("1.pdf");

 for(int ipix=0;ipix<64;ipix++){
	c1->cd(ipix+1);
	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-30, xmaxbin+170);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-30, ymaxbin[ipix]+170);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->Print("1.pdf");

 c1->Print("1.pdf]");

 return 0;
}

