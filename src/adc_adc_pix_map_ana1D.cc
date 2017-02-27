#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TChain.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TH1I.h>
#include"feio.h"

int main(int argc, char** argv)
{
 UShort_t fadc[RICHfrontend::NCHANNELS];

 TChain* tt = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);
 tt->SetBranchAddress("fadc", fadc);
 int channel = 42;

//////////////////////////////////////////////////////////////
 TH1I *hadc0 = new TH1I("hadc0", "adc for pix", 300,0.5,3000.5);
 TH1I *zadc0 = new TH1I("zadc0", "adc for pix", 1000,0.5,1000.5);

 long unsigned int nen = tt->GetEntries();
 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);

	hadc0->Fill(fadc[channel]);
	zadc0->Fill(fadc[channel]);
 }


//////////////////////////////////////////////////////////////
 TString pdfname(Form("adc4pixel%02d.pdf",RICHfrontend::chan2pix[channel]));
 TCanvas* c1 = new TCanvas("c1","c1", 800,1100);
 c1->Divide(1,2,0.00001,0.00001);
 c1->Print(pdfname+"[");

 c1->cd(1)->SetLogy();
 hadc0->Draw();
 c1->cd(2)->SetLogy();
 zadc0->Draw();
 c1->Print(pdfname);

 c1->Print(pdfname+"]");

 delete hadc0, zadc0;
 delete tt;

 return 0;
}

