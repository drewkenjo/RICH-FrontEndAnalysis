#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TFile.h>
#include<TTree.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TH2I.h>
#include"feio.h"

int main(int argc, char** argv)
{
 UInt_t key, val;
 UChar_t chan;
 TFile* ff = new TFile(argv[1]);
 TTree* tt = (TTree*) ff->Get("tmap");
 tt->SetBranchAddress("key", &key);
 tt->SetBranchAddress("val", &val);
 long unsigned int nen = tt->GetEntries();
 TString title(tt->GetTitle());
 int pixel = TString(title(TRegexp("[0-9]*$"))).Atoi();

//////////////////////////////////////////////////////////////
 TH2I *hadc2adc[64], *zadc2adc[64];
 TH1I *zadc0 = new TH1I("zadc0", Form("adc for pix %02d",pixel), 300,200.5,800.5);
 TH1I *zadc1 = new TH1I("zadc1", Form("adc for pix %02d",pixel), 300,200.5,800.5);
 for(int ipix=0;ipix<64;ipix++){
	hadc2adc[ipix] = new TH2I(Form("hadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix+1), 300,0.5,3000.5, 300,0.5,3000.5);
	zadc2adc[ipix] = new TH2I(Form("zadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix+1), 750,0.5,1500.5, 750,0.5,1500.5);
 }

 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);

	UShort_t adc1 = key & 0xfff;
	UShort_t adc0 = (key>>12) & 0xfff;
	UChar_t ipix = ((key>>24) & 0xff)-1;

	if(ipix==pixel){
		zadc0->Fill(adc0, val);
	}

	hadc2adc[ipix]->Fill(adc0, adc1, val);
	zadc2adc[ipix]->Fill(adc0, adc1, val);
 }


//////////////////////////////////////////////////////////////
 gStyle->SetLineScalePS(.7);
 gStyle->SetOptStat(0);
 gStyle->SetPadBottomMargin(0.0001);
 gStyle->SetPadTopMargin(0.0001);
 gStyle->SetPadLeftMargin(0.0001);
 gStyle->SetPadRightMargin(0.0001);
 gStyle->SetOptLogz();

 TString pdfname(Form("adc2adc4pixel%02d.pdf",pixel));
 TCanvas* c1 = new TCanvas("c1","c1", 1300,1300);
 c1->Divide(8,8,0.00001,0.00001);
 c1->Print(pdfname+"[");

 for(int ipix=0;ipix<64;ipix++){
	c1->cd(ipix+1);
	hadc2adc[ipix]->Draw("colz");
 }
 c1->Print(pdfname);

 double xmaxbin = zadc2adc[pixel-1]->ProjectionX()->GetMaximumBin();
 double ymaxbin[64];
 for(int ipix=0;ipix<64;ipix++){
	c1->cd(ipix+1);
	ymaxbin[ipix] = zadc2adc[ipix]->ProjectionY()->GetMaximumBin();
	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-50, xmaxbin+500);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-50, ymaxbin[ipix]+500);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->Print(pdfname);

 for(int ipix=0;ipix<64;ipix++){
	c1->cd(ipix+1);
	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-30, xmaxbin+170);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-30, ymaxbin[ipix]+170);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->Print(pdfname);

 for(int ipix=0;ipix<64;ipix++){
	c1->cd(ipix+1);
	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-20, xmaxbin+80);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-20, ymaxbin[ipix]+80);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->Print(pdfname);

 c1->cd()->SetMargin(0.1,0.1,.1,.1);
 c1->SetLogy();
 zadc0->SetLineWidth(3);
 zadc0->Draw("hist");
 c1->Print(pdfname);

 c1->Print(pdfname+"]");

 for(int ipix=0;ipix<64;ipix++)
	delete hadc2adc[ipix], zadc2adc[ipix];
 delete zadc0, zadc1, c1;
 delete tt, ff;

 return 0;
}

