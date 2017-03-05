#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TFile.h>
#include<TTree.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TH2I.h>
#include<TLine.h>
#include"feio.h"

using RICHfrontend::NCHANNELS;
using RICHfrontend::chan2pix;

int main(int argc, char** argv)
{
 if(argc!=2 || (! TString(argv[1]).Contains(TRegexp(".root$")))){
	std::cerr<<"USAGE: "<<argv[0]<<" \"root file with adc map for one channel only\""<<std::endl;
	exit(222);
 }

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
 TH2I *hadc2adc[NCHANNELS], *zadc2adc[NCHANNELS];
 for(int ipix=0;ipix<NCHANNELS;ipix++){
	hadc2adc[ipix] = new TH2I(Form("hadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix%64+1), 300,0.5,3000.5, 300,0.5,3000.5);
	zadc2adc[ipix] = new TH2I(Form("zadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix%64+1), 750,0.5,1500.5, 750,0.5,1500.5);
 }

 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);

	UShort_t adc1 = key & 0xfff;
	UShort_t adc0 = (key>>12) & 0xfff;
	UChar_t ich = (key>>24) & 0xff;

	UChar_t ipix = chan2pix[ich%64]-1 + (ich/64)*64;

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

 double xmaxbin = zadc2adc[pixel-1]->ProjectionX()->GetMaximumBin();
 delete zadc2adc[pixel-1]->ProjectionX();

 double ymaxbin[NCHANNELS];
 for(int ipad=0;ipad<NCHANNELS;ipad++){
	int ipos = (2-(ipad%24)/8);
	int ipix = ipad%8 + (ipad/24)*8 + ipos*64;

	ymaxbin[ipix] = zadc2adc[ipix]->ProjectionY()->GetMaximumBin();
	delete zadc2adc[ipix]->ProjectionY();
 }

 TString pdfname(Form("adc2adc4pixel%02d.pdf",pixel));
 TCanvas* c1 = new TCanvas("c1","c1", 3000,1000);
 c1->Divide(24,8,0.00001,0.00001);
 c1->Print(pdfname+"[");

 TLine ll;
 ll.SetLineWidth(5);
 ll.SetLineColor(kRed);

/*
 for(int ipad=0;ipad<NCHANNELS;ipad++){
	c1->cd(ipad+1);
	int ipos = (2-(ipad%24)/8);
	int ipix = ipad%8 + (ipad/24)*8 + ipos*64;

	hadc2adc[ipix]->Draw("colz");
 }
 c1->Print(pdfname);

 for(int ipad=0;ipad<NCHANNELS;ipad++){
	c1->cd(ipad+1);
	int ipos = (2-(ipad%24)/8);
	int ipix = ipad%8 + (ipad/24)*8 + ipos*64;

	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-50, xmaxbin+500);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-50, ymaxbin[ipix]+500);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->Print(pdfname);
*/

 for(int ipad=0;ipad<NCHANNELS;ipad++){
	c1->cd(ipad+1);
	int ipos = (2-(ipad%24)/8);
	int ipix = ipad%8 + (ipad/24)*8 + ipos*64;

	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-30, xmaxbin+170);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-30, ymaxbin[ipix]+170);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->cd();
 ll.DrawLineNDC(0.3333,0,0.33333,1);
 ll.DrawLineNDC(0.6666,0,0.66666,1);
 c1->Print(pdfname);

/*
 for(int ipad=0;ipad<NCHANNELS;ipad++){
	c1->cd(ipad+1);
	int ipos = (2-(ipad%24)/8);
	int ipix = ipad%8 + (ipad/24)*8 + ipos*64;

	zadc2adc[ipix]->GetXaxis()->SetRange(xmaxbin-20, xmaxbin+80);
	zadc2adc[ipix]->GetYaxis()->SetRange(ymaxbin[ipix]-20, ymaxbin[ipix]+80);
	zadc2adc[ipix]->Draw("colz");
 }
 c1->Print(pdfname);
*/

 c1->Print(pdfname+"]");

 for(int ipix=0;ipix<NCHANNELS;ipix++)
	delete hadc2adc[ipix], zadc2adc[ipix];
 delete tt, ff;

 return 0;
}

