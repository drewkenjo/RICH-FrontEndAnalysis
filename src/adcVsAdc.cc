#include<iostream>
#include<list>
#include<TFile.h>
#include<TTree.h>
#include<TF1.h>
#include<TH2F.h>
#include<TH3F.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(int);
	~goodRICHEvent();
	void Fill(rawEvent&);

  private:
	TH2F* h_adc_adc[NCHANNELS];
	TH2F* z_adc_adc[NCHANNELS];
	TH2F* z_adc_pix[NCHANNELS];
	int nasic;
	uint8_t ipix[64];
};


goodRICHEvent::goodRICHEvent(int _nasic):nasic(_nasic){
	std::copy(chan2pix, chan2pix+64, ipix);

	for(int ich=0;ich<NCHANNELS;ich++){
		h_adc_adc[ich]  = new TH2F(Form("h_adc_adc_%03d",ich), Form("Channel %d, pixel %d;ADC;max ADC in MAPMT",ich,ipix[ich%64]),
			360,350.5,2150.5,		180,350.5,2150.5);
		z_adc_adc[ich]  = new TH2F(Form("z_adc_adc_%03d",ich), Form("Channel %d, pixel %d;ADC;max ADC in MAPMT",ich,ipix[ich%64]),
			500,350.5,850.5,		180,350.5,2150.5);
//			170,450.5,2150.5,		170,450.5,2150.5);
		z_adc_pix[ich]  = new TH2F(Form("z_adc_pix_%03d",ich), Form("Channel %d, pixel %d;ADC;pix# with max ADC in MAPMT",ich,ipix[ich%64]),
			500,350.5,850.5,		64,0.5,64.5);
//			64,0.5,64.5,			400,450.5,850.5);
	}
}


void goodRICHEvent::Fill(rawEvent &rev)
{
	return;
	RICHEvent::Fill(rev);

	struct adcmax{UShort_t fadc; int chan;};
	adcmax fmax[3][2]={{{0,0}, {0,0}}, {{0,0}, {0,0}}, {{0,0}, {0,0}}};
	for(int ichan=0; ichan<NCHANNELS; ichan++){
		int iasic = ichan/64;
		if(fadc[ichan] > fmax[iasic][0].fadc){
			fmax[iasic][1] = fmax[iasic][0];
			fmax[iasic][0] = { fadc[ichan], ichan };
		}
		else if(fadc[ichan] > fmax[iasic][1].fadc)
			fmax[iasic][1] = { fadc[ichan], ichan };
	}

	for(int ichan=0; ichan<NCHANNELS; ichan++){
		int iasic = ichan/64;
		if(ichan!=fmax[iasic][0].chan){
			h_adc_adc[ichan]->Fill(fadc[ichan], fmax[iasic][0].fadc);
			z_adc_adc[ichan]->Fill(fadc[ichan], fmax[iasic][0].fadc);
			z_adc_pix[ichan]->Fill(fadc[ichan], ipix[fmax[iasic][0].chan%64]);
		}
		else{
			h_adc_adc[ichan]->Fill(fadc[ichan], fmax[iasic][1].fadc);
			z_adc_adc[ichan]->Fill(fadc[ichan], fmax[iasic][1].fadc);
			z_adc_pix[ichan]->Fill(fadc[ichan], ipix[fmax[iasic][1].chan%64]);
		}
	}
}


goodRICHEvent::~goodRICHEvent(){
	return;

	TString pdfname("adc_vs_maxadc.pdf");
	TCanvas* c1 = new TCanvas("c1","c1",2000,800);
	c1->Divide(3,2,.0001,.0001);
	c1->Print(pdfname+"[");
	for(int ich=0; ich<NCHANNELS; ich++){
		if(nasic==3 || ich<64 || ich>127){
			TH1* h1 = h_adc_adc[ich]->ProjectionX();

			c1->cd(1)->SetGrid();
			gPad->SetTopMargin(0);
			gPad->SetLogz();
			h_adc_adc[ich]->GetXaxis()->SetRange(h1->FindFirstBinAbove(10)-5, h1->FindLastBinAbove(10));
			h_adc_adc[ich]->Draw("colz");

			c1->cd(4)->SetGrid();
			gPad->SetTopMargin(0);
			h1->GetXaxis()->SetRange(h1->GetMaximumBin()+50, h1->GetNbinsX());
			h1->SetMaximum(h1->GetMaximum()*1.1);
			h1->GetXaxis()->SetRange(h1->FindFirstBinAbove(10)-5, h1->FindLastBinAbove(10));
			h1->Draw();

/////////////////////////////////////////////////////////////////
			TH1* h2 = z_adc_adc[ich]->ProjectionX();

			c1->cd(2)->SetGrid();
			gPad->SetTopMargin(0);
			gPad->SetLogz();
			z_adc_adc[ich]->GetYaxis()->SetRange(z_adc_adc[ich]->FindFirstBinAbove(10,2)-1, z_adc_adc[ich]->FindLastBinAbove(10,2));
			z_adc_adc[ich]->GetXaxis()->SetRange(h2->FindFirstBinAbove(10)-1, h2->FindFirstBinAbove(10)+150);
			z_adc_adc[ich]->Draw("colz");

			c1->cd(5)->SetGrid();
			gPad->SetTopMargin(0);
			gPad->SetLogy();
			h2->GetXaxis()->SetRange(h2->FindFirstBinAbove(10)-1, h2->FindFirstBinAbove(10)+150);
			double mm = h2->GetBinCenter(h2->GetMaximumBin());
			double ss = mm - h2->GetBinCenter(h2->GetXaxis()->GetFirst());
			TF1* f2 = new TF1("f2", "gaus", mm-ss, mm+ss);
			h2->Fit(f2,"QR");
			h2->Draw();

/////////////////////////////////////////////////////////////////
			c1->cd(3)->SetGrid();
			gPad->SetTopMargin(0);
			gPad->SetLogz();
			z_adc_pix[ich]->GetXaxis()->SetRange(h2->FindFirstBinAbove(10)-1, h2->FindFirstBinAbove(10)+150);
			z_adc_pix[ich]->Draw("colz");

			c1->cd(6)->SetGrid();
			gPad->SetTopMargin(0);
			TH1* h3 = (TH1*) h2->Clone("h3");
			for(int ib=1;ib<=h2->GetNbinsX();ib++)
				if(f2->Eval(h3->GetBinCenter(ib))>0.1)
					h3->SetBinContent(ib, h2->GetBinContent(ib)/f2->Eval(h3->GetBinCenter(ib)));
			h3->Draw();

			c1->Print(pdfname);

			delete h1, h2, h3;
		}
		delete h_adc_adc[ich], z_adc_adc[ich], z_adc_pix[ich];
	}
	c1->Print(pdfname+"]");
}


//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev(TString(argv[1]).Contains("2ASIC") ? 2 : 3);

 for(int iarg=1;iarg<argc;iarg++){
	if(!TString(argv[iarg]).Contains(TRegexp(".root$"))) continue;
	TString bname(TString(argv[iarg])(TRegexp("thr[0-9]*")));

	TFile* ff = new TFile(argv[iarg]);
	TTree* tt = (TTree*) ff->Get("h22");

	tt->SetBranchAddress("trigID", &rawEv.trigID);
	tt->SetBranchAddress("timeStamp", &rawEv.timeStamp);

	tt->SetBranchAddress("fadc", rawEv.fadc);

	tt->SetBranchAddress("fnedge", &rawEv.fnedge);
	tt->SetBranchAddress("ftdc", rawEv.ftdc);
	tt->SetBranchAddress("fchan", rawEv.fchan);
	tt->SetBranchAddress("fpolar", rawEv.fpolar);

	int nen = tt->GetEntries();
	for(int ien=0; ien<nen; ien++){
		tt->GetEntry(ien);
		ev.Fill(rawEv);
	}

	delete tt, ff;
 }

 return 0;
}

