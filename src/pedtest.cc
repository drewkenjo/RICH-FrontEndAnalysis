#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<TStyle.h>
#include<TCanvas.h>
#include<TH1I.h>

int main(int argc, char** argv)
{
 TH1I* hh[argc][192];
 for(int iarg=1;iarg<argc;iarg++)
 for(int ich=0;ich<192;ich++)
	hh[iarg][ich] = new TH1I(Form("hh%d%03d",iarg,ich), Form("chan %d", ich), 3000,.5,3000.5);

 UShort_t fadc[192];
 for(int iarg=1;iarg<argc;iarg++){
	TFile* ff = new TFile(argv[iarg]);
	TTree* tt = (TTree*) ff->Get("h22");
	tt->SetBranchAddress("fadc",fadc);
	int nen = tt->GetEntries();
	for(int ien=0;ien<nen;ien++){
		tt->GetEntry(ien);
		for(int ich=0;ich<192;ich++)
			hh[iarg][ich]->Fill(fadc[ich]);
	}
	delete tt, ff;
 }

 gStyle->SetOptFit(1111);
 TCanvas* c1 = new TCanvas("c1","c1",1100,800);
 c1->SetLogy();
 c1->Print("ped.pdf[");
 for(int ich=0;ich<192;ich++){
	double xmax=0, xmin=3330, ymax=0;
	for(int iarg=1;iarg<argc;iarg++){
		ymax = std::max(ymax, (double) hh[iarg][ich]->GetMaximum());
		xmax = std::max(xmax, (double) hh[iarg][ich]->FindLastBinAbove(0));
		xmin = std::min(xmin, (double) hh[iarg][ich]->FindFirstBinAbove(0));
	}

	//c1->DrawFrame(xmin-4, .15, xmax+4, ymax*1.1);
	hh[1][ich]->GetXaxis()->SetRange(xmin-4, xmax+4);
	hh[1][ich]->Fit("gaus");
	for(int iarg=2;iarg<argc;iarg++){
		hh[iarg][ich]->SetLineStyle(2);
		hh[iarg][ich]->Draw("same");
	}
	c1->Print("ped.pdf");
 }
 c1->Print("ped.pdf]");
}
