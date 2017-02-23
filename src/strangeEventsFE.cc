#include<iostream>
#include<TChain.h>
#include<TH2F.h>
#include<TH1F.h>
#include<TMultiGraph.h>
#include<TGraph.h>
#include<TStyle.h>
#include<TCanvas.h>
#include"feio.h"

using namespace RICHfrontend;

class strangeRICHEvent:public RICHEvent{
  public:
	strangeRICHEvent();
	~strangeRICHEvent();
	void Fill(rawEvent&);
  private:
	TCanvas* c1;
};

strangeRICHEvent::strangeRICHEvent(){
	gStyle->SetMarkerSize(4);
	c1 = new TCanvas("c1","c1",800,600);
	c1->Print("strange_events.pdf[");
}

void strangeRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);

	for(int ichan=0; ichan<NCHANNELS; ichan++){
		struct edgeMarker{
			int start;
			int end;
		};
		std::vector<edgeMarker> marks;
		bool isbad = false;


		for(int iedge=0; iedge<ftdc[ichan].size(); iedge++){
			if(iedge<ftdc[ichan].size()-1 && ftdc[ichan][iedge] > ftdc[ichan][iedge+1])
				isbad = true;
			
	 		if(fpolar[ichan][iedge]==fLeadingEdge){
				if(iedge<ftdc[ichan].size()-1	&& fpolar[ichan][iedge+1]==fFallingEdge){
					marks.push_back({ftdc[ichan][iedge], ftdc[ichan][iedge+1]});
					iedge++;
				}
				else{
					isbad = true;
					marks.push_back( {ftdc[ichan][iedge], -1} );
				}
			}
			else{
				isbad = true;
				marks.push_back( {-1, ftdc[ichan][iedge]} );
			}
		}

		if(isbad && (ichan==22 || ichan==85 || ichan==170)){
			TGraph goodgr, leadgr, fallgr;

			for(int imark=0;imark<marks.size();imark++){
				if(marks[imark].end<0)
					leadgr.SetPoint(leadgr.GetN(), marks[imark].start, imark);
				else if(marks[imark].start<0)
					fallgr.SetPoint(fallgr.GetN(), marks[imark].end, imark);
				else{
					goodgr.SetPoint(goodgr.GetN(), marks[imark].start, imark);
					goodgr.SetPoint(goodgr.GetN(), marks[imark].end, imark);
				}
			}

			TMultiGraph mgr;
			if(goodgr.GetN()>0){
				goodgr.SetMarkerStyle(20);
				goodgr.SetLineWidth(3);
				mgr.Add(&goodgr, "P");
			}
			if(leadgr.GetN()>0){
				leadgr.SetMarkerStyle(22);
				leadgr.SetMarkerColor(kRed);
				mgr.Add(&leadgr, "P");
			}
			if(fallgr.GetN()>0){
				fallgr.SetMarkerStyle(23);
				fallgr.SetMarkerColor(kBlue);
				mgr.Add(&fallgr, "P");
			}
			mgr.Draw("A*");
			mgr.SetTitle(Form("channel %d, triggerID %d", ichan, trigID));

			c1->Print("strange_events.pdf");
		}
	}
}

strangeRICHEvent::~strangeRICHEvent(){
	c1->Print("strange_events.pdf]");
}




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

 strangeRICHEvent ev;

 int nen = tt->GetEntries();
 for(int ien=0; ien<nen; ien++){
	tt->GetEntry(ien);
	ev.Fill(rawEv);
 }

 return 0;
}

