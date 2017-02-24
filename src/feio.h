#ifndef feio_hh
#define feio_hh


namespace RICHfrontend{
 enum{	NCHANNELS = 192,
		MAXEDGES = 2000};

 struct rawEvent{
	Int_t trigID;
	ULong64_t timeStamp;
	UChar_t hold;

	UShort_t fadc[NCHANNELS];

	UShort_t fnedge;
	UShort_t ftdc[MAXEDGES];
	UShort_t fchan[MAXEDGES];
	Bool_t fpolar[MAXEDGES];
 };


 class RICHEvent{
  protected:
	bool fLeadingEdge;
	bool fTrailingEdge;

	ULong64_t timeStamp;
	Int_t trigID, hold;
	UShort_t fadc[NCHANNELS];

	std::vector<UShort_t> ftdc[NCHANNELS];
	std::vector<Bool_t> fpolar[NCHANNELS];

	std::vector<UShort_t> ftime[NCHANNELS];
	std::vector<UShort_t> fdur[NCHANNELS];

  public:
	RICHEvent():fLeadingEdge(true),fTrailingEdge(false){};

	void Clear(){
		trigID=0;
		timeStamp=0;

		for(int ich=0;ich<NCHANNELS;ich++){
			ftdc[ich].clear();
			fpolar[ich].clear();
			ftime[ich].clear();
			fdur[ich].clear();
		}
	};

	void Fill(rawEvent &rev){
		Clear();

		trigID = rev.trigID;
		timeStamp = rev.timeStamp;
		hold = rev.hold;

		std::copy(rev.fadc, rev.fadc+NCHANNELS, fadc);

		for(int iedge=0; iedge<rev.fnedge; iedge++){
			int ichan = rev.fchan[iedge];
			fpolar[ichan].push_back(rev.fpolar[iedge]);
			ftdc[ichan].push_back(rev.ftdc[iedge]);
		}
	};

	~RICHEvent(){};
 };
}

#endif
