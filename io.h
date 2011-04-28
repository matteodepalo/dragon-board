#define anaChannel 5
#define trigChannel 0
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#define NCALIBFILES 11

//======== Routine to open pedestals file ========= //

TList* OpenPedestals(char *PedFile) {
	TString FileName = PedFile;
	TList *grPedList = new TList();
	TGraphErrors *g;
	TFile *f = new TFile(FileName);
	for (int h = 0; h < DOMINO_NCH; h++) {
		TString key = "PedGraph";
		key += h;
		g = (TGraphErrors *) f->Get(key);
		grPedList->Add(g);
	}
	f->Close();
	return grPedList;
}

//============== Routine to open file =============//
FILE *OpenDataFile(char * infile) {
	FILE *fdata;
	int i, j;

	cout << " Input file: " << infile << endl;

	fdata = fopen(infile, "r");
	return fdata;
}

//============== Routine to open period file =============//

Double_t* OpenPeriodFile(char *PeriodFile) {
	FILE *fperiod;
	Double_t* arr = new Double_t[2];
	fperiod = fopen(PeriodFile, "r");
	fread(arr, sizeof(Double_t), 2, fperiod);
	fclose(fperiod);
	return arr;
}

//======== Routine to open calibration file ========= //

TList* OpenCalibFile(char *CalibFile) {
	TString FileName = CalibFile;
	TList *grCellCalibList = new TList();
	TGraphErrors *g;
	TFile *f = new TFile(FileName);
	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		TString key = "CalibDataCell";
		key += ch;
		g = (TGraphErrors *) f->Get(key);
		grCellCalibList->Add(g);
	}
	f->Close();
	return grCellCalibList;
}

//======== Routine to get time residuals graph from TimeResFile file ========= //

TGraphErrors* OpenTimeResFileGraph(char *TimeResGraphFile) {
	TString FileName = TimeResGraphFile;
	TGraphErrors *g;
	TFile *f = new TFile(FileName);
	TString key = "TimeResGraph";
	g = (TGraphErrors *) f->Get(key);
	f->Close();
	return g;
}

//======== Routine to get time residuals graph from TimeResFile file ========= //

TList* OpenTimeResFileHList(char *TimeResHistoFile) {
	TString FileName = TimeResHistoFile;
	TList *hList = new TList();
	TFile *f = new TFile(FileName);
	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		TString key = "TimeResHCell";
		key += ch;
		hCell = (TH1F *) f->Get(key);
		hList->Add(hCell);
	}
	f->Close();
	return hList;
}
