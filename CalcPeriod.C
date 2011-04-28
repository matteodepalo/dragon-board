#define anaChannel 5
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#include "datastructure.h"
#include "io.h"
#include "functions.h"

const Double_t Pi = TMath::Pi();

void CalcPeriod(char *DataFile = "drs4_peds_5buffers.dat", Int_t nevt,
		Int_t startEv = 1, char *PedFile) {

	// create progress bar
	TGHProgressBar *gProgress = ProgressBar("Calcolo periodo");

	// Redefine DOMINO Depth in ADC counts
	const Float_t DominoDepthADC = pow(2, DOMINO_DEPTH);

	// open file

	FILE *fdata = OpenDataFile(DataFile);
	struct channel_struct *p;
	struct channel_struct *dep;

	// create list of graphs for pedestals
	TList *grPedList = new TList();
	TGraphErrors *grPed;

	// create period histogram

	TString title = "Period histogram";
	TH1 *hPeriod = new TH1F(title,title, 2*((Int_t) DOMINO_NCELL), (Double_t) -DOMINO_NCELL, (Double_t) DOMINO_NCELL);

	// calculate or read pedestals from file
	grPedList = OpenPedestals(PedFile);
	grPed = (TGraphErrors *) grPedList->At(anaChannel);


	// Count number of events in data file
	int nevtDataMax = 0;
	while (!feof(fdata)) {
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		nevtDataMax++;
	}
	printf("nevtDataMax: %d\n", nevtDataMax - 1);

	if (nevt > (nevtDataMax - startEv) || nevt == 0)
		nevt = nevtDataMax - startEv;
	cout << endl << "==>> Processing " << nevt << " events from file "
			<< DataFile << endl;

	rewind(fdata);


	Int_t ievt = 1;
	// go to first event (startEv)
	while (ievt < startEv) {
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		if (feof(fdata))
			break;
		ievt++;
	}

	ievt = 1;
	Int_t flagEnd = 0;
	Int_t fitusati = 0;
	Double_t chtmp;
	Double_t PedVal, itmp;
	Double_t mean, rms;
	Double_t ratio;



	//debug canvas
	TCanvas *cfitTest = new TCanvas("cfitTest", "fit tests", 1200, 780);
	cfitTest->Divide(1,nevt);

	// loop on events

	gProgress->Reset();
	gProgress->SetMax(nevt);

	gSystem->ProcessEvents();

	while (ievt <= nevt && !flagEnd) {
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		if (feof(fdata))
			flagEnd = 1;

		p = (struct channel_struct *) &event_data.ch[0]; // read bunch of data
		dep = (struct channel_struct *) &event_data.ch[1]; // read bunch of data
		// goes to channel to analyze
		p += anaChannel;

		// read data, subtract pedestals values and save results in grAnaChDataTemp graph with
		// fixed error for each point (x = 0.5 and y = 2.1). Also generate an array with Domino
		// X and Y values

		TGraphErrors *grAnaChDataTemp = new TGraphErrors(DOMINO_NCELL);

		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			// Read pedestal value for this cell
			grPed->GetPoint(ch, itmp, PedVal);
			chtmp = (Double_t)(p->data[ch]); // data value
			chtmp = chtmp - PedVal;
			grAnaChDataTemp->SetPoint(ch, (Double_t) ch, chtmp);
			grAnaChDataTemp->SetPointError(ch, 0.5, 2.1);
		}
		// create fit functions
		TF1 *fsin = new TF1("fsin", sigSin, 0., 1024., 4);
		fsin->SetParameters(600., 255., 150., 150.);
		fsin->SetParNames("amplitude", "Period", "Phase", "DC-Offset");

		grAnaChDataTemp->Fit("fsin", "Q");
		TF1 *fsinFit = grAnaChDataTemp->GetFunction("fsin");
		fsinFit->SetParNames("amplitude", "Period", "Phase", "DC-Offset");

		// debug
		cfitTest->cd(ievt);
		grAnaChDataTemp->SetMarkerStyle(20);
		grAnaChDataTemp->SetMarkerSize(0.3);
		grAnaChDataTemp->GetYaxis()->SetLabelSize(0.12);
		grAnaChDataTemp->GetXaxis()->SetLabelSize(0.12);
		grAnaChDataTemp->Draw("APE");

		Double_t fitPeriod, fitAmplitude, chisquare;
		fitPeriod = fsinFit->GetParameter("Period");
		fitAmplitude = TMath::Abs(fsinFit->GetParameter("amplitude"));
		chisquare = fsinFit->GetChisquare();

		cout << "period: " << fitPeriod << " amplitude: " << fitAmplitude << " chisquare: " << chisquare << endl;

		if(chisquare > 0.1e+06) {
			gProgress->Increment(1);
			gSystem->DispatchOneEvent(kTRUE);
			ievt++;
			continue;
		}

		gProgress->Increment(1);
		gSystem->DispatchOneEvent(kTRUE);

		hPeriod->Fill(fitPeriod);
		fitusati++;

		ievt++;

	}

	cout << "fit scartati :" << nevt - fitusati << endl;
	//draw
	TString Title = "Period distribution for nevt events";
	TCanvas *cPeriod = new TCanvas("cPeriod", Title, 700, 700);
	hPeriod->Draw();
	hPeriod->Fit("gaus");

	TF1 *fgausFit = hPeriod->GetFunction("gaus");
	//mean = fgausFit->GetParameter(1);
	//	rms = fgausFit->GetParameter(2);

	mean = hPeriod->GetMean();
	rms = hPeriod->GetRMS();

	TString OutFile = "Period";
	OutFile += nevt;
	OutFile += "events.dat";
	FILE *f = fopen(OutFile.Data(), "w");
	fwrite(&mean, sizeof(mean), 1, f);
	fwrite(&rms, sizeof(rms), 1, f);

	((TGMainFrame *) gProgress->GetParent())->CloseWindow();
	fclose(f);

	cout << "mean: " << mean << " rms: " << rms << endl;

	fclose(fdata);

}

