// -----------------------
//  Dragon project
// Analysis tool for test purposes
// read DOMINO4 data;
//  calculate pedestals 
//  pedestal subtraction
//
//  Input:
//   - pedestal file
//   - data file
//  Output:
//   -  a lot of plots
//
//   Structure:
//   1. read pedestal file
//   2. fill charge distribution plot for each cell
//   3. compute pedestal: gaussian fit for each cell -> mean=pedestal
//      - save pedestals on a file
//   4. read data file, subtract pedestal
//   5. signal analysis:
//      - compute average charge for each cell
//      - signal time evolution (?)
//      - signal reconstruction (?)
//
//
// a.stamerra sep.2009
//

#define DOMINO_NCH 8
#define DOMINO_NCELL 1024
#define DOMINO_DEPTH 12  // number of bit
#include "datastructure.h"
#include "io.h"
#include "functions.h"

//============== Routine to open file =============//
FILE *OpenDataFile(char * infile) {
	FILE *fdata;
	int i, j;

	cout << " Input file: " << infile << endl;

	fdata = fopen(infile, "r");
	return fdata;
}

//============== Function Gauss =============//
Double_t Gauss(Double_t* x, Double_t* par) {
	return par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2)) / (par[2]
	                                                                    * sqrt(2 * TMath::Pi()));
}

//============== Signal Function Sin =============//
Double_t sigSin(Double_t* x, Double_t* par) {
	Double_t A = par[0], periodo = par[1], fase = par[2], offset = par[3];
	return offset + A * sin(2 * TMath::Pi() * (x[0] + fase) / periodo);
}

//============== Routine to compute pedestals =============//
TList* CalcPedestals(char *PedFile, Int_t nevtped) {
	// Redefine DOMINO Depth in ADC counts
	const Float_t DominoDepthADC = pow(2, DOMINO_DEPTH);

	// check if pedfile is a data (.dat) or a ROOT (.root) file.
	// If .dat -> compute pedestals, fill graphs and save them on a root file
	// if .root -> read file, get graphs and pedestal values.
	//        the key name in the root file is "PedGraph"
	//
	TString FileName = PedFile;
	TList *grPedList = new TList();
	TGraphErrors *g;
	if (FileName.EndsWith(".root")) {
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

	// ------to be added
	//first check file
	//    if (is a root file)
	//        return NULL;

	// open binary file with pedestal run
	FILE *fped = OpenDataFile(PedFile);
	struct channel_struct *p;

	// read pedFile and count number of events. Set nevtped to this number
	int nevtpedMax = 0;
	while (!feof(fped)) {
		fread((void *) &event_data, 1, sizeof(event_data), fped);
		nevtpedMax++;
	}
	printf("nevtPedMax: %d\n", nevtpedMax);

	if (nevtped > nevtpedMax || nevtped == 0)
		nevtped = nevtpedMax;

	rewind(fped);

	// Create histograms
	TList *hCellPedList = new TList();
	TH1F *hCellPed; // charge distribution for each cell (DOMINO_NCELL x DOMINO_NCH histos)

	for (int h = 0; h < DOMINO_NCH; h++) {
		TString title = "Pedestal ch";
		title += h;
		g = new TGraphErrors(DOMINO_NCELL);
		grPedList->Add(g);
		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			TString title = "Ped ch";
			title += h;
			title += " cell";
			title += ch;
			hCellPed = new TH1F(title, title, DominoDepthADC / 2, 0.,
					DominoDepthADC);
			hCellPedList->Add(hCellPed);
		}
	}

	// Read PedFile and calculate Pedestals for each cell

	// filling
	cout << "----- Pedestals ----" << endl;
	int ievt = 1;
	int flagEnd = 0;
	//    cout << endl << "1st ->";
	// loop on events
	while (ievt <= nevtped && !flagEnd) {
		fread((void *) &event_data, 1, sizeof(event_data), fped);
		if (feof(fped))
			flagEnd = 1;

		if (ievt % (nevtped / 10 + 1) == 0)
			cout << "*" << endl;
		//        cout << "event:" << ievt << " nevt:" << nevtped <<endl;
		//        printf("tcell=%d Depth=%d\n",event_data.tcell,event_data.depth);

		p = (struct channel_struct *) &event_data.ch[0]; // read bunch of data
		// loop on channels
		for (int h = 0; h < DOMINO_NCH; h++) {
			//            printf("%d ",p->data[0]);
			// loop on cells
			for (int i = 0; i < DOMINO_NCELL; i++) {
				//charge distribution for each cell
				((TH1 *) hCellPedList->At(h * DOMINO_NCELL + i))->Fill(
						p->data[i]);
			}
			p++; // next channel
		}
		ievt++; // next event
	}
	cout << endl;

	// Pedestals (mean and RMS obtained from distributions of single cells)
	// cout << endl << " ---- compute pedestals from data" << endl;
	//TF1 *fgauss= new TF1("fgauss",Gauss,-.,10.,3);
	//fgauss->SetParLimits(0, 0.1, 10000.);
	//fgauss->SetParLimits(1, 0., DominoDepthADC);
	//fgauss->SetParLimits(2, 0.1, 20.);
	Float_t mean, rms;
	TH1F *hcelltmp;
	cout << " --- Compute pedestals -----" << endl;
	for (int h = 0; h < DOMINO_NCH; h++) {
		cout << " channel:" << h << endl;
		for (int i = 0; i < DOMINO_NCELL; i++) {
			hcelltmp = (TH1F *) hCellPedList->At(h * DOMINO_NCELL + i);
			//fgauss->SetParameter(0, (Double_t) nevtped / 4.);
			//fgauss->SetParameter(1, hcelltmp->GetMean());
			//fgauss->SetParameter(2, hcelltmp->GetRMS());
			hcelltmp->Fit("gaus", "Q");
			mean = hcelltmp->GetFunction("gaus")->GetParameter(1);
			rms = hcelltmp->GetFunction("gaus")->GetParameter(2);
			//
			((TGraphErrors *) grPedList->At(h))->SetPoint(i, i, mean);
			((TGraphErrors *) grPedList->At(h))->SetPointError(i, 0, rms);
			// cout << mean[h*DOMINO_NCELL+i] << "," << rms[h*DOMINO_NCELL+i] << endl;
		}
	}

	TCanvas *c = new TCanvas("c", "testPed", 800, 800);
	c->Divide(3, 3);
	for (int l = 1; l < 10; l++) {
		c->cd(l);
		((TH1 *) hCellPedList->At(l * 512 + 219))->DrawCopy();
	}

	fclose(fped);
	// save root file with graphs containing the pedestals.
	TString OutFile = PedFile;
	OutFile.Remove(OutFile.Length() - 4);
	OutFile += ".root";
	TFile *f = new TFile(OutFile, "RECREATE");
	for (int h = 0; h < DOMINO_NCH; h++) {
		TString key = "PedGraph";
		key += h;
		((TGraphErrors*) grPedList->At(h))->Write(key);
	}
	f->Close();
	cout << " ---- Write pedestals on file " << OutFile << endl;

	hCellPedList->Delete();

	return grPedList;
}

void AnalyzePed(char *PedFile, Int_t nevtped = 100) {

	const Float_t DominoDepthADC = pow(2, DOMINO_DEPTH);

	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(11);

	TString FileName = PedFile;
	if (!FileName.EndsWith(".root") && !FileName.EndsWith(".dat")) {
		cout << "Pedestal file must end with .dat or .root" << endl;
		return;
	}

	TList *grPedList = new TList();
	TGraphErrors *grPed; // for each channel, pedestal value and RMS for each cell is plotted

	grPedList = CalcPedestals(PedFile, nevtped);


	TCanvas *c1 = new TCanvas("c1", "Pedestals histograms", 1000,1000);
	c1->Divide(3,3);
	Double_t x, y;
	for (int h = 0; h < DOMINO_NCH; h++) {
		TString title = "Pedestals channel ";
		title += h;
		TH1F *hCellPedTmp = new TH1F(title, title, 1000, 1500.,
				2500.);
		grPed = ((TGraphErrors *) grPedList->At(h));
		for (ch = 0; ch < DOMINO_NCELL; ch++) {
			grPed->GetPoint(ch, x, y);
			hCellPedTmp->Fill(y);
		}
		c1->cd(h + 1);
		hCellPedTmp->Fit("gaus");
		((TF1 *) hCellPedTmp->GetFunction("gaus"))->SetLineColor(4);
		((TAxis *) hCellPedTmp->GetXaxis())->SetTitle("ADC Counts");
		hCellPedTmp->DrawCopy();
	}

	cout << "Draw Pedestals" << endl;
	TString Title = "Pedestals";
	TCanvas *c2 = new TCanvas("c2", Title, 1050, 780);
	c2->SetBorderMode(0);
	c2->SetBorderSize(0.);
	c2->Divide(1, 8);
	//    gStyle->SetCanvasBorderMode(0.);
	//    gStyle->SetCanvasBorderSize(0.);

	for (int h = 0; h < DOMINO_NCH; h++) {
		c2->cd(h + 1);
		grPed = ((TGraphErrors *) grPedList->At(h));
		TString Title = "Pedestal channel ";
		Title += h;
		grPed->SetTitle(Title);
		grPed->SetMarkerStyle(20);
		grPed->SetMarkerSize(0.5);
		grPed->GetYaxis()->SetLabelSize(0.12);
		grPed->GetXaxis()->SetLabelSize(0.12);
		//        cout <<  " err:" << grPed->GetErrorY(102) << " " ;
		//        cout << x << "--" << y << endl;
		grPed->Draw("APE");
	}

}

