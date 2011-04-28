#define anaChannel 5
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#include "datastructure.h"
#include "io.h"
#include "functions.h"


void AnalyzeData(char *DataFile = "drs4_peds_5buffers.dat", Int_t nevt,
		Int_t startEv = 1, char *PedFile, Int_t DrawExtraGraphs = 0) {


	// Redefine DOMINO Depth in ADC counts
	const Float_t DominoDepthADC = pow(2, DOMINO_DEPTH);

	// open file

	FILE *fdata = OpenDataFile(DataFile);
	struct channel_struct *p;
	struct channel_struct *dep;

	// create histograms
	// create list of histograms for channels and distribution

	TList *DistChList = new TList();
	TH1F *distch; // histo with distribution of cell-charge, for each channel

	TList *DistChSubList = new TList();
	TH1F *distchsub; // histo with distribution of cell-charge, pedestals subtracted, for each channel

	TList *DistCh0SubList = new TList();
	TH1F *distch0sub; // histo with distribution of cell-charge, pedestals subtracted,
	// channel 0 subtracted for each channel

	TList *grPedList = new TList();
	TGraphErrors *grPed; // for each channel, pedestal value and RMS for each cell is plotted

	TList *hCellList = new TList();
	TH1F *hCell; // charge distribution for each cell (DOMINO_NCELL x DOMINO_NCH histos)
	TList *hCellSubList = new TList();
	TH1F *hCellSub; // charge distribution for each cell (DOMINO_NCELL x DOMINO_NCH histos), pedestal subtracted

	TList *hRMSList = new TList();
	TH1F *hRMSdist; // histo with RMS distribution (statistical RMS of distribution)
	TList *hRMSFitList = new TList();
	TH1F *hRMSFitdist; // histo with RMS distribution (RMS of Gaussian fit)

	TList *grDataList = new TList();
	TGraphErrors *grData; // charge-cell and RMS for each cell is plotted

	TList *grDataSubList = new TList();
	TGraphErrors *grDataSub; // pedestal subtracted charge-cell and RMS for each cell is plotted


	for (int h = 0; h < DOMINO_NCH; h++) {
		//
		TString title = "Data Dist channel";
		title += h;
		distch = new TH1F(title, title, DominoDepthADC, 0., DominoDepthADC);
		DistChList->Add(distch);
		//
		TString title = "Data Dist Ped Sub channel";
		title += h;
		distchsub = new TH1F(title, title, DominoDepthADC, -DominoDepthADC/2, DominoDepthADC/2);
		DistChSubList->Add(distchsub);
		//
		TString title = "Data Dist Ped Ch0 Sub channel";
		title += h;
		distch0sub = new TH1F(title, title, DominoDepthADC, -DominoDepthADC/2, DominoDepthADC/2);
		DistCh0SubList->Add(distch0sub);
		//
		TString title = "Pedestal ch";
		title += h;
		grPed = new TGraphErrors(DOMINO_NCELL);
		grPed->SetTitle(title);
		grPedList->Add(grPed);
		//
		TString title = "Data ch";
		title += h;
		grData = new TGraphErrors(DOMINO_NCELL);
		grData->SetTitle(title);
		grDataList->Add(grData);
		//
		// Mean data and RMS for each channel and cell
		TString title = "Data PedSubtracted ch";
		title += h;
		grDataSub = new TGraphErrors(DOMINO_NCELL);
		grDataSub->SetTitle(title);
		grDataSubList->Add(grDataSub);
		//
		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			// data distribution histos
			TString title = "Data ch";
			title += h;
			title += " cell";
			title += ch;
			hCell = new TH1F(title, title, DominoDepthADC, 0., DominoDepthADC);
			hCellList->Add(hCell);
			// data (ped subtracted) distribution histos
			TString title = "Data PedSub ch";
			title += h;
			title += " cell ";
			title += ch;
			hCellSub = new TH1F(title, title, 2 * DominoDepthADC, -1
					* DominoDepthADC, DominoDepthADC);
			hCellSubList->Add(hCellSub);
		}
		// Data-RMS distribution histos
		TString title = "RMSDist channel";
		title += h;
		hRMSdist = new TH1F(title, title, 100, 0, 20.);
		hRMSList->Add(hRMSdist);
		// Data-RMS (calculated through a fit) distribution histos
		TString title = "RMSFitDist channel";
		title += h;
		hRMSFitdist = new TH1F(title, title, 100, 0, 20.);
		hRMSFitList->Add(hRMSFitdist);
	}
	//--------------
	//
	// calculate or read pedestals from file
	grPedList = OpenPedestals(PedFile);

	//    return;
	//
	// ====== Read data file and subtract the pedestals
	//
	// Count number of events in data file
	int nevtDataMax = 0;
	while (!feof(fdata)) {
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		nevtDataMax++;
	}
	printf("nevtDataMax: %d\n", nevtDataMax);

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
	// filling
	ievt = 1;
	Int_t flagEnd = 0;
	Double_t chtmp;
	Double_t PedVal, itmp, Ch0Val;
	// loop on events
	cout << endl << " --- read DATA file:" << fdata << endl;
	while (ievt <= nevt && !flagEnd) {
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		if (feof(fdata))
			flagEnd = 1;
		if (ievt % (nevt / 10 + 1) == 0)
			cout << "*" << endl;
		p = (struct channel_struct *) &event_data.ch[0]; // read bunch of data
		dep = (struct channel_struct *) &event_data.ch[1]; // read bunch of data

		TGraphErrors *grCh0 = new TGraphErrors(DOMINO_NCELL);

		// loop on channels
		for (int h = 0; h < DOMINO_NCH; h++) {
			// loop on cells
			distch = (TH1F *) DistChList->At(h);
			distchsub = (TH1F *) DistChSubList->At(h);
			grPed = (TGraphErrors *) grPedList->At(h);
			distch0sub = (TH1F *) DistCh0SubList->At(h);
			if(h==0) {
				for(i = 0; i < DOMINO_NCELL;i++) {
					grPed->GetPoint(i, itmp, PedVal);
					chtmp = (Double_t)(p->data[i]);
					chtmp = chtmp - PedVal;
					grCh0->SetPoint(i,itmp, chtmp);
				}
			}
			for (int i = 0; i < DOMINO_NCELL; i++) {
				// Read pedestal value for this cell
				grPed->GetPoint(i, itmp, PedVal);
				grCh0->GetPoint(i, itmp, Ch0Val);
				//                cout << itmp << ", " << PedVal << endl;
				// Read calibration correction for this cell
				//                CalFact =

				//charge distribution for each cell, pedestal subtracted
				chtmp = (Double_t)(p->data[i]); // data value
				//                cout  << "tcell, tcell, depth: " << chtmp << ","  << p->data[i] << "," << deptmp << endl;
				distch->Fill(chtmp);
				// Check data value: must be within DOMINO Depth
				//                if(chtmp > DominoDepthADC)
				//                    cout << " === WARNING!!! Channel " << h << " Cell " << i << " has value " << chtmp << endl;
				//                cout << "Charge: " << p->data[i] << endl;
				((TH1 *) hCellList->At(h * DOMINO_NCELL + i))->Fill(chtmp);
				// Now the pedestal is subtracted
				chtmp = chtmp - PedVal;
				distchsub->Fill(chtmp);
				((TH1 *) hCellSubList->At(h * DOMINO_NCELL + i))->Fill(chtmp);
				chtmp = chtmp - Ch0Val;
				distch0sub->Fill(chtmp);
			}
			p++; // next channel
		}
		ievt++; // next event
	}
	cout << endl;

	// now mean and RMS for each cell are computed and save in histos and graphs
	cout << " --- filling data histos and grphs " << endl;
	TF1 *fgauss = new TF1("fgauss", Gauss, -10., 10., 3);
	fgauss->SetParLimits(0, 0.1, 10000.);
	fgauss->SetParLimits(1, 0., 4096.);
	fgauss->SetParLimits(2, 0.1, 20.);
	Float_t mean, rms, meansub, rmssub;
	for (int h = 0; h < DOMINO_NCH; h++) {
		//        for (int h=5; h<6; h++){
		cout << " Channel:" << h << endl;
		hRMSdist = (TH1F *) hRMSList->At(h);
		hRMSFitdist = (TH1F *) hRMSFitList->At(h);
		grData = (TGraphErrors *) grDataList->At(h);
		grDataSub = (TGraphErrors *) grDataSubList->At(h);
		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			// data distribution histos
			//            cout << "cell:" << ch << " index:" << h*DOMINO_NCELL+ch << " Mean,RMS:"<<hCell->GetMean()<< "," << hCell->GetRMS()<<endl;
			hCell = (TH1F *) hCellList->At(h * DOMINO_NCELL + ch);
			mean = hCell->GetMean();
			rms = hCell->GetRMS();
			hCellSub = (TH1F *) hCellSubList->At(h * DOMINO_NCELL + ch);
			meansub = hCellSub->GetMean();
			rmssub = hCellSub->GetRMS();
			fgauss->SetParameter(0, (Double_t) nevt / 4.);
			fgauss->SetParameter(1, mean);
			fgauss->SetParameter(2, rms);
			//            hCell->Fit("fgauss","QN0");
			grData->SetPoint(ch, ch, mean);
			grData->SetPointError(ch, 0, rms);
			grDataSub->SetPoint(ch, ch, meansub);
			//            grDataSub->SetPointError(ch,0.5,rmssub);
			grDataSub->SetPointError(ch, 0.5, 2.1);
			hRMSdist->Fill(rms);
			hRMSFitdist->Fill(fgauss->GetParameter(2));
			//           cout << "cell:" << ch << " index:" << h*DOMINO_NCELL+ch << " Mean,RMS:"<< mean << "," << rms<<endl;
		}
	}

	Double_t x, y, chtmp, x1, x2, y1, y2;

	/*TList *grCellCalibList = OpenCalibFile("CalibrationData1000events.root");

	TGraphErrors *grCellCalib;
	TGraphErrors *grDataSubCalib = new TGraphErrors(DOMINO_NCELL);
	grDataSubCalib->SetTitle("Data after calibration correction");
	grDataSub = (TGraphErrors *) grDataSubList->At(anaChannel);


	for(ch = 0; ch < DOMINO_NCELL; ch++) {
		grCellCalib = ((TGraphErrors *) grCellCalibList->At(ch));
		grCellCalib->Fit("pol3", "Q");
		TF1 *pol3fit = ((TF1 *) grCellCalib->GetFunction("pol3"));
		grDataSub->GetPoint(ch, x, y);
		chtmp = y - (Double_t)(pol3fit->Eval(y/3.25));
		grDataSubCalib->SetPoint(ch, x, chtmp);
	}

	TCanvas *cGrTest = new TCanvas("grTest", "test per vedere i dati", 1000,1000);

	grDataSubCalib->Draw("APEL");*/


	TString Title = "Charge Distribution per channel";
	gStyle->SetOptFit(111);
	TCanvas *cdistch = new TCanvas("cdistch", Title, 1000, 1000);
	cdistch->Divide(3, 3);
	for (int i = 0; i < DOMINO_NCH; i++) {
		cdistch->cd(i + 1);
		TH1 *dhist = (TH1 *) DistChList->At(i);
		dhist->DrawCopy();
		dhist->SetLineWidth(1);
		dhist->Fit("gaus", "Q");
		dhist->GetFunction("gaus")->SetLineColor(4);
		dhist->GetFunction("gaus")->SetLineWidth(2);
	}

	TString Title = "Charge Distribution Pedestals Subtracted per channel";
	TCanvas *cdistchsub = new TCanvas("cdistchsub", Title, 1000, 1000);
	cdistchsub->Divide(3, 3);
	for (int i = 0; i < DOMINO_NCH; i++) {
		cdistchsub->cd(i + 1);
		TH1 *dsubhist = (TH1 *) DistChSubList->At(i);
		dsubhist->DrawCopy();
		dsubhist->SetLineWidth(1);
		dsubhist->Fit("gaus", "Q");
		dsubhist->GetFunction("gaus")->SetLineColor(4);
		dsubhist->GetFunction("gaus")->SetLineWidth(2);
	}

	TString Title = "Charge Distribution Pedestals and Ch0 Subtracted per channel";
	TCanvas *cdistch0sub = new TCanvas("cdistch0sub", Title, 1000, 1000);
	cdistch0sub->Divide(3, 3);
	for (int i = 0; i < DOMINO_NCH; i++) {
		cdistch0sub->cd(i + 1);
		TH1 *dch0subhist = (TH1 *) DistCh0SubList->At(i);
		dch0subhist->DrawCopy();
		dch0subhist->SetLineWidth(1);
		dch0subhist->Fit("gaus", "Q");
		dch0subhist->GetFunction("gaus")->SetLineColor(4);
		dch0subhist->GetFunction("gaus")->SetLineWidth(2);
	}

	TCanvas *cDataSubTest = new TCanvas("cDataSubTest", "Data after pedestal subtraction", 1000, 1000);
	cDataSubTest->Divide(1,8);
	for (h = 0; h< DOMINO_NCH; h++) {
		grDataSub = (TGraphErrors *) grDataSubList->At(h);
		cDataSubTest->cd(h+1);
		grDataSub->GetYaxis()->SetLabelSize(0.06);
		grDataSub->GetXaxis()->SetLabelSize(0.06);
		grDataSub->Draw("APE");
	}

	TCanvas *cDataSubTestCh5 = new TCanvas("cDataSubTestCh5", "Data after pedestal subtraction Ch5", 1200, 800);
	grDataSub = (TGraphErrors *) grDataSubList->At(anaChannel);
	grDataSub->GetYaxis()->SetLabelSize(0.06);
	grDataSub->GetYaxis()->SetTitle("ADC Counts");
	grDataSub->GetXaxis()->SetTitle("Cell");
	grDataSub->GetXaxis()->SetLabelSize(0.06);
	TLine *refval = new TLine(0,350,1024,350);
	refval->SetLineWidth(3);
	refval->SetLineStyle(2);
	refval->SetLineColor(2);
	TLine *i1 = new TLine(121,-50,121,800);
	i1->SetLineStyle(2);
	TLine *i2 = new TLine(291,-50,291,800);
	i2->SetLineStyle(2);
	TLine *i3 = new TLine(461,-50,461,800);
	i3->SetLineStyle(2);
	TLine *i4 = new TLine(632,-50,632,800);
	i4->SetLineStyle(2);
	TLine *i5 = new TLine(803,-50,803,800);
	i5->SetLineStyle(2);
	TLine *i6 = new TLine(975,-50,975,800);
	i6->SetLineStyle(2);
	TLine *ireal1 = new TLine(121+20,600,121+20,800);
	ireal1->SetLineWidth(3);
	ireal1->SetLineColor(4);
	TLine *ireal2 = new TLine(291-20,600,291-20,800);
	ireal2->SetLineWidth(3);
	ireal2->SetLineColor(4);
	TLine *ireal3 = new TLine(461+20,600,461+20,800);
	ireal3->SetLineWidth(3);
	ireal3->SetLineColor(4);
	TLine *ireal4 = new TLine(632-20,600,632-20,800);
	ireal4->SetLineWidth(3);
	ireal4->SetLineColor(4);
	TLine *ireal5 = new TLine(803+20,600,803+20,800);
	ireal5->SetLineWidth(3);
	ireal5->SetLineColor(4);
	TLine *ireal6 = new TLine(975-20,600,975-20,800);
	ireal6->SetLineWidth(3);
	ireal6->SetLineColor(4);
	grDataSub->Draw("APE");
	refval->Draw("SAME");
	i1->Draw("SAME");
	i2->Draw("SAME");
	i3->Draw("SAME");
	i4->Draw("SAME");
	i5->Draw("SAME");
	i6->Draw("SAME");
	ireal1->Draw("SAME");
	ireal2->Draw("SAME");
	ireal3->Draw("SAME");
	ireal4->Draw("SAME");
	ireal5->Draw("SAME");
	ireal6->Draw("SAME");


	TCanvas *cDataTest = new TCanvas("cDataTest", "Raw Data", 1000,1000);
	cDataTest->Divide(1,8);
	for(h = 0; h < DOMINO_NCH; h++) {
		cDataTest->cd(h+1);
		grData = (TGraphErrors *) grDataList->At(h);
		grData->SetMarkerStyle(20);
		grData->SetMarkerSize(0.5);
		grData->Draw("APE");

	}

	// save root file with graph containing channel 5 data after pedestals subtraction.
	/*
	cout << "test" << endl;

	TString OutFile = DataSubFile;
	TFile *f = new TFile(OutFile,"RECREATE");
	int h = anaChannel;
	TString key="DataSubGraph";
	key += h;
	((TGraphErrors*)grDataSubList->At(h))->Write(key);
	f->Close();
	cout << " ---- Write data on file " << endl;
	 */

	// =======================================================//
	// =====================Matteo's Code=====================//
	// =======================================================//
	/*
	Int_t cht, incCht, decCht, xflag, nPeriods, iMax, iMin;
	Double_t xdiff, incDiff, decDiff, incDiffTemp, decDiffTemp, incXDiff, decXDiff;
	Double_t fitMax, fitMin, fitPeriod, chisquare;
	Double_t DominoXval[DOMINO_NCELL];
	Double_t DominoYval[DOMINO_NCELL];
	Double_t FitXval[DOMINO_NCELL];
	Double_t FitYval[DOMINO_NCELL];


        // opens grDataSub.root

        TString FileName = DataSubFile;
        TGraphErrors *grDataSub;
        int h = anaChannel;
        TFile *f = new TFile(FileName);
        TString key = "DataSubGraph";
        key += h;
        grDataSub = (TGraphErrors *) f->Get(key);
        f->Close();


	// Create a new graph with channel 5 data
	TGraphErrors *grDataSubAnaCh;
	int h = anaChannel;
	grDataSubAnaCh = (TGraphErrors *) grDataSubList->At(h);

	TGraphErrors *grDataSubFix = grDataSubAnaCh->Clone();
	TGraphErrors *grRes = new TGraphErrors(DOMINO_NCELL);
	TList *grResPeriodList = new TList();



	Double_t xtemp, ytemp, DominoMax, DominoMin;

	for (int ch = 0; ch < DOMINO_NCELL; ch++){
		// get domino-output point and save in array
		grDataSubAnaCh->GetPoint(ch, DominoXval[ch], DominoYval[ch]);
	}

	// find the domino point with max y-value
	iMax = 0;
	for(int ch = 0; ch < DOMINO_NCELL; ch++) {
		if(DominoYval[ch] > DominoYval[iMax]) {
			DominoMax = DominoYval[ch];
			iMax = ch;
		}
	}

	cout << "DominoMax e': " << DominoMax << endl;

	// find the domino point with min y-value
	iMin = 0;
	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		if(DominoYval[ch] < DominoYval[iMin]) {
			DominoMin = DominoYval[ch];
			iMin = ch;
		}
	}

	cout << "DominoMin e': " << DominoMin << endl;

	// remove points from the graph that will be used for fit
	for (int ch = 0; ch < DOMINO_NCELL; ch++){
		grDataSubFix->GetPoint(ch, xtemp, ytemp);
		if(ytemp > 0.8*DominoMax || ytemp < 0.2*DominoMin)
			grDataSubFix->RemovePoint(ch);
	}


	TF1 *fsin = new TF1("fsin", sigSin, 0., 1024., 4);
	fsin->SetParameters(600., DOMINO_NCELL / 4., 150., 150.);
	fsin->SetParNames("amplitude", "Period", "Phase", "DC-Offset");
	grDataSubFix->Fit("fsin");
	TF1 *fsinFit = grDataSubFix->GetFunction("fsin");
	fsinFit->SetParNames("amplitude", "Period", "Phase", "DC-Offset");
	chisquare = grDataSub->Chisquare(fsinFit);
	cout << "il chi quadro della funzione di fit e' : " << chisquare << endl;

	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		// get Fit-value and save in array
		FitXval[ch] = DominoXval[ch];
		FitYval[ch] = fsinFit->Eval(FitXval[ch]);
	}

	fitPeriod = fsinFit->GetParameter("Period");
	cout << "il periodo della funzione e': " << fitPeriod << endl;

	nPeriods = (Int_t) (DOMINO_NCELL/fitPeriod);
	cout << "il numero di periodi della funzione e': " << nPeriods << endl;

	fitMax = fsinFit->GetMaximum();
	cout << "il massimo della funzione e': " << fitMax << endl;

	fitMin = fsinFit->GetMinimum();
	cout << "il minimo della funzione e': " << fitMin << endl;




	// computes the y difference between the ch-domino point and the i-fit point
	// and stops when the difference changes sign
	//
	// first and last points are not included in the cicle
	//
	// if the fit point y-value is bigger or smaller than the fit function max*0.8 or min*0.2
	// the point is removed

	for (int ch = 1; ch < DOMINO_NCELL - 1; ch++) {

		if(FitYval[ch] > 0.8*fitMax || FitYval[ch] < 0.2*fitMin) {
			grRes->RemovePoint(ch);
			continue;
		}

		incDiff = DominoYval[ch] - FitYval[ch];
		incDiffTemp = DominoYval[ch] - FitYval[ch + 1];

		decDiff = DominoYval[ch] - FitYval[ch];
		decDiffTemp = DominoYval[ch] - FitYval[ch - 1];

		if(abs(incDiffTemp) < abs(incDiff) || (sign(incDiff) != sign(incDiffTemp) && abs(decDiffTemp) > abs(decDiff))) {
			for (int i = ch; i < DOMINO_NCELL; i++, incDiff = incDiffTemp) {
				incDiffTemp = DominoYval[ch] - FitYval[i];

				if (sign(incDiff) != sign(incDiffTemp)) {
					if(abs(incDiffTemp) < abs(incDiff))
						incCht = i;
					else
						incCht = i - 1;
					break;
				}
			}
			xflag = 1;
		}
		else if(abs(decDiffTemp) < abs(decDiff) || (sign(decDiff) != sign(decDiffTemp) && abs(incDiffTemp) > abs(incDiff))) {
			for (int j = ch; j >= 0 ; j--, decDiff = decDiffTemp) {
				decDiffTemp = DominoYval[ch] - FitYval[j];

				if (sign(decDiff) != sign(decDiffTemp)) {
					if(abs(decDiffTemp) < abs(decDiff))
						decCht = j;
					else
						decCht = j + 1;
					break;
				}
			}
			xflag = -1;
		}

		if(xflag == 1)
			xdiff = FitXval[incCht] - DominoXval[ch];
		else
			xdiff = FitXval[decCht] - DominoXval[ch];

		grRes->SetPoint(ch, (Double_t) ch, xdiff);
	}

	cout << "Draw Time Residuals" << endl;
	TString Title = "Time Residuals";
	TCanvas *timeres = new TCanvas("timeres", Title, 1200, 780);
	grRes->SetMarkerStyle(20);
	grRes->SetMarkerSize(0.3);
	grRes->GetYaxis()->SetLabelSize(0.12);
	grRes->GetXaxis()->SetLabelSize(0.12);
	grRes->Draw("APE");


	// The previous graph is now split in N graphs, where N is the number of fit periods

	// this will be needed to set the function phase
	//
    //    iMax = 0;
	//
    //    for(ch = 0; ch < fitPeriod - 1; ch++) {
    //            if(FitYval[ch] > FitYval[iMax]) iMax = ch;
    //    }

	cout << "il primo massimo ha l'indice : " << iMax << endl;

	for (i = 0; i < nPeriods; i++) {
		TGraphErrors *grResPeriod = new TGraphErrors((Int_t) fitPeriod);
		grResPeriodList->Add(grResPeriod);

		for(ch = i*fitPeriod + 1; ch < fitPeriod + (i*fitPeriod); ch++) {

			if(FitYval[ch] > 0.8*fitMax || FitYval[ch] < 0.2*fitMin) {
				grResPeriod->RemovePoint(ch);
				continue;
			}

			incDiff = DominoYval[ch] - FitYval[ch];
			incDiffTemp = DominoYval[ch] - FitYval[ch + 1];

			decDiff = DominoYval[ch] - FitYval[ch];
			decDiffTemp = DominoYval[ch] - FitYval[ch - 1];

			if(abs(incDiffTemp) < abs(incDiff) || (sign(incDiff) != sign(incDiffTemp) && abs(decDiffTemp) > abs(decDiff))) {
				for (int k = ch; k < k*fitPeriod + fitPeriod; k++, incDiff = incDiffTemp) {
					incDiffTemp = DominoYval[ch] - FitYval[k];

					if (sign(incDiff) != sign(incDiffTemp)) {
						if(abs(incDiffTemp) < abs(incDiff))
							incCht = k;
						else
							incCht = k - 1;
						break;
					}
				}
				xflag = 1;
			}
			else if(abs(decDiffTemp) < abs(decDiff) || (sign(decDiff) != sign(decDiffTemp) && abs(incDiffTemp) > abs(incDiff))) {
				for (int j = ch; j > i*fitPeriod; j--, decDiff = decDiffTemp) {
					decDiffTemp = DominoYval[ch] - FitYval[j];

					if (sign(decDiff) != sign(decDiffTemp)) {
						if(abs(decDiffTemp) < abs(decDiff))
							decCht = j;
						else
							decCht = j + 1;
						break;
					}
				}
				xflag = -1;
			}

			if(xflag == 1)
				xdiff = FitXval[incCht] - DominoXval[ch];
			else
				xdiff = FitXval[decCht] - DominoXval[ch];

			grResPeriod->SetPoint(ch - i*fitPeriod, (Double_t) (ch - i*fitPeriod), xdiff);
		}
	}

	TCanvas *timeresperiod = new TCanvas("timeresperiod", "Time Residuals Period", 1200, 780);
	for(i = 0; i < nPeriods; i++) {
		grResPeriod = ((TGraphErrors *) grResPeriodList->At(i));
		grResPeriod->SetMarkerStyle(20);
		grResPeriod->SetMarkerSize(0.3);
		grResPeriod->GetYaxis()->SetLabelSize(0.12);
		grResPeriod->GetXaxis()->SetLabelSize(0.12);
		grResPeriod->Draw("APEsame");
	}

	cout << "Draw Data - Pedestals Subtracted" << endl;
	TString Title = "Average Charge - Pedestal subtracted";
	TCanvas *csubdata = new TCanvas("csubdata", Title, 1200, 780);
	grDataSubAnaCh->SetMarkerStyle(20);
	grDataSubAnaCh->SetMarkerSize(0.3);
	grDataSubAnaCh->GetYaxis()->SetLabelSize(0.12);
	grDataSubAnaCh->GetXaxis()->SetLabelSize(0.12);
	grDataSubAnaCh->Draw("APE");
	fsinFit->Draw("same");
	 */
	// draw extra graphs
	if (DrawExtraGraphs == 1) {
		cout << " ----- DRAW Results ------" << endl;
		//================ DRAW Results ==================

		TCanvas *c = new TCanvas("ctmp", "test", 800, 800);
		c->Divide(3, 3);
		for (int pad = 1; pad < 10; pad++) {
			c->cd(pad);
			((TH1 *) hCellList->At(pad * 512 + 219))->DrawCopy();
			hCellSub = (TH1F *) hCellSubList->At(pad * 512 + 219);
			hCellSub->SetLineColor(2);
			hCellSub->DrawCopy("same");
		}

		cout << "Draw RMS distributions" << endl;
		TString Title = "RMS distributions per channel";
		TCanvas *c4 = new TCanvas("c4", Title, 700, 700);
		c4->Divide(3, 3);
		for (int i = 0; i < DOMINO_NCH; i++) {
			c4->cd(i + 2);
			hRMSdist = (TH1F *) hRMSList->At(i);
			hRMSFitdist = (TH1F *) hRMSFitList->At(i);
			hRMSFitdist->SetLineColor(2);
			hRMSFitdist->DrawCopy();
			hRMSdist->DrawCopy("same");
		}


		TList *grDataCh0SubList = new TList();
		TGraphErrors *grDataCh0Sub;
		for(h = 0; h< DOMINO_NCELL; h++) {
			grDataCh0Sub = new TGraphErrors(DOMINO_NCELL);
			grDataCh0SubList->Add(grDataCh0Sub);
		}

		TGraphErrors *grDataSubCh0 = (TGraphErrors *) grDataSubList->At(6);
		for(h = 0; h < DOMINO_NCH; h++) {
			grDataSub = (TGraphErrors *) grDataSubList->At(h);
			grDataCh0Sub = (TGraphErrors *) grDataCh0SubList->At(h);
			for(ch = 0; ch < DOMINO_NCELL; ch++) {
				grDataSubCh0->GetPoint(ch, x1, y1);
				grDataSub->GetPoint(ch, x2, y2);
				grDataCh0Sub->SetPoint(ch, x1 , y2 - y1);
			}
		}

		TCanvas *cDataCH0Sub = new TCanvas("cDataCH0Sub","cDataCH0Sub", 1000,1000);
		cDataCH0Sub->Divide(1,8);
		for(h = 0; h < DOMINO_NCH; h++) {
			cDataCH0Sub->cd(h+1);
			grDataCh0Sub = (TGraphErrors *) grDataCh0SubList->At(h);
			grDataCh0Sub->GetYaxis()->SetLabelSize(0.12);
			grDataCh0Sub->GetXaxis()->SetLabelSize(0.12);
			grDataCh0Sub->Draw("APEL");
		}



		cout << "Draw Data - Pedestals Subtracted" << endl;
		TString Title = "Average Charge - Pedestal subtracted";
		TCanvas *csubdata = new TCanvas("csubdata", Title, 1000, 1000);
		csubdata->Divide(3,3);

		for(h = 0; h < DOMINO_NCH; h++) {
			csubdata->cd(h+1);
			TString title = "DataSub channel ";
			title += h;
			TH1F *hCellDataSub = new TH1F(title, title, 100, -20, 20);
			grDataSub = (TGraphErrors *) grDataSubList->At(h);
			for(ch = 0; ch < DOMINO_NCELL; ch++) {
				grDataSub->GetPoint(ch, x, y);
				hCellDataSub->Fill(y);
			}
			hCellDataSub->Fit("gaus", "Q");
			hCellDataSub->GetXaxis()->SetTitle("ADC Counts");
			hCellDataSub->GetFunction("gaus")->SetLineColor(4);
			hCellDataSub->DrawCopy();
		}

		cout << "breakpoint" << endl;
		TCanvas *csubdata2 = new TCanvas("csubdata2", "DataSub for every channel", 1000, 1000);
		TString title = "DataSub every channel ";
		TH1F *hCellChDataSubTot = new TH1F(title, title, 100, -20, 20);
		for(h = 0; h < DOMINO_NCH; h++) {
			grDataSub = (TGraphErrors *) grDataSubList->At(h);
			for(ch = 0; ch < DOMINO_NCELL; ch++) {
				grDataSub->GetPoint(ch, x, y);
				hCellChDataSubTot->Fill(y);
			}
			hCellChDataSubTot->Fit("gaus", "Q");
			hCellChDataSubTot->GetXaxis()->SetTitle("ADC Counts");
			hCellChDataSubTot->GetFunction("gaus")->SetLineColor(4);
			hCellChDataSubTot->Draw();
		}

		cout << "Draw Pedestals" << endl;
		TString Title = "Pedestals";
		TCanvas *c2 = new TCanvas("c2", Title, 1050, 780);
		c2->SetBorderMode(0);
		c2->SetBorderSize(0.);
		c2->Divide(1, 8);
		//    gStyle->SetCanvasBorderMode(0.);
		//    gStyle->SetCanvasBorderSize(0.);
		Double_t x, y;
		for (int h = 0; h < DOMINO_NCH; h++) {
			c2->cd(h + 1);
			grPed = ((TGraphErrors *) grPedList->At(h));
			grPed->SetMarkerStyle(20);
			grPed->SetMarkerSize(0.5);
			grPed->GetYaxis()->SetLabelSize(0.12);
			grPed->GetXaxis()->SetLabelSize(0.12);
			//        cout <<  " err:" << grPed->GetErrorY(102) << " " ;
			//        cout << x << "--" << y << endl;
			grPed->Draw("APE");
		}

		cout << "Draw Data - Average charge" << endl;
		TString Title = "Average_Charge";
		TCanvas *cdata = new TCanvas("cdata", Title, 1050, 780);
		cdata->Divide(1, 8);
		Double_t x, y;
		for (int h = 0; h < DOMINO_NCH; h++) {
			cdata->cd(h + 1);
			grData = ((TGraphErrors *) grDataList->At(h));
			grData->SetMarkerStyle(20);
			grData->SetMarkerSize(0.3);
			grData->GetYaxis()->SetLabelSize(0.12);
			grData->GetXaxis()->SetLabelSize(0.12);
			grData->GetPoint(10, x, y);
			//        cout << x << "-" << y << endl;
			grData->Draw("APE");
		}

		cout << "Draw Data - Pedestals Subtracted" << endl;
		TString Title = "Average Charge - Pedestal subtracted";
		TCanvas *csubdata = new TCanvas("csubdata", Title, 1200, 780);
		csubdata->Divide(1, 8);
		TF1 *fsin = new TF1("fsin", sigSin, 0., 1024., 4);
		TH1D *resDist = new TH1D("resDist", "Residuals Signal", 100, -100., 100.);

		cout << "Draw Data - Pedestals Subtracted" << endl;
		TString Title = "Residuals";
		TCanvas *residuals = new TCanvas("residuals", Title, 1200, 780);
		resDist->DrawCopy();

	}

	fclose(fdata);

	hCellList->Delete();
	hCellSubList->Delete();
	hRMSList->Delete();
	hRMSFitList->Delete();
}

