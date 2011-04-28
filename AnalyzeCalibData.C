#define anaChannel 7
#define trigChannel 0
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#define NCALIBFILES 11
#include "datastructure.h"
#include "io.h"
#include "functions.h"

void AnalyzeCalibData(char *CalibFile="CalibrationDataNew1000events.root") {

	TList *grCellCalibList = OpenCalibFile(CalibFile);

	TGraphErrors *grCellCalib;

	Double_t ADCcountTmp, xTmp, ADCcountTmpErr;

	// Redefine DOMINO Depth in ADC counts
	const Float_t DominoDepthADC = pow(2, DOMINO_DEPTH);

	/*
	// create histogram for 200 mV calibration data to use as conversion factor
	TString title = "ADC counts for 200 mv";
	TH1F *h200mV = new TH1F(title,title,DominoDepthADC, -DominoDepthADC, DominoDepthADC);

	// create histogram for 50 mV calibration data to use as conversion factor
	TString title = "ADC counts for 50 mv";
	TH1F *h50mV = new TH1F(title,title,DominoDepthADC, -DominoDepthADC, DominoDepthADC);

	TH1F *h0mV = new TH1F("h0mV","h0mV",400,-100.,100.);
	// get value for 200mV calibration data and fill histo
	for(ch = 0; ch < DOMINO_NCELL; ch++) {
		((TGraphErrors *) grCellCalibList->At(ch))->GetPoint(9, xTmp, ADCcountTmp);
		h200mV->Fill(ADCcountTmp);
	}

	// get value for 50mV calibration data and fill histo
	for(ch = 0; ch < DOMINO_NCELL; ch++) {
		((TGraphErrors *) grCellCalibList->At(ch))->GetPoint(6, xTmp, ADCcountTmp);
		h50mV->Fill(ADCcountTmp);
	}

	TCanvas *c200mV = new TCanvas("200mV ADC counts","200mV ADC counts",800,600);
	h200mV->Fit("gaus");
	h200mV->DrawCopy();

	TCanvas *c50mV = new TCanvas("50mV ADC counts","50mV ADC counts",800,600);
	h50mV->Fit("gaus");
	h50mV->DrawCopy();

	TGraphErrors *gr100mV = new TGraphErrors(DOMINO_NCELL);
	for(ch =0; ch < DOMINO_NCELL; ch++) {
		((TGraphErrors *) grCellCalibList->At(ch))->GetPoint(7, xTmp, ADCcountTmp);
		ADCcountTmpErr = ((TGraphErrors *) grCellCalibList->At(ch))->GetErrorY(7);
		gr100mV->SetPoint(ch, ch, ADCcountTmp);
		gr100mV->SetPointError(ch, 0., ADCcountTmpErr);
	}
	TGraphErrors *gr0mV = new TGraphErrors(DOMINO_NCELL);
	for(ch =0; ch < DOMINO_NCELL; ch++) {
		((TGraphErrors *) grCellCalibList->At(ch))->GetPoint(5, xTmp, ADCcountTmp);
		ADCcountTmpErr = ((TGraphErrors *) grCellCalibList->At(ch))->GetErrorY(7);
		gr0mV->SetPoint(ch, ch, ADCcountTmp);
		gr0mV->SetPointError(ch, 0., ADCcountTmpErr);
		h0mV->Fill(ADCcountTmp);
	}

	TCanvas *cConvmV = new TCanvas("ConvMv","ConvMv", 1000,800);
	cConvmV->Divide(2,2);
	cConvmV->cd(1);
	gr100mV->Draw("APE");
	cConvmV->cd(3);
	gr0mV->Draw("APE");
	cConvmV->cd(4);
	h0mV->Fit("gaus");
	h0mV->Draw(""); */

	// fill histo with calibration fit slopes to find ADC Counts - mV conversion factor

	TH1F *hADCmV = new TH1F("ADCmV","ADC Counts - mV conversion factor",400, 0., 10.);
	//TGraphErrors *grADCmV = new TGraphErrors(DOMINO_NCELL);

	for(int ch = 0; ch < DOMINO_NCELL; ch++) {
		((TGraphErrors *) grCellCalibList->At(ch))->Fit("pol3", "0Q+");
		((TGraphErrors *) grCellCalibList->At(ch))->Fit("pol1", "0Q+", "", -200, 100);
		Float_t p1 = ((TF1 *)((TGraphErrors *) grCellCalibList->At(ch))->GetFunction("pol1"))->GetParameter(1);
		Float_t errp1 = ((TF1 *)((TGraphErrors *) grCellCalibList->At(ch))->GetFunction("pol1"))->GetParError(1);
		hADCmV->Fill(p1);
		//grADCmV->SetPoint(ch, ch, p1);
		//grADCmV->SetPointError(ch, 0., errp1);
	}

	// draw histo with conversion factor

	TCanvas *cADCmV = new TCanvas("ADCmV","ADC Counts - mV conversion factor", 1000,800);
	//cADCmV->Divide(1,2);
	//cADCmV->cd(1);
	hADCmV->Fit("gaus", "Q");
	hADCmV->GetXaxis()->SetTitle("Calibration fits slopes");
	hADCmV->GetFunction("gaus")->SetLineColor(4);
	hADCmV->Draw();
	//cADCmV->cd(2);
	//grADCmV->SetTitle("Calibration fit slope for every cell");
	//grADCmV->GetXaxis()->SetTitle("Cell");
	//grADCmV->GetXaxis()->SetTitle("Fit Slope");
	//grADCmV->Draw("APE");

	// draw every calibration graph in the same canvas

	TCanvas *cCalibAll = new TCanvas("cCalibAll", "All the calibration graphs", 1200,800);

	grCellCalib = (TGraphErrors *) grCellCalibList->At(0);
	grCellCalib->SetTitle("Calibration data for every cell");
	grCellCalib->SetMarkerStyle(20);
	grCellCalib->SetMarkerSize(0.5);
	grCellCalib->GetXaxis()->SetTitle("mV");
	grCellCalib->GetYaxis()->SetTitle("ADC Counts");
	grCellCalib->GetXaxis()->SetLabelSize(0.05);
	grCellCalib->GetYaxis()->SetLabelSize(0.05);
	grCellCalib->Draw("APX");
	for (ch = 1; ch < DOMINO_NCELL; ch++) {
		grCellCalib = (TGraphErrors *) grCellCalibList->At(ch);
		grCellCalib->SetTitle("Calibration data for every cell");
		grCellCalib->SetMarkerStyle(20);
		grCellCalib->SetMarkerSize(0.5);
		grCellCalib->Draw("PX");
	}


	// calibration residuals for channel 500 with pol1 fit

	TGraphErrors *grCalibResCh500 = new TGraphErrors(NCALIBFILES);

	for(int i = 0; i < NCALIBFILES; i++) {
		TF1 *pol1fit = ((TGraphErrors *) grCellCalibList->At(500))->GetFunction("pol1");
		((TGraphErrors *) grCellCalibList->At(500))->GetPoint(i, xTmp, ADCcountTmp);
		ADCcountTmpErr = ((TGraphErrors *) grCellCalibList->At(500))->GetErrorY(i);
		grCalibResCh500->SetPoint(i, xTmp, (Double_t)(ADCcountTmp - pol1fit->Eval(xTmp)));
		grCalibResCh500->SetPointError(i, 0., ADCcountTmpErr);
	}

	TCanvas *cResidualsCh500 = new TCanvas("CalibResCh500", "Calibration Residuals Cell 500",1000,800);
	grCalibResCh500->SetMarkerStyle(20);
	grCalibResCh500->SetMarkerSize(1);
	grCalibResCh500->GetXaxis()->SetTitle("mV");
	grCalibResCh500->GetYaxis()->SetTitle("Calibration residuals");
	grCalibResCh500->SetTitle("Calibration residuals pol1 cell 500");
	grCalibResCh500->Draw("APEL");

	// calibration residuals for channel 500 with pol3 fit

	TGraphErrors *grCalibRes2Ch500 = new TGraphErrors(NCALIBFILES);

	for(int i = 0; i < NCALIBFILES; i++) {
		TF1 *pol3fit = ((TGraphErrors *) grCellCalibList->At(500))->GetFunction("pol3");
		((TGraphErrors *) grCellCalibList->At(500))->GetPoint(i, xTmp, ADCcountTmp);
		ADCcountTmpErr = ((TGraphErrors *) grCellCalibList->At(500))->GetErrorY(i);
		grCalibRes2Ch500->SetPoint(i, xTmp, (Double_t)(ADCcountTmp - pol3fit->Eval(xTmp)));
		grCalibRes2Ch500->SetPointError(i, 0., ADCcountTmpErr);
	}

	TCanvas *cResiduals2Ch500 = new TCanvas("CalibRes2Ch500", "Calibration Residuals 2 Cell 500",1000,800);
	grCalibRes2Ch500->SetMarkerStyle(20);
	grCalibRes2Ch500->SetMarkerSize(1);
	grCalibRes2Ch500->GetXaxis()->SetTitle("mV");
	grCalibRes2Ch500->GetYaxis()->SetTitle("Calibration residuals");
	grCalibRes2Ch500->SetTitle("Calibration residuals pol3 cell 500");
	grCalibRes2Ch500->Fit("pol0", "Q+");
	grCalibRes2Ch500->GetFunction("pol0")->SetLineColor(4);
	grCalibRes2Ch500->Draw("APEL");

	// calibration residuals for every cell with pol1 fit

	TList *grCalibResList = new TList();
	for(ch = 0; ch < DOMINO_NCELL; ch++) {
		TGraphErrors *grCalibRes = new TGraphErrors(NCALIBFILES);
		grCalibRes->SetTitle("Calibration residuals for every cell");
		grCalibResList->Add(grCalibRes);
	}

	for (ch = 0; ch < DOMINO_NCELL; ch++) {
		TF1 *pol1fit = ((TGraphErrors *) grCellCalibList->At(ch))->GetFunction("pol1");
		for (i = 0; i < NCALIBFILES; i++) {
			((TGraphErrors *) grCellCalibList->At(ch))->GetPoint(i, xTmp, ADCcountTmp);
			grCalibRes = (TGraphErrors *) grCalibResList->At(ch);
			grCalibRes->SetPoint(i, xTmp, (Double_t)(ADCcountTmp - pol1fit->Eval(xTmp)));
		}
	}

	TCanvas *cCalibResAll = new TCanvas("cCalibResAll", "Calibration residuals for every cell pol1",1000,800);

	grCalibRes = (TGraphErrors *) grCalibResList->At(0);
	grCalibRes->SetMarkerStyle(20);
	grCalibRes->SetMarkerSize(0.5);
	grCalibRes->GetXaxis()->SetTitle("mV");
	grCalibRes->GetYaxis()->SetTitle("Calibration residuals");
	grCalibRes->Draw("AP");
	for (ch = 1; ch < DOMINO_NCELL; ch++) {
		grCalibRes = (TGraphErrors *) grCalibResList->At(ch);
		grCalibRes->SetMarkerStyle(20);
		grCalibRes->SetMarkerSize(0.5);
		grCalibRes->GetXaxis()->SetTitle("mV");
		grCalibRes->GetYaxis()->SetTitle("Calibration residuals");
		grCalibRes->Draw("P");
	}

	// calibration residuals for every cell with pol3 fit

	for (ch = 0; ch < DOMINO_NCELL; ch++) {
		TF1 *pol3fit = ((TGraphErrors *) grCellCalibList->At(ch))->GetFunction("pol3");
		for (i = 0; i < NCALIBFILES; i++) {
			((TGraphErrors *) grCellCalibList->At(ch))->GetPoint(i, xTmp, ADCcountTmp);
			grCalibRes = (TGraphErrors *) grCalibResList->At(ch);
			grCalibRes->SetPoint(i, xTmp, (Double_t)(ADCcountTmp - pol3fit->Eval(xTmp)));
		}
	}

	TCanvas *cCalibResAll2 = new TCanvas("cCalibResAll2", "Calibration residuals for every cell pol3",1000,800);

	grCalibRes = (TGraphErrors *) grCalibResList->At(0);
	grCalibRes->SetTitle("Calibration residuals for every cell 2");
	grCalibRes->SetMarkerStyle(20);
	grCalibRes->SetMarkerSize(0.5);
	grCalibRes->GetXaxis()->SetTitle("mV");
	grCalibRes->GetYaxis()->SetTitle("Calibration residuals");
	grCalibRes->Draw("AP");
	for (ch = 1; ch < DOMINO_NCELL; ch++) {
		grCalibRes = (TGraphErrors *) grCalibResList->At(ch);
		grCalibRes->SetMarkerStyle(20);
		grCalibRes->SetMarkerSize(0.5);
		grCalibRes->GetXaxis()->SetTitle("mV");
		grCalibRes->GetYaxis()->SetTitle("Calibration residuals");
		grCalibRes->Draw("P");
	}


	// calibration data for 5 "random" cells draw

	/*
	TString title = "Calibration graphs";
	TCanvas *c = new TCanvas(title, title, 800, 600);
	c->Divide(1, 5);
	for (int ch = 0; ch < 5; ch++) {
		c->cd(ch + 1);
		TString title = "Calibration Data Cell";
		title += ch*100;
		grCellTmp = ((TGraphErrors *) grCellCalibList->At(ch*100));
		TF1 *pol3fit = ((TGraphErrors *) grCellCalibList->At(ch*100))->GetFunction("pol3");
		grCellTmp->SetTitle(title);
		grCellTmp->SetMarkerStyle(20);
		grCellTmp->SetMarkerSize(0.5);
		grCellTmp->GetYaxis()->SetLabelSize(0.12);
		grCellTmp->GetXaxis()->SetLabelSize(0.12);
		grCellTmp->Draw("APE");
		pol3fit->Draw("same");
	}
	 */
	// calibration curves for every cell pol3 fit

	TCanvas *cCalibCurves = new TCanvas("cCalibCurves",
			"calibration curves for every cell", 1200, 800);
	grCellTmp = ((TGraphErrors *) grCellCalibList->At(0));
	TF1 *pol3fit = grCellTmp->GetFunction("pol3");
	TString title = "Calibration Curves";
	pol3fit->SetTitle(title);
	pol3fit->SetLineWidth(1);
	pol3fit->GetXaxis()->SetTitle("mV");
	pol3fit->GetYaxis()->SetTitle("ADC Counts");
	pol3fit->Draw();
	for (int ch =1; ch < DOMINO_NCELL;ch++) {
		grCellTmp = ((TGraphErrors *) grCellCalibList->At(ch));
		pol3fit = grCellTmp->GetFunction("pol3");
		TString title = "Calibration Curves";
		pol3fit->SetLineWidth(1);
		pol3fit->Draw("SAME");
	}

	// calibration curves for every cell pol1 fit

	TCanvas *cCalibLines = new TCanvas("cCalibLines",
			"calibration lines for every cell", 1200, 800);
	grCellTmp = ((TGraphErrors *) grCellCalibList->At(0));
	TF1 *pol1fit = grCellTmp->GetFunction("pol1");
	pol1fit->SetRange(-300,300);
	TString title = "Calibration Lines";
	pol1fit->SetTitle(title);
	pol1fit->SetLineWidth(1);
	pol1fit->GetXaxis()->SetTitle("mV");
	pol1fit->GetYaxis()->SetTitle("ADC Counts");
	pol1fit->Draw();
	for (int ch =1; ch < DOMINO_NCELL;ch++) {
		grCellTmp = ((TGraphErrors *) grCellCalibList->At(ch));
		pol1fit = grCellTmp->GetFunction("pol1");
		pol1fit->SetRange(-300,300);
		TString title = "Calibration Lines";
		pol1fit->SetLineWidth(1);
		pol1fit->Draw("SAME");
	}

	// calibration lines for every cell with slope = first parameter (p1) of pol3 fit

	TCanvas *cCalibLines2 = new TCanvas("cCalibLines2",
			"calibration lines for every cell 2", 1200, 800);
	grCellTmp = ((TGraphErrors *) grCellCalibList->At(0));
	pol3fit = grCellTmp->GetFunction("pol3");
	TString title = "Calibration Lines 2";
	Float_t p0 = pol3fit->GetParameter(0);
	Float_t p1 = pol3fit->GetParameter(1);
	TF1 *p1line = new TF1("p1line", "[0]*x + [1]",-300,300);
	p1line->SetParameter(0,p1);
	p1line->SetParameter(1,p0);
	p1line->SetTitle(title);
	p1line->SetLineWidth(1);
	p1line->GetXaxis()->SetTitle("mV");
	p1line->GetYaxis()->SetTitle("ADC Counts");
	p1line->Draw();
	for (int ch =1; ch < DOMINO_NCELL;ch++) {
		grCellTmp = ((TGraphErrors *) grCellCalibList->At(ch));
		pol3fit = grCellTmp->GetFunction("pol3");
		TString title = "Calibration Lines 2";
		p0 = pol3fit->GetParameter(0);
		p1 = pol3fit->GetParameter(1);
		TF1 *p1line = new TF1("p1line", "[0]*x + [1]",-300,300);
		p1line->SetParameter(0,p1);
		p1line->SetParameter(1,p0);
		p1line->SetLineWidth(1);
		p1line->Draw("SAME");
	}

	// calibration curve for channel 500

	TCanvas *cCH500CalibCurves = new TCanvas("cCH500CalibCurves",
			"cell 500 calibration curves", 1200, 800);
	cCH500CalibCurves->SetGrid();
	grCellTmp = ((TGraphErrors *) grCellCalibList->At(320));
	pol3fit = grCellTmp->GetFunction("pol3");
	pol1fit = grCellTmp->GetFunction("pol1");
	pol1fit->SetRange(-300,300);
	TString title = "Calibration Curve - Cell 500";
	grCellTmp->SetTitle(title);
	grCellTmp->SetMarkerStyle(21);
	grCellTmp->SetMarkerSize(1);
	grCellTmp->GetYaxis()->SetTitle("ADC Counts");
	grCellTmp->GetXaxis()->SetTitle("mV");
	pol3fit->SetLineWidth(2);
	pol3fit->SetLineColor(4);
	grCellTmp->Draw("APE");
	pol3fit->Draw("same");
	pol1fit->SetLineStyle(7);
	pol1fit->SetLineWidth(2);
	pol1fit->SetLineColor(2);
	pol1fit->Draw("same");


}
