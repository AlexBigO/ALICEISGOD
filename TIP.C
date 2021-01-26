#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdbool.h>
#include <vector>
#include <iostream>
#include <vector>
#include <math.h> 
#include <TRandom.h>



// Exponential law
// par[0] = C
// par[1] = m_0
// par[2] = T
double expo_law(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = par[0]*TMath::Power( par[2]*(par[1]+par[2]), -1)*x[0]*exp(- (m_T - par[1]) / par[2]);
	return res;
}


// Boltzmann distribution
// par[0] = C
// par[1] = m_0
// par[2] = T
double boltzmann(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	// double res = par[0]*par[1]*par[1]*exp(- (m_T - par[1]) / par[2]);    // version avec m_T au carré et pas p_T au carré
	double res = par[0]*x[0]*x[0]*exp(- (m_T - par[1]) / par[2]);
	return res;
}


// Levy-Tsallis distribution
// par[0] = C
// par[1] = m
// par[2] = T
// par[3] = n
double levy(double *x , double *par)
{
	double m_T = sqrt(x[0]*x[0] + par[1]*par[1]);
	double res = par[0]*(par[3]-1)*(par[3]-2)/(par[3]*par[2]*(par[3]*par[2]+par[1]*(par[3]-2))) *x[0]*TMath::Power(1+(m_T - par[1]) / (par[3]*par[2]) , - par[3]);
	return res;
}


// Power law from [5] in biblio
// par[0] = C
// par[1] = p_0
// par[2] = n
double power_law_Five(double *x , double *par)
{
	double res = par[0]*4*(par[2]-1)*(par[2]-2)*TMath::Power(par[2]-3,-2)*TMath::Power(par[1],-2)*x[0]*TMath::Power(1+ x[0] / (par[1]* (par[2]-3)/2 ) , - par[2]);
	return res;
}


// "Power law" distribution  from [7] in biblio
// par[0] = C
// par[1] = p0
// par[2] = n
double power_law_Seven(double *x, double*par) 
{
	double res = par[0]*x[0]*TMath::Power( 1 + TMath::Power(x[0] / par[1] , 2) , - par[2] );
	return res ;
}



int main(){
cout << " début " << endl ;

		


		// Création du fichier root de sortie et récupération du fichier root données
TString outputfilename="result.root" ;
TFile* OutputHisto = new TFile(outputfilename, "RECREATE");
TFile *myFile = new TFile("HEPData-1569102768-v1-root.root");


		
		// Récupération des différents histogrammes 
TDirectoryFile* dirFile = (TDirectoryFile*)myFile->Get("Table 1");


	TH1F* H1_pp=(TH1F*)dirFile->Get("Hist1D_y1");
	TH1F* H2_pp=(TH1F*)dirFile->Get("Hist1D_y2");
	TH1F* H3_pp=(TH1F*)dirFile->Get("Hist1D_y3");
	TH1F* H4_pp=(TH1F*)dirFile->Get("Hist1D_y4");
	TH1F* H5_pp=(TH1F*)dirFile->Get("Hist1D_y5");
	TH1F* H6_pp=(TH1F*)dirFile->Get("Hist1D_y6");
	TH1F* H7_pp=(TH1F*)dirFile->Get("Hist1D_y7");
	TH1F* H8_pp=(TH1F*)dirFile->Get("Hist1D_y8");
	TH1F* H9_pp=(TH1F*)dirFile->Get("Hist1D_y9");
	TH1F* H10_pp=(TH1F*)dirFile->Get("Hist1D_y10");

	TH1F* H1_pp_e1=(TH1F*)dirFile->Get("Hist1D_y1_e1");
	TH1F* H2_pp_e2=(TH1F*)dirFile->Get("Hist1D_y2_e1");
	TH1F* H3_pp_e3=(TH1F*)dirFile->Get("Hist1D_y3_e1");
	TH1F* H4_pp_e4=(TH1F*)dirFile->Get("Hist1D_y4_e1");
	TH1F* H5_pp_e5=(TH1F*)dirFile->Get("Hist1D_y5_e1");
	TH1F* H6_pp_e6=(TH1F*)dirFile->Get("Hist1D_y6_e1");
	TH1F* H7_pp_e7=(TH1F*)dirFile->Get("Hist1D_y7_e1");
	TH1F* H8_pp_e8=(TH1F*)dirFile->Get("Hist1D_y8_e1");
	TH1F* H9_pp_e9=(TH1F*)dirFile->Get("Hist1D_y9_e1");
	TH1F* H10_pp_e10=(TH1F*)dirFile->Get("Hist1D_y10_e1");
	
	TH1F* H1_pp_E1=(TH1F*)dirFile->Get("Hist1D_y1_e2");
	TH1F* H2_pp_E2=(TH1F*)dirFile->Get("Hist1D_y2_e2");
	TH1F* H3_pp_E3=(TH1F*)dirFile->Get("Hist1D_y3_e2");
	TH1F* H4_pp_E4=(TH1F*)dirFile->Get("Hist1D_y4_e2");
	TH1F* H5_pp_E5=(TH1F*)dirFile->Get("Hist1D_y5_e2");
	TH1F* H6_pp_E6=(TH1F*)dirFile->Get("Hist1D_y6_e2");
	TH1F* H7_pp_E7=(TH1F*)dirFile->Get("Hist1D_y7_e2");
	TH1F* H8_pp_E8=(TH1F*)dirFile->Get("Hist1D_y8_e2");
	TH1F* H9_pp_E9=(TH1F*)dirFile->Get("Hist1D_y9_e2");
	TH1F* H10_pp_E10=(TH1F*)dirFile->Get("Hist1D_y10_e2");
	
	TH1F* H1_pp_SU1=(TH1F*)dirFile->Get("Hist1D_y1_e3");
	TH1F* H2_pp_SU2=(TH1F*)dirFile->Get("Hist1D_y2_e3");
	TH1F* H3_pp_SU3=(TH1F*)dirFile->Get("Hist1D_y3_e3");
	TH1F* H4_pp_SU4=(TH1F*)dirFile->Get("Hist1D_y4_e3");
	TH1F* H5_pp_SU5=(TH1F*)dirFile->Get("Hist1D_y5_e3");
	TH1F* H6_pp_SU6=(TH1F*)dirFile->Get("Hist1D_y6_e3");
	TH1F* H7_pp_SU7=(TH1F*)dirFile->Get("Hist1D_y7_e3");
	TH1F* H8_pp_SU8=(TH1F*)dirFile->Get("Hist1D_y8_e3");
	TH1F* H9_pp_SU9=(TH1F*)dirFile->Get("Hist1D_y9_e3");
	TH1F* H10_pp_SU10=(TH1F*)dirFile->Get("Hist1D_y10_e3");
	
		// Changement des noms d'axes et titres
	H1_pp->SetNameTitle(" Distribution Pt 1 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H1_pp->SetXTitle("p_{T} [GeV/c]");
	H1_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H2_pp->SetNameTitle(" Distribution Pt 2 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H2_pp->SetXTitle("p_{T} [GeV/c]");
	H2_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H3_pp->SetNameTitle(" Distribution Pt 3 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H3_pp->SetXTitle("p_{T} [GeV/c]");
	H3_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H4_pp->SetNameTitle(" Distribution Pt 4 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H4_pp->SetXTitle("p_{T} [GeV/c]");
	H4_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H5_pp->SetNameTitle(" Distribution Pt 5 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H5_pp->SetXTitle("p_{T} [GeV/c]");
	H5_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H6_pp->SetNameTitle(" Distribution Pt 6 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H6_pp->SetXTitle("p_{T} [GeV/c]");
	H6_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H7_pp->SetNameTitle(" Distribution Pt 7 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H7_pp->SetXTitle("p_{T} [GeV/c]");
	H7_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H8_pp->SetNameTitle(" Distribution Pt 8 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H8_pp->SetXTitle("p_{T} [GeV/c]");
	H8_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H9_pp->SetNameTitle(" Distribution Pt 9 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H9_pp->SetXTitle("p_{T} [GeV/c]");
	H9_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	H10_pp->SetNameTitle(" Distribution Pt 10 " , " PT distributiond des pions lors de collisions Pb-PB " );
	H10_pp->SetXTitle("p_{T} [GeV/c]");
	H10_pp->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	
	
		// Compilations des différentes erreurs et ajout sur l'histogramme PT

int Nbinx = H1_pp->GetNbinsX();
cout << " number of bin = " << Nbinx << endl;


for(int i = 0; i <= Nbinx  ; i++){
	H1_pp->SetBinError(i, H1_pp_e1->GetBinContent(i) + H1_pp_E1->GetBinContent(i) + H1_pp_SU1->GetBinContent(i) ) ;
	H2_pp->SetBinError(i, H2_pp_e2->GetBinContent(i) + H2_pp_E2->GetBinContent(i) + H2_pp_SU2->GetBinContent(i) ) ;
	H3_pp->SetBinError(i, H3_pp_e3->GetBinContent(i) + H3_pp_E3->GetBinContent(i) + H3_pp_SU3->GetBinContent(i) ) ;
	H4_pp->SetBinError(i, H4_pp_e4->GetBinContent(i) + H4_pp_E4->GetBinContent(i) + H4_pp_SU4->GetBinContent(i) ) ;
	H5_pp->SetBinError(i, H5_pp_e5->GetBinContent(i) + H5_pp_E5->GetBinContent(i) + H5_pp_SU5->GetBinContent(i) ) ;
	H6_pp->SetBinError(i, H6_pp_e6->GetBinContent(i) + H6_pp_E6->GetBinContent(i) + H6_pp_SU6->GetBinContent(i) ) ;
	H7_pp->SetBinError(i, H7_pp_e7->GetBinContent(i) + H7_pp_E7->GetBinContent(i) + H7_pp_SU7->GetBinContent(i) ) ;
	H8_pp->SetBinError(i, H8_pp_e8->GetBinContent(i) + H8_pp_E8->GetBinContent(i) + H8_pp_SU8->GetBinContent(i) ) ;
	H9_pp->SetBinError(i, H9_pp_e9->GetBinContent(i) + H9_pp_E9->GetBinContent(i) + H9_pp_SU9->GetBinContent(i) ) ;
	H10_pp->SetBinError(i, H10_pp_e10->GetBinContent(i) + H10_pp_E10->GetBinContent(i) + H10_pp_SU10->GetBinContent(i) ) ;
	
}


		//Scaling en prépartion du fit

double Scale_Histo[11] ;
Scale_Histo[1]=H1_pp->Integral() ;
Scale_Histo[2]=H2_pp->Integral() ;
Scale_Histo[3]=H3_pp->Integral() ;
Scale_Histo[4]=H4_pp->Integral() ;
Scale_Histo[5]=H5_pp->Integral() ;
Scale_Histo[6]=H6_pp->Integral() ;
Scale_Histo[7]=H7_pp->Integral() ;
Scale_Histo[8]=H8_pp->Integral() ;
Scale_Histo[9]=H9_pp->Integral() ;
Scale_Histo[10]=H10_pp->Integral() ;


H1_pp->Scale(1/H1_pp->Integral());
H2_pp->Scale(1/H2_pp->Integral());
H3_pp->Scale(1/H3_pp->Integral());
H4_pp->Scale(1/H4_pp->Integral());
H5_pp->Scale(1/H5_pp->Integral());
H6_pp->Scale(1/H6_pp->Integral());
H7_pp->Scale(1/H7_pp->Integral());
H8_pp->Scale(1/H8_pp->Integral());
H9_pp->Scale(1/H9_pp->Integral());
H10_pp->Scale(1/H10_pp->Integral());


		// Changement de range sur l'axe X
		
H1_pp->SetAxisRange(0,3,"X");
H2_pp->SetAxisRange(0,5,"X");
H3_pp->SetAxisRange(0,5,"X");
H4_pp->SetAxisRange(0,5,"X");
H5_pp->SetAxisRange(0,5,"X");
H6_pp->SetAxisRange(0,5,"X");
H7_pp->SetAxisRange(0,5,"X");
H8_pp->SetAxisRange(0,5,"X");
H9_pp->SetAxisRange(0,5,"X");
H10_pp->SetAxisRange(0,5,"X");


		// Creating the different fitting functions
		
TF1 *func1 = new TF1("expo_law",expo_law,0,20,3);
func1->SetParameter(0,1);
func1->SetParameter(1,0.1);
func1->SetParameter(2,0.1);
func1->SetParNames("C","m","T");
func1->SetLineColor(kViolet);

TF1 *func2 = new TF1("boltzmann",boltzmann,0,20,3);
func2->SetParameter(0,0.1);
func2->SetParameter(1,0.1);
func2->SetParameter(2,0.1);
func2->SetParNames("C", "m", "T");	
func2->SetLineColor(kGreen);

TF1 *func3 = new TF1("levy",levy,0,20,4);
func3->SetParameter(0,0.1);
func3->SetParameter(1,0.1);
func3->SetParameter(2,0.1);
func3->SetParameter(3,2);
func3->SetParNames("C","m","T","n");	
func3->SetLineColor(kCyan);

TF1 *func4 = new TF1("power_law_Five", power_law_Five,0,20,3);
func4->SetParameter(0,0.1);
func4->SetParameter(1,0.1);
func4->SetParameter(2,1);
func4->SetParNames("C","p0","n");	
func4->SetLineColor(kMagenta);

TF1 *func5 = new TF1("power_law_Seven", power_law_Seven ,0,20,3);
func5->SetParameter(0,0.1);
func5->SetParameter(1,0.1);
func5->SetParameter(2,1);
func5->SetParNames("C","p0","n");	
func5->SetLineColor(kRed);

cout << "check line" << " " << __LINE__ << endl;
		// Sauvegarde des canvas dans l'output-histo-file càd notre result.root

		// Création des canvas
TCanvas* T1 = new TCanvas("Canvas 1" , " CANVAS_1_PP1 " , 200 , 200 );
/*
TCanvas* T2 = new TCanvas("Canvas 2" , " CANVAS_2_PP1 " , 200 , 200 );
TCanvas* T3 = new TCanvas("Canvas 3" , " CANVAS_3_PP1 " , 200 , 200 );
TCanvas* T4 = new TCanvas("Canvas 4" , " CANVAS_4_PP1 " , 200 , 200 );
TCanvas* T5 = new TCanvas("Canvas 5" , " CANVAS_5_PP1 " , 200 , 200 );
TCanvas* T6 = new TCanvas("Canvas 6" , " CANVAS_6_PP1 " , 200 , 200 );
TCanvas* T7 = new TCanvas("Canvas 7" , " CANVAS_7_PP1 " , 200 , 200 );
TCanvas* T8 = new TCanvas("Canvas 8" , " CANVAS_8_PP1 " , 200 , 200 );
TCanvas* T9 = new TCanvas("Canvas 9" , " CANVAS_9_PP1 " , 200 , 200 );
TCanvas* T10 = new TCanvas("Canvas 10" , " CANVAS_10_PP1 " , 200 , 200 );
*/
		// Sauvegarde des canvas dans l'output-histo-file càd notre result.root
		// Appels des différents fits
gStyle->SetOptFit(0);
gStyle->SetOptStat(1111);
	
OutputHisto->cd();
cout << "check line" << " " << __LINE__ << endl;
T1->cd();
H1_pp->Draw("][ E1");
// We fit the data with the different fitting functions

H1_pp->Fit("expo_law","L + I R","SAME HIST",0,20);   
H1_pp->Fit("boltzmann","L + I R","SAME HIST",0,20);
H1_pp->Fit("levy","L + I R","SAME HIST",0,20);
//H1_pp->Fit("power_law_Five","L + I R","SAME HIST",0,20);
H1_pp->Fit("power_law_Seven","L + I R","SAME HIST",0,20);
cout << "check line" << " " << __LINE__ << endl;

H1_pp->GetFunction("expo_law")->SetLineColor(kViolet);
H1_pp->GetFunction("boltzmann")->SetLineColor(kGreen);
H1_pp->GetFunction("levy")->SetLineColor(kCyan);
//H1_pp->GetFunction("power_law_Five")->SetLineColor(kMagenta);
H1_pp->GetFunction("power_law_Seven")->SetLineColor(kRed);

cout << "check line" << " " << __LINE__ << endl;

cout << "check line" << " " << __LINE__ << endl;

		// Configuration de la légende
auto legend = new TLegend(0.4,0.7,0.6,0.9);
//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend->AddEntry(H1_pp,"Data","lep");
legend->AddEntry(func1,"Exponential law fit","l");
legend->AddEntry(func2,"Boltzmann law fit","l");
legend->AddEntry(func3,"Levy-Tsallis fit","l");
legend->AddEntry(func4, "Power law fit from [5]", "l");
legend->AddEntry(func5, "Power law fit from [7]", "l");
legend->Draw();
cout << "check line" << " " << __LINE__ << endl;
T1->Write();
cout << "check line" << " " << __LINE__ << endl;

/*
T2->cd();
H2_pp->Draw("][ E1");
H2_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H2_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H2_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H2_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H2_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H2_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T2->Write();

T3->cd();
H3_pp->Draw("][ E1");
H3_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H3_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H3_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H3_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H3_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H3_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T3->Write();
 
T4->cd();
H4_pp->Draw("][ E1");
H4_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H4_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H4_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H4_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H4_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H4_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T4->Write();

T5->cd();
H5_pp->Draw("][ E1");
H5_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H5_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H5_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H5_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H5_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H5_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T5->Write();

T6->cd();
H6_pp->Draw("][ E1");
H6_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H6_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H6_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H6_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H6_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H6_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T6->Write();

T7->cd();
H7_pp->Draw("][ E1");
H7_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H7_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H7_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H7_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H7_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H7_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T7->Write();

T8->cd();
H8_pp->Draw("][ E1");
H8_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H8_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H8_pp->Fit("Plaw","L + I R B","SAME HIST", 0 , 20);
H8_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H8_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H8_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T8->Write();

T9->cd();
H9_pp->Draw("][ E1");
H9_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H9_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H9_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H9_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H9_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H9_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T9->Write();

T10->cd();
H10_pp->Draw("][ E1");
H10_pp->Fit("expo", "L" ,"SAME HIST",0 ,20);
H10_pp->Fit("boltzmann", "L + I M R B" ,"SAME HIST",0 ,20);
H10_pp->Fit("Plaw","L + I R B","SAME HIST",0,20);
H10_pp->GetFunction("Plaw")->SetLineColor(kGreen);
H10_pp->GetFunction("boltzmann")->SetLineColor(kViolet);
H10_pp->GetFunction("expo")->SetLineColor(kCyan);
legend->Draw();

T10->Write();
*/

OutputHisto->Close();




cout << " fin " << endl ;
}
