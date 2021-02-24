// Illustration of the difference/ratio between the pseudo-rapidity and the rapidity of identified particles, y = f(mPDG, pT, eta) 
// The basic formula can be found under : https://en.wikipedia.org/wiki/Pseudorapidity


// Origin : Antonin Maire, November 2017
// First implementation : 16 nov. 2017, Antonin Maire, antonin.maire@cern.ch
// Last modification : 

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "Riostream.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraphBentErrors.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStopwatch.h"
#include "TColor.h"
#include "TASImage.h"

#endif




static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);



enum gkColor{
    red     = 0,
    orange  = 1,
    yellow  = 2,
    kaki    = 3,
    green   = 4,
    cyan    = 5,
    azure   = 6,
    blue    = 7,
    violet  = 8,
    magenta = 9,
    black   = 10,
    nbColor
};



// Preferred colors and markers
const Int_t fillColors[] = {kRed-10,     kOrange-9,   kYellow-9,     kSpring+1,    kGreen-7,    kCyan-8,    kAzure-8,    kBlue-9,     kViolet-9,    kMagenta-9,   kGray+1     }; // for syst bands
const Int_t colors[]     = {kRed+1 ,     kOrange+1,   kOrange-4,     kSpring+9,    kGreen+2,    kCyan+2,    kAzure+2,    kBlue+1,     kViolet-5,    kMagenta+1,   kBlack      };
const Int_t markers[]    = {kFullCircle, kOpenSquare, kFullDotMedium,kOpenSquare,  kOpenSquare, kFullCross, kFullCircle, kFullCircle, kOpenCircle,  kOpenDiamond, kFullCircle };


enum gkParticle{
    kPDG_piPm       =  0,
    kPDG_Kpm        =  1,
    kPDG_K0s        =  2,
    kPDG_K892       =  3,
    kPDG_p          =  4,
    kPDG_phi1020    =  5,
    kPDG_Lambda     =  6,
    kPDG_XiPm       =  7,
    kPDG_Omega      =  8,
    kPDG_D0         =  9,
    kPDG_Dpm        = 10,
    kPDG_Dstar      = 11,
    kPDG_Ds         = 12,
    kPDG_LambdaC    = 13,
    kPDG_Jpsi       = 14,
    kPDG_Psi2S      = 15,
    kPDG_Bplus      = 16,
    kPDG_B0         = 17,
    kPDG_Bs         = 18,
    kPDG_Upsilon1S  = 19,
    kPDG_Upsilon2S  = 20,
    kPDG_Upsilon3S  = 21,
    kPDG_Wpm        = 22,
    kPDG_Z0         = 23,
    kPDG_deuterium  = 24,
    kPDG_tritium    = 25,
    kPDG_He3        = 26,
    kPDG_He4        = 27,
    nbPart
};


const TString Str_TLatexPartName[] = {"#pi^{#pm}","K^{#pm}", "K^{0}#kern[-1.5]{_{s}}", "K*(#scale[0.8]{892})^{0}", "p"        , "#phi(#scale[0.8]{1020})", "#Lambda", "#Xi^{-}", "#Omega^{-}", "D^{0}", "D^{#pm}", "D*(#scale[0.8]{2010})^{#pm}", "D^{#pm}#kern[-1.5]{_{s}}", "#Lambda^{+}#kern[-1.0]{_{c}}", "J/#psi", "#psi(2S)", "B^{+}",  "B^{0}", "B^{0}#kern[-1.0]{_{s}}", "Y(1S)", "Y(2S)", "Y(3S)",  "W^{#pm}", "Z^{0}", "d"        , "t"      , "^{3}He", "^{4}He" };
const TString Str_TxtPartName[]    = { "piPm",    "Kpm",     "K0s",                    "K892",                     "proton"   , "phi1020",                 "Lambda",  "XiPm",    "Omega"     , "D0"   , "Dpm",     "Dstar",                       "Ds"                      , "LambdaC",                      "Jpsi",   "Psi2S",    "Bplus",  "B0",    "Bs",                     "Y1S",   "Y2S",   "Y3S",    "Wpm",     "Z0"   , "deuterium", "tritium", "He3"   , "He4"    };
const Double_t mPDG[]              = { 0.13957018, 0.493677, 0.497611,                 0.89581,                    0.938272081, 1.019461,                  1.115683,  1.31486,   1.67245     , 1.86483, 1.86958,   2.01026,                       1.96827                   , 2.28646,                        3.096900, 3.686097,   5.27931,  5.27962, 5.36682,                  9.64030, 10.02326, 10.3552, 80.385,    91.1876, 1.876123   , 2.809432 , 2.809413, 3.728401 };
    // PDG 2016


Double_t yAsFuncEta(Double_t *eta, Double_t *par);

void myLegendSetUp(TLegend *currentLegend = 0,
                    float currentTextSize = 0.07);

void myPadSetUp(   TPad *currentPad, 
                    float currentLeft   = 0.11, 
                    float currentTop    = 0.04, 
                    float currentRight  = 0.04, 
                    float currentBottom = 0.15);

void myHistoSetUp( TH1 *currentHisto             = 0x0, 
                    Size_t  currentMarkerSize    = 1.0,
                    Style_t currentMarkerStyle   = 20,
                    Color_t currentMarkerColor   = 0,
                    Width_t currentLineWidth     = 1, 
                    Style_t currentLineStyle     = 1, 
                    Color_t currentLineColor     = 0,
                    Style_t currentFillStyle     = 1001,
                    Color_t currentFillColor     = kGray);


void myGraphSetUp( TGraph   *currentGraph      = 0, 
                    Size_t  currentMarkerSize  = 1.0,
                    Style_t currentMarkerStyle = 20,  
                    Color_t currentMarkerColor = 0,
                    Width_t currentLineWidth   = 1,
                    Style_t currentLineStyle   = 1,
                    Color_t currentLineColor   = 0,
                    Style_t currentFillStyle   = 1001,
                    Color_t currentFillColor   = kGray);

void myFuncSetUp(  TF1* f1, 
                    Color_t lColor      = kGray+1, 
                    Style_t lLineStyle  = 1, 
                    Width_t lLineWidth  = 1);

void myOptions(     Int_t lStat = 0);












void Root_ComputeRapidityVsPseudoRapidity(Short_t kPDGpart = kPDG_piPm,  TString rWrite = "0", Bool_t rDrawRatio = 1){
  // rWrite     = "eps", "pdf", "png"
    
    
  myOptions();
  gROOT->ForceStyle();

  TDatime now;
  int iDate = now.GetDate();
  int iYear=iDate/10000;
  int iMonth=(iDate%10000)/100;
  int iDay=iDate%100;
  TString cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"};


    TCanvas *myCanvas[30]       = {0x0};
    TH1F    *myBlankHisto[30]   = {0x0};
    
    vector<TLegend*> myLegend;   myLegend.resize(nbPart);


  
    const Int_t lNbCanvas = 3;
    const Double_t lEtaMin = -2.5;
    const Double_t lEtaMax = 4;
    
      // Display N lines of 2 columns of canvas
    for(Int_t iCan = 0; iCan < lNbCanvas; iCan++){
            myCanvas[iCan] = new TCanvas( Form("myCanvas%d", iCan ), Form("Canvas %02d - [%i %s %i]",iCan, iDay, cMonth[iMonth-1].Data(), iYear), (iCan%3)*550, TMath::Floor(iCan/3)*440, 540,410 );
            myCanvas[iCan]->ToggleEventStatus();
            myCanvas[iCan]->ToggleEditor();
            myCanvas[iCan]->ToggleToolBar();
            myPadSetUp(myCanvas[iCan],0.15,0.04, 0.04,0.16);
            
            myCanvas[iCan]->Draw();
            myCanvas[iCan]->ToggleEventStatus();
            myCanvas[iCan]->ToggleEditor();
            myCanvas[iCan]->ToggleToolBar();
    }


myCanvas[0]->cd();
  myBlankHisto[0] = new TH1F("myBlankHisto[0]","Blank Histogram",650, lEtaMin,lEtaMax);
  myBlankHisto[0]->SetMaximum( lEtaMax );
  myBlankHisto[0]->SetMinimum( lEtaMin );
  myBlankHisto[0]->SetXTitle("Pseudo-rapidity, #eta");
  myBlankHisto[0]->SetYTitle("Rapidity, #it{y}");
  myBlankHisto[0]->SetNdivisions(515,"y");
  myBlankHisto[0]->Draw();
  
  if(rDrawRatio){
    myCanvas[1]->cd();  
        myBlankHisto[1] = new TH1F("myBlankHisto[1]","Blank Histogram",650, lEtaMin,lEtaMax);
        myBlankHisto[1]->SetMaximum(1.5);
        myBlankHisto[1]->SetMinimum(0.01);
        myBlankHisto[1]->SetXTitle("Pseudo-rapidity, #eta");
        myBlankHisto[1]->SetYTitle("Ratio #scale[0.7]{#it{y} / #it{#eta}}");
        myBlankHisto[1]->SetNdivisions(515,"y");
        myBlankHisto[1]->Draw();
        
    myCanvas[2]->cd();  
        myBlankHisto[2] = new TH1F("myBlankHisto[2]","Blank Histogram",650, lEtaMin,lEtaMax);
        myBlankHisto[2]->SetMaximum(2);
        myBlankHisto[2]->SetMinimum(-2);
        myBlankHisto[2]->SetXTitle("Pseudo-rapidity, #eta");
        myBlankHisto[2]->SetYTitle("Difference #scale[0.7]{#it{y} - #it{#eta}}");
        myBlankHisto[2]->SetNdivisions(510,"y");
        myBlankHisto[2]->Draw();    
  }
 
  
  // Compute functions
  
    Int_t lNbPartType = nbPart;
    Int_t lNbPtCases  = 6;
    vector<TString> str_infoMassHyp;   str_infoMassHyp.resize(nbPart);
    
    //-------  
    vector<vector<Double_t> > pTPart;    

        // NOTE : set up the size of the 2D momentum array  (lNbPartType x lNbPtCases)
        pTPart.resize( lNbPartType );
        for (Int_t iPart = 0; iPart < lNbPartType ; ++iPart) {
            pTPart[iPart].resize( lNbPtCases );
        }// end set-up size of the 2D array
    
    
        pTPart[kPDGpart][0] =   0.1;  // in Gev/c
        pTPart[kPDGpart][1] =   0.5;
        pTPart[kPDGpart][2] =   1.0;
        pTPart[kPDGpart][3] =   5.0;
        pTPart[kPDGpart][4] =  10.0;
        pTPart[kPDGpart][5] =  50.0;
    
    
    //-------
    vector<vector<TF1*> > lYfunc;
        // NOTE : set up the size of the 2D TF1 array  (lNbPartType x lNbPtCases)
        lYfunc.resize( lNbPartType );
        for (Int_t iPart = 0; iPart < lNbPartType ; ++iPart) {
            lYfunc[iPart]  .resize( lNbPtCases );
        }// end set-up size of the 2D array    
        
        
        
        
        
  myLegend[kPDGpart] = new TLegend(0.43,0.22,0.90,0.39);
  myLegendSetUp(myLegend[kPDGpart],0.03);
  // myLegend[kPDGpart]->SetHeader( Form("#scale[1.2]{%s}", str_infoMassHyp[kPDGpart].Data() ) );


   for(Int_t iPtIdx = 0; iPtIdx < lNbPtCases; ++iPtIdx ){
        lYfunc[kPDGpart][iPtIdx] = new TF1( Form("yAsFuncEta_ParticleId%02d_ptId%02d", kPDGpart, iPtIdx), yAsFuncEta, -10, 10, 2);
        lYfunc[kPDGpart][iPtIdx]->SetParNames("m0", "pT");
        
        lYfunc[kPDGpart][iPtIdx]->FixParameter(0, mPDG[ kPDGpart ]   ); // mass in GeV/c²
        lYfunc[kPDGpart][iPtIdx]->FixParameter(1, pTPart[kPDGpart ][iPtIdx]); // pT in GeV/c
        
        if(iPtIdx == 0) myFuncSetUp( lYfunc  [kPDGpart][iPtIdx], colors[ red    ], 1, 2);
        if(iPtIdx == 1) myFuncSetUp( lYfunc  [kPDGpart][iPtIdx], colors[ orange ], 1, 2);
        if(iPtIdx == 2) myFuncSetUp( lYfunc  [kPDGpart][iPtIdx], colors[ yellow ], 1, 2);
        if(iPtIdx == 3) myFuncSetUp( lYfunc  [kPDGpart][iPtIdx], colors[ green  ], 1, 2);
        if(iPtIdx == 4) myFuncSetUp( lYfunc  [kPDGpart][iPtIdx], colors[ azure  ], 1, 2);
        if(iPtIdx == 5) myFuncSetUp( lYfunc  [kPDGpart][iPtIdx], colors[ violet ], 1, 2);
        
        myCanvas[0]->cd(); 
        lYfunc[kPDGpart][iPtIdx]->Draw("same");
                
        myLegend[kPDGpart]->AddEntry(lYfunc[kPDGpart][iPtIdx],Form("#it{p}_{T} = %4.1f GeV/#it{c}", pTPart[kPDGpart ][iPtIdx]), "L"  );
   }
    


    myCanvas[0]->cd(); 
    myLegend[kPDGpart]->Draw(); 
    

    
    str_infoMassHyp[kPDGpart].Append( Form("m_{PDG} = m(%s) =  %7.5f GeV/#it{c}^{2}", Str_TLatexPartName[kPDGpart].Data(), mPDG[ kPDGpart ] ) );
    
        
    TLatex *infoMassHyp = new TLatex();
        infoMassHyp->SetNDC();
        //  infoMassHyp->SetTextFont(42);
        infoMassHyp->SetTextSize(0.04);
        infoMassHyp->SetTextColor( kBlack );
        infoMassHyp->DrawLatex( 0.2,0.88, str_infoMassHyp[kPDGpart].Data() );

    
    TString Str_Formula("y = ln #left[ #frac{ #sqrt{m_{PDG}^{2} + p_{#it{T}}^{2}.#it{cosh}^{2} #eta} + p_{#it{T}}.#it{sinh} #eta }{ #sqrt{ m_{PDG}^{2} +  p_{#it{T}}^{2}} }  #right]");
    TLatex *lComment2 = new TLatex();
            lComment2->SetNDC();     
            lComment2->SetTextColor(kGray+3);
            lComment2->SetTextSize(0.03);
            lComment2->SetTextFont(52);
            lComment2->DrawLatex(0.2,0.76, Str_Formula.Data());
  
    myCanvas[0]->cd();   
        TF1 *lf1Diagonal = new TF1( "yEqualToEta", "x", -10, 10 );
        myFuncSetUp( lf1Diagonal, kGray+1, 2, 2); // myFuncSetUp( TF1* f1, Color_t lColor, Style_t lLineStyle, Width_t lLineWidth)
        myLegend[kPDGpart]->AddEntry(lf1Diagonal,        "E >> m_{PDG} #rightarrow  y #approx #eta",      "L");
        lf1Diagonal->Draw("same");            
            
    
    
    if(rDrawRatio){
        
        myCanvas[1]->cd();
            infoMassHyp->DrawLatex( 0.2,0.88, str_infoMassHyp[kPDGpart].Data() );
            lComment2  ->DrawLatex( 0.2,0.76, Str_Formula.Data());
                    
        
        myCanvas[2]->cd();
            infoMassHyp->DrawLatex( 0.47,0.88, str_infoMassHyp[kPDGpart].Data() );
            lComment2  ->DrawLatex( 0.47,0.76, Str_Formula.Data());
        
        
    
        Int_t lNbEtaVal = 1250;
    
        vector<Double_t> lEta;
            lEta.resize( lNbEtaVal );
            for(Int_t iEtaIdx = 0; iEtaIdx < lNbEtaVal; ++iEtaIdx ){
                lEta[iEtaIdx] = lEtaMin + (lEtaMax - lEtaMin)*iEtaIdx/lNbEtaVal;
            }
        
                
        // Vectors for ratio and difference      
        vector<vector<vector<Double_t> > > lRatio; // lRatio[iPart][iPtIdx][iEtaIdx]
        
            lRatio.resize( lNbPartType );
            // NOTE : set up the size of the 3D double array  (lNbPartType x lNbPtCases x lNbEtaVal)
            for (Int_t iPart = 0; iPart < lNbPartType; ++iPart){        
                lRatio[iPart].resize( lNbPtCases );
                for (Int_t iPtIdx = 0; iPtIdx < lNbPtCases ; ++iPtIdx) {
                    lRatio[iPart][iPtIdx].resize( lNbEtaVal );
                }// end set-up sizes
            } // end set-up sizes
            
        vector<vector<vector<Double_t> > > lDiff; // lDiff[iPart][iPtIdx][iEtaIdx]
        
            lDiff.resize( lNbPartType );
            // NOTE : set up the size of the 3D double array  (lNbPartType x lNbPtCases x lNbEtaVal)
            for (Int_t iPart = 0; iPart < lNbPartType; ++iPart){        
                lDiff[iPart].resize( lNbPtCases );
                for (Int_t iPtIdx = 0; iPtIdx < lNbPtCases ; ++iPtIdx) {
                    lDiff[iPart][iPtIdx].resize( lNbEtaVal );
                }// end set-up sizes
            } // end set-up sizes    
            

        
        // Compute the value of ratios y/eta
        for(Int_t iPtIdx = 0; iPtIdx < lNbPtCases; ++iPtIdx){
            Int_t lEtaPbIdx = -1;        
                    
            for(Int_t iEtaIdx = 0; iEtaIdx < lNbEtaVal; ++iEtaIdx){
                
                lDiff [kPDGpart][iPtIdx][iEtaIdx] =  lYfunc[kPDGpart][iPtIdx]->Eval( lEta[iEtaIdx] ) -  lEta[iEtaIdx] ;
                
                // Spot Pbtiq case where eta = 0 --> division by 0 ... should occur at most 1 time per fill
                if( TMath::Abs(lEta[iEtaIdx]) < 1e-13){
                    lEtaPbIdx = iEtaIdx;
                    Printf("Warning : issue with division by 0 : EtaIdx = %04d", iEtaIdx);
                    continue;
                }                        
                // FIXME : rap = TMath::ASinH((lPt[iPtIdx] / mt)  *  TMath::SinH( lEta[iEtaIdx] ));
                lRatio[kPDGpart][iPtIdx][iEtaIdx] =  lYfunc[kPDGpart][iPtIdx]->Eval( lEta[iEtaIdx] ) /  lEta[iEtaIdx] ;
                

            }// end loop iEtaIdx
            
            if( lEtaPbIdx != -1 ) lRatio[kPDGpart][iPtIdx][lEtaPbIdx] = 0.5*( lRatio[kPDGpart][iPtIdx][lEtaPbIdx-1] + lRatio[kPDGpart][iPtIdx][lEtaPbIdx+1] ) ;
        }
            

        

        // Display ratios + diff with TGraph  
        vector<vector<TGraph*> > graphRatio;    
        graphRatio.resize( lNbPartType );    
            for (Int_t iPart = 0; iPart < lNbPartType ; ++iPart) {
                graphRatio[iPart]  .resize( lNbPtCases );
            }// end set-up size of the 2D array    
        
        vector<vector<TGraph*> > graphDiff;    
        graphDiff.resize( lNbPartType );    
            for (Int_t iPart = 0; iPart < lNbPartType ; ++iPart) {
                graphDiff[iPart]  .resize( lNbPtCases );
            }// end set-up size of the 2D array  
        
        
        
        Double_t* lArrEta     = 0;
        Double_t* lArrRatio   = 0;
        Double_t* lArrDiff    = 0;

        
        for(Int_t iPtIdx = 0; iPtIdx < lNbPtCases; ++iPtIdx){
            
            lArrEta   = &lEta[0];
            lArrRatio = &lRatio[kPDGpart][iPtIdx][0];
            lArrDiff  = &lDiff [kPDGpart][iPtIdx][0];
            
            
            graphRatio[kPDGpart][iPtIdx] = new TGraph(lNbEtaVal, lArrEta, lArrRatio );
            graphRatio[kPDGpart][iPtIdx]->SetName( Form("graphRatio_%02d_%02d", kPDGpart, iPtIdx) );
            
            graphDiff[kPDGpart][iPtIdx] = new TGraph(lNbEtaVal, lArrEta, lArrDiff );
            graphDiff[kPDGpart][iPtIdx]->SetName( Form("graphDiff_%02d_%02d", kPDGpart, iPtIdx) );
            
            myGraphSetUp( graphRatio[kPDGpart][iPtIdx], 
                    1,                                          // currentMarkerSize
                    kDot,                                       // currentMarkerStyle,  
                    lYfunc[kPDGpart][iPtIdx]->GetLineColor(),   // Color_t currentMarkerColor,
                    2,                                          // Width_t currentLineWidth  ,           
                    1,                                          // Style_t currentLineStyle  ,     
                    lYfunc[kPDGpart][iPtIdx]->GetLineColor()    // Color_t currentLineColor  
                    );            
            Printf("-- Particle [%d] / Momentum case [%d] / ratio graph name : %s --", kPDGpart, iPtIdx, graphRatio[kPDGpart][iPtIdx]->GetName() ) ;
            
            myGraphSetUp( graphDiff[kPDGpart][iPtIdx], 
                    1,                                          // currentMarkerSize
                    kDot,                                       // currentMarkerStyle,  
                    lYfunc[kPDGpart][iPtIdx]->GetLineColor(),   // Color_t currentMarkerColor,
                    2,                                          // Width_t currentLineWidth  ,           
                    1,                                          // Style_t currentLineStyle  ,     
                    lYfunc[kPDGpart][iPtIdx]->GetLineColor()    // Color_t currentLineColor  
                    );            
            Printf("-- Particle [%d] / Momentum case [%d] / diff  graph name : %s --", kPDGpart, iPtIdx, graphDiff[kPDGpart][iPtIdx]->GetName() ) ;
            
            myCanvas[1]->cd();
            graphRatio[kPDGpart][iPtIdx]->Draw("Csame");
            
            myCanvas[2]->cd();
            graphDiff [kPDGpart][iPtIdx]->Draw("Csame");
            
        }
        

    }// end if rDrawRatio            
                

                
  TString Str_FileName(Form ("Fig_RapVsEta_%02d_%s", kPDGpart, Str_TxtPartName[kPDGpart].Data()) );
  Printf("  --> Generic FileName : %s", Str_FileName.Data() );
  Printf(" ");

  if (rWrite.Contains("png") ) for(Int_t iCan = 0; iCan < lNbCanvas; ++iCan)  myCanvas[iCan]->SaveAs( Form("%s-Can%02d.png", Str_FileName.Data(), iCan) );
  if (rWrite.Contains("pdf") ) for(Int_t iCan = 0; iCan < lNbCanvas; ++iCan)  myCanvas[iCan]->SaveAs( Form("%s-Can%02d.pdf", Str_FileName.Data(), iCan) );
  if (rWrite.Contains("eps") ) for(Int_t iCan = 0; iCan < lNbCanvas; ++iCan)  myCanvas[iCan]->SaveAs( Form("%s-Can%02d.eps", Str_FileName.Data(), iCan) );
  Printf(" ");
    
}



void myLegendSetUp(TLegend *currentLegend,float currentTextSize){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  //currentLegend->SetFillStyle(0);
  currentLegend->SetFillStyle(1001);
  currentLegend->SetFillColor(0);
  currentLegend->SetLineWidth(1);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  currentLegend->SetNColumns(2);  
  
  return;
}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  currentPad->SetFillColor(0);
  currentPad->SetBorderMode(0);
  currentPad->SetBorderSize(2);
   //currentPad->SetLogx();
  // currentPad->SetTicks(); // in myOptions()
  currentPad->SetFrameBorderMode(0);
  currentPad->SetGridx();
  currentPad->SetGridy();
  
  return;
}


void myHistoSetUp( TH1 *currentHisto         , 
                   Size_t  currentMarkerSize ,
                   Style_t currentMarkerStyle,
                   Color_t currentMarkerColor,
                   Width_t currentLineWidth  , 
                   Style_t currentLineStyle  , 
                   Color_t currentLineColor  ,
                   Style_t currentFillStyle  ,
                   Color_t currentFillColor
                 ){
    
    currentHisto->SetStats(0);
    
    currentHisto->SetMarkerSize(currentMarkerSize);
    currentHisto->SetMarkerStyle(currentMarkerStyle);
    currentHisto->SetMarkerColor(currentMarkerColor);
    currentHisto->SetLineWidth(currentLineWidth);
    currentHisto->SetLineStyle(currentLineStyle);
    currentHisto->SetLineColor(currentLineColor);
    currentHisto->SetFillStyle(currentFillStyle);
    currentHisto->SetFillColor(currentFillColor);
    return;
}


void myGraphSetUp(TGraph *currentGraph, 
                  Size_t currentMarkerSize  ,
                  Style_t currentMarkerStyle,  
                  Color_t currentMarkerColor,
                  Width_t currentLineWidth  ,           
                  Style_t currentLineStyle  ,     
                  Color_t currentLineColor  ,
                  Style_t currentFillStyle  ,
                  Color_t currentFillColor  
                 ){
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineWidth(currentLineWidth);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  currentGraph->SetFillStyle(currentFillStyle);
  currentGraph->SetFillColor(currentFillColor);
  return;
}

void myFuncSetUp(TF1* f1, Color_t lColor, Style_t lLineStyle, Width_t lLineWidth){
    f1->SetLineColor(lColor);
    f1->SetLineStyle(lLineStyle);
    f1->SetLineWidth(lLineWidth);    
}

void myOptions(Int_t lStat){
  // Set gStyle
  int font = 42;
  Bool_t kGrayPalette = 0;
  
  // From plain
  gStyle->Reset("Plain");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetGridColor(kGray);

  if(kGrayPalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.0,"xyz");  
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleFillColor(kWhite);  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat("KSiouRMe");
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}



Double_t yAsFuncEta(Double_t *eta, Double_t *par)
{
    
  Double_t lMass    = par[0]; 
  Double_t lpT      = par[1];
  
  // https://en.wikipedia.org/wiki/Pseudorapidity, validated analytically by Antonin on Nov 16 2017
  // See also http://dx.doi.org/10.1155/2013/710534, On current conversion between particle rapidity and pseudorapidity distributions in High Energy Collisions
  
  return TMath::Log( // Neper Log ln(1) =0 
                (
                    TMath::Sqrt( lMass*lMass + lpT*lpT * TMath::Power( TMath::CosH(eta[0]), 2) )
                +   
                    lpT * TMath::SinH(eta[0])
                )
                /
                (  TMath::Sqrt( lMass*lMass + lpT*lpT ) )
            ); 

}
