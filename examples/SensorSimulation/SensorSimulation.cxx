#include <iostream>
#include <sstream>



#include <ratio>
#include <cmath>

#include "TF1.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "tbb/parallel_for.h"

const double c = 3e8;
const double Na = 6.022e23;
const double re = 2.817e-13;
const double alpha = 7.297353e-3;
const std::vector<double> si = {14, 28.0855, 2.33, 173, 0.1492, 0.2014, 2.87, 3.25, -4.44};
const std::vector<double> al = {13, 26.981, 2.7, 166, 0.0802, 0.1708, 3.01, 3.63, -4.24};
const std::vector<double> h2o = {17, 33.006, 1, 75, 0.0911, 0.24, 2.8, 3.48, -3.5};
/* Generate energy spectrum (*)
   (Energy loss through al) ()
   Does it scatter? Mass attenuation (*)
   Generate angle (*)
   Calculate energy (*)
   Store in TTree (*)
   Calculate radius and losses (*)
   Does it hit? ()

 */

class Event
{
public:
	double incPhotonEnergy;
	double electronEnergy;
	double photonAngle;
	double electronAngle;
	double electronKEnergy;
	double electronRadius;

	Event()
	{
		incPhotonEnergy=(0);
		electronEnergy=(0);
		photonAngle=0;
		electronAngle=0;
		electronKEnergy=0;
		electronRadius=0;
	}
};


class BetheValues
{
public:
	double Z;
	double A;
	double rho;
	double I;
	double a;
	double X0;
	double X1;
	double m;
	double C0;

	BetheValues()
	{
		Z=(0);
		A=(0);
		rho=(0);
		I=(0);
		a=(0);
		X0=(0);
		X1=(0);
		m=(0);
		C0=0;
		betheCoefficientMap =  { {"si",si}, {"al",al}, {"h2o",h2o} };
	}


	void SetParameters(std::string material);

private:

	std::map<std::string,  std::vector<double>> betheCoefficientMap =  { {"si",si}, {"al",al}, {"h2o",h2o}  };

};



void BetheValues::SetParameters(std::string material)
{

	Z = betheCoefficientMap[material].at(0);
	A = betheCoefficientMap[material].at(1);
	rho = betheCoefficientMap[material].at(2);
	I = betheCoefficientMap[material].at(3);
	a = betheCoefficientMap[material].at(4);
	X0 = betheCoefficientMap[material].at(5);
	X1 = betheCoefficientMap[material].at(6);
	m = betheCoefficientMap[material].at(7);
	C0= betheCoefficientMap[material].at(8);
}
Double_t bremsstrahlung(Double_t *x, Double_t *par)
{
	Double_t Z = par[0];             //Z number of target
	Double_t spectrumEnergy = par[1];//Energy of spectrum (MeV)
	Double_t energy = x[0];

	return std::sqrt(Z) * ( (spectrumEnergy - energy ) / energy ) * ( -73.9 - 1.2336 * energy + 36.502 * std::log(Z) + ( ( 148.5 * std::pow(spectrumEnergy, 0.1293 ) / Z ) ) ) * ( 1 + ( -0.006624 + 0.0002906 * spectrumEnergy ) * ( Z/energy ) );
	//ROOT function to produce a bremsstrahlung spectrum for the incoming photons.
}

Double_t kleinNishna(Double_t *x, Double_t *par)
{
	Double_t energy = par[0];             //Energy of photon (MeV)
	Double_t theta = x[0];

	Double_t P = 1 / ( 1 + ( energy/0.511 ) * ( 1 - std::cos(theta) ) );

	return  std::pow( 1/( 137.04 ) , 2 ) * std::pow( 0.38616e-12 , 2 ) * P * P * ( P + (1/P) - 1 + std::pow( std::cos(theta) , 2 ) ) / 2;

}

Double_t kleinNishna2(Double_t *x, Double_t *par)
{
	Double_t energy = par[0];             //Energy of photon (MeV)
	Double_t theta = x[0];

	Double_t P = 1 / ( 1 + ( energy/0.511 ) * ( 1 - std::cos(theta) ) );

	return  par[1] * std::pow( 1/( 137.04 ) , 2 ) * std::pow( 0.38616e-12 , 2 ) * P * P * ( P + (1/P) - 1 + std::pow( std::cos(theta) , 2 ) ) / 2;

}

inline double siMassAttenuation(double photonEnergy){return 0.1462*std::pow(photonEnergy, -0.469);} //Ouputs mass attenuation coefficent for si using photon energy. (Valid only between 0.08-6 MeV)
inline double h2oMassAttenuation(double photonEnergy){return 0.0666*std::pow(photonEnergy, -0.446);} //Ouputs mass attenuation coefficent for si using photon energy. (Valid only between 0.08-6 MeV)
inline double peMassAttenuation(double photonEnergy){return 0.0627*std::pow(photonEnergy, -0.446);} //Ouputs mass attenuation coefficent for si using photon energy. (Valid only between 0.08-6 MeV)
inline double alMassAttenuation(double photonEnergy){return 0.1624*std::pow(photonEnergy, -0.462);}
inline double airMassAttenuation(double photonEnergy){return 7e-5*std::pow(photonEnergy, -0.444);}

void photonGenerator(std::vector<double> &photonData, TF1 *bremFunc, double numOfEvents); //Generates a given number of photons from the spectrum.
void electronGenerator(TRandom *r, TFile &file, std::vector<double> &photonData, Event &tempEventData, TTree &eventTree, std::string material, double materialThickness, TF1 *kleinNishnaFunc); //Generates electrons based on the input of photon energies, material and material thickness. Outputs TTree of electron data.
void hitCounter(TTree &tree, double heightFromSensor, double thickness, std::string material, double &cumilativeHits, bool BU , TF1 * landau);
void materialBlock(TFile &file, std::vector<double> &photonData, std::string material, double thickness, TF1 *kleinNishnaFunc, TRandom3 *r);
void setLandauPar(double mpv, double sigma, TF1 *landau);
double electronEnergy(double photonEnergy, double photonAngle);
double electronAngle(double photonEnergy, double photonAngle);
double energyLosses(double electronKEnergy, std::string material, double materialThickness, TF1 * landau); //Material thickness in (cm)
double delta(BetheValues constants, double eta);
double F(double tau, double beta);
double f(double Z);
double electronRadius(double electronEnergy, double electronAngle);
double segmentLength(double radius, double startHeight, double endHeight);
double sigma(double electronKEnergy, std::string material, double materialThickness);
BetheValues SetParameters(std::string material);
TH1F * createAndFillHist(std::string name, std::vector<double> input, double bins, double min, double max);
TF1 * initialiseBremsstrahlung(double Z, double E0);
TRandom3 * initialiseRandomSeed();
const char* toChar(std::string input);
std::string toChar(double input);


int main(int argc, const char** argv)
{


	bool BU=false;
	double airGap = 0;
	double buildUp= 0;
	double runs(0);
	int photons(0);
	double alThickness(0);
	int bu(1);

	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >>  photons >> bu;//write the inputs

	if (bu==1)
	{
		BU=true;
	}
	else if (bu==0)
	{
		BU=false;
	}

	TF1 * bremFunc = initialiseBremsstrahlung(74,6);
	TRandom3 *r = initialiseRandomSeed();

	double cumilativeHits(0);
	// TH1F *electronHits = new TH1F(mainName, mainName, 20, 500, 1500);
	// TTree histTree("electronHits", "electronHits"); //Generate a TTree.
	// histTree.Branch("hits", &cumilativeHits, "Hits/D");//Generating branch.
	TF1 *kleinNishnaFunc2 = new TF1("kleinNishna2", kleinNishna2, 0, 1.55, 2);
	TF1 *kleinNishnaFunc= new TF1("kleinNishna", kleinNishna, 0, 3.14159 , 1);
	kleinNishnaFunc2->SetParameter(0, 1.4);
	kleinNishnaFunc2->SetParameter(1, 1e30);
	TF1 *landau = new TF1("landau", "TMath::Landau(x,[0],[1])", -0, 200);



	//kleinNishnaFunc2->Write();




	for(double airGap = 0.001; airGap<0.006; airGap+=0.001)
	{
		for (double alThickness = 0.02; alThickness < 0.11; alThickness+=0.01)
		{



			TGraphErrors *dataGraph = new TGraphErrors(bu+4);
			dataGraph->SetName("Data");
			if (bu==0)
			{
				dataGraph->SetPoint(1,0.009, 0.265);
				dataGraph->SetPoint(0,0.003, 0.22);
				dataGraph->SetPoint(2,0.012, 0.2556);
				dataGraph->SetPoint(3,0.015, 0.2566);
				dataGraph->SetPointError(0,0, 3.6e-3);
				dataGraph->SetPointError(1,0, 4.086e-3);
				dataGraph->SetPointError(2,0, 4.38e-3);
				dataGraph->SetPointError(3,0, 4.38e-3);
			}
			else if (bu==1)
			{

				dataGraph->SetPoint(0,0.003, 0.206);
				dataGraph->SetPoint(1,0.006, 0.202);
				dataGraph->SetPoint(2,0.009, 0.199);
				dataGraph->SetPoint(3,0.012, 0.198);
				dataGraph->SetPoint(4,0.015, 0.195);
				dataGraph->SetPointError(0,0, 3.02e-3);
				dataGraph->SetPointError(1,0, 3.08e-3);
				dataGraph->SetPointError(2,0, 3.139e-3);
				dataGraph->SetPointError(3,0, 3.139e-3);
				dataGraph->SetPointError(3,0, 3.0715e-3);
			}
			TGraphErrors simGraph(5);

			int n(0);

			for (double  buildUp = 0.003; buildUp <= 0.015; buildUp+=0.003)
			{







				std::string name;
				if(BU)
				{

					name = "BUnew/electronData/" + toChar(airGap) + toChar("_") + toChar(buildUp) + toChar("_") + toChar(alThickness) + "_" + toChar(photons) + "_" + toChar("_loops_BU_electron_data.root");
				}
				else if(!BU)
				{

					name ="noBUnew/electronData/" + toChar(airGap) + toChar("_") + toChar(buildUp) + toChar("_") + toChar(alThickness) + "_" + toChar(photons) + "_" + toChar("_loops_electron_data.root");
				}
				std::cout<<name<<std::endl;
				const char* outputName = toChar(name);

				TFile file(outputName, "RECREATE");
				std::vector<double> photonData;
				double height = buildUp+airGap;

				photonGenerator(photonData, bremFunc, photons);
				TH1F *initialSpectrum = createAndFillHist("initialSpectrum", photonData ,60, 0, 6);
				materialBlock(file, photonData, "al", alThickness*100, kleinNishnaFunc,r);

				cumilativeHits=0;

				for (double i = 0 ; i < ((height)/0.001)+1 ; ++i)
				{

					double tempHeight = height - (i*0.001);
					std::string s = toChar(tempHeight) + "(m)";
					const char *treeName = s.c_str();
					Event tempEventData; //Instantiate temp event class to store data.
					TTree tree(treeName, treeName); //Generate a TTree.
					tree.Branch("EventData", &tempEventData, "incPhotonEnergy/D:electronEnergy/D:photonAngle/D:electronAngle/D:electronKEnergy/D:electronRadius/D");//Generating branch.

					if (i<(int)(buildUp/0.001))
					{
						electronGenerator(r, file,  photonData, tempEventData, tree, "h2o", 0.1, kleinNishnaFunc);
						hitCounter(tree, height-(i*0.001), (buildUp) - (i*0.001), "bulk",cumilativeHits, BU, landau);
					}
					else if (i==(int)(height/0.001))
					{
						std::cout<<"Silico"<<std::endl;
						electronGenerator(r, file,  photonData, tempEventData, tree, "si", 0.0015, kleinNishnaFunc);
						for (int i = 1; i <= tree.GetEntries(); ++i)
						{	
							double tempEnergy = tree.GetLeaf("EventData","electronEnergy")->GetValue(0);
							cumilativeHits +=energyLosses(tempEnergy, "si", 0.0015 , landau);
						}

					}
					else if(i==(int)((height/0.001)-1) && BU)
					{
						std::cout<<"PE layer"<<std::endl;
						std::cout<<height-(i*0.001)<<std::endl;
						electronGenerator(r, file,  photonData, tempEventData, tree, "pe", 0.1, kleinNishnaFunc);
						hitCounter(tree, height-(i*0.001), 0.001, "pe",cumilativeHits, BU, landau);
					}
					else
					{
						electronGenerator(r, file,  photonData, tempEventData, tree, "air", 0.1, kleinNishnaFunc);
						hitCounter(tree, height-(i*0.001), (height-buildUp) - (i*0.001), "air",cumilativeHits, BU, landau);
					}
					tree.Write();
					std::cout<<((height/0.001)-1)<<" "<<i<<std::endl;
					std::cout<<"Layer: "<<tempHeight<<" Hits: "<<cumilativeHits<<" bu "<<BU<<std::endl;
				}




				std::cout<<" Total Hits: "<<cumilativeHits<<" bu "<<buildUp<<" n"<<n<<std::endl;
				if (bu==0)
				{
					if (n!=1)
					{
						simGraph.SetPoint(n,buildUp,cumilativeHits);
					}
					
				}
				else if (bu==1)
				{
					simGraph.SetPoint(n,buildUp,cumilativeHits);
				}


				n++;
				// histTree.Fill();
				TH1F *spectrum = createAndFillHist("spectrum", photonData, 60, 0, 6);
				TH1F *difference = new TH1F("diff","diff",60,0,6);
				for (int i = 1; i <= spectrum->GetNbinsX(); ++i)
				{
					double temp =  spectrum->GetBinContent(i) - initialSpectrum->GetBinContent(i) ;
					difference->SetBinContent(i,temp);
				}
				spectrum->Write();
				initialSpectrum->Write();
				difference->Write();
				difference->Draw("hist");

				file.Close();

			}
			std::string mainname, pdfname;

			if(bu==1)//File naming 
			{
				BU=true;
				airGap+=0.001;
				mainname = "BUnew/sensorData/" + toChar(airGap)  + "_" + toChar(alThickness) + "_" + toChar(photons) + "_loops_BU_hit_data.root";
				pdfname = "BUnew/sensorData/" + toChar(airGap)  + "_" + toChar(alThickness) + "_" + toChar(photons) + "_loops_BU_hit_data.pdf";
			}
			else if(bu==0)
			{	
				BU=false;

				mainname = "noBUnew/sensorData/" + toChar(airGap)  + "_" + toChar(alThickness) + "_" + toChar(photons) + "_loops_hit_data.root";
				pdfname = "noBUnew/sensorData/" + toChar(airGap)  + "_" + toChar(alThickness) + "_" + toChar(photons) + "_loops_hit_data.pdf";
			}
			const char * mainName = toChar(mainname);
			const char * pdfName = toChar(pdfname);
			std::cout<<mainName<<std::endl;
			TFile *outFile = new TFile(mainName, "RECREATE");
			TGraphErrors *normalisedSignal = new TGraphErrors(5);
			normalisedSignal->SetTitle("Simulation");

			double totalSum(0);
			if (bu==0)
			{
				for (int i = 0; i < 5; ++i)
				{
					double tempx(0), tempy(0);
					simGraph.GetPoint(i,tempx, tempy);
					if (n!=1)
					{
						totalSum+=tempy;
					}

				}

				for (int i = 0; i < 5; ++i)
				{
					if (n!=1)
					{
						double tempx(0), tempy(0);
						simGraph.GetPoint(i,tempx, tempy);
						normalisedSignal->SetPoint(i, tempx, tempy/totalSum );
						normalisedSignal->SetPointError(i,0,(tempy/totalSum)*0.05);
					}

				}
			}

			else
			{
				for (int i = 0; i < 6; ++i)
				{
					double tempx(0), tempy(0);
					simGraph.GetPoint(i,tempx, tempy);
					totalSum+=tempy;
				}

				for (int i = 0; i < 6; ++i)
				{
				
					double tempx(0), tempy(0);
					simGraph.GetPoint(i,tempx, tempy);
					normalisedSignal->SetPoint(i, tempx, tempy/totalSum );
					normalisedSignal->SetPointError(i,0,(tempy/totalSum)*0.05);
					
				}
			}


			dataGraph->SetMarkerStyle(5);
			dataGraph->SetTitle("Data");

			normalisedSignal->SetMarkerStyle(5);
			normalisedSignal->SetMarkerColor(kRed);
			dataGraph->Draw("AP");
			normalisedSignal->Draw("AP");
			TCanvas canvas("canvas", "canvas", 600, 400);
			canvas.SetGrid();
			canvas.SetTitle(mainName);
			TMultiGraph *multi = new TMultiGraph();
			multi->Add(dataGraph);
			multi->Add(normalisedSignal);
			multi->SetTitle(mainName);
			multi->Draw("AP");
			multi->GetYaxis()->SetTitle("Normalised Signal");
			multi->GetXaxis()->SetTitle("Build Up Thickness (m)");

			multi->Write();

			//canvas.BuildLegend();
			//canvas.Update();
			canvas.Write();
			canvas.SaveAs(pdfName);
			normalisedSignal->Write();
			simGraph.Write();
			outFile->Close();
		}
	}

	// TFile * mainoutput = new TFile(mainName, "RECREATE");



	// histTree.Write();
	// mainoutput->Close();
}

TF1 * initialiseBremsstrahlung(double Z, double E0)
{
	TF1 *bremFunc = new TF1("bremfunc", bremsstrahlung, 0.35,6,2); //Sets up TF1 for bremsstrahlung spectrum.
	bremFunc->SetParameter(0, Z); //Setting parameters
	bremFunc->SetParameter(1, E0);
	return bremFunc;
}
TRandom3 * initialiseRandomSeed()
{
	TRandom3 *r = new TRandom3(0);//Setting random number seeds
	gRandom = r;
	return r;
}

void hitCounter(TTree &tree, double heightFromSensor, double thickness, std::string material, double &cumilativeHits, bool BU , TF1 * landau)
{

	double hit(0);
	if (material=="bulk")
	{
		for (int i = 1; i <= tree.GetEntries(); ++i)
		{
			tree.GetEntry(i);
			double tempEnergy = tree.GetLeaf("EventData","electronEnergy")->GetValue(0); //Getting leaf values for electron KE and radius
			double tempRadius = tree.GetLeaf("EventData","electronRadius")->GetValue(0);
			double tempAngle = tree.GetLeaf("EventData","electronAngle")->GetValue(0);
			double tempKEnergy = tree.GetLeaf("EventData","electronKEnergy")->GetValue(0);
			double tempLoss(0);
			//std::cout<<"Temp angle "<<tempAngle<<std::endl;
			//Bulk
			if(!BU)
			{
				tempEnergy = tempEnergy - energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, 0, thickness) * 100 , landau); //Calculating loss through bulk
				//Calculating loss through sensor
				if ((2 * tempRadius) >  (heightFromSensor + 1.5e-5) && tempEnergy > 0.511)
				{
					tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
					if (tempEnergy-tempLoss > 0.511)
					{
						hit+=tempLoss;
					}
					else
					{
						hit+=(tempEnergy-0.511);
					}
					
					tempEnergy = tempEnergy - tempLoss;
					if (tempEnergy<0.511)
						{
						break;
						}

					double loops(0);
					while(loops==0)
					{
						tempRadius = electronRadius(tempEnergy, tempAngle);
						tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1e-4)) * 100 , landau);
						tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, 0, thickness) * 100, landau ); //Calculating loss through bulk
						if ((2 * tempRadius) >  (heightFromSensor + 1.5e-5) && tempEnergy > 0.511)
						{
							tempLoss = energyLosses(tempKEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
							if (tempEnergy-tempLoss > 0.511)
							{
								hit+=tempLoss;
							}
							else
							{
								hit+=(tempEnergy-0.511);
							}
					
							tempEnergy = tempEnergy - tempLoss;
							if (tempEnergy<0.511)
							{
							break;
							}
						}	
						else
						{
							loops=1;
						}
					}

				}
			}
			else
			{
				tempEnergy = tempEnergy - energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, 0, thickness) * 100 , landau); //Calculating loss through bulk
				tempEnergy = tempEnergy - energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, (heightFromSensor-0.001), heightFromSensor) * 100 , landau);

				if ((2 * tempRadius) >  (heightFromSensor + 1.5e-5) && tempEnergy > 0.511) //Test to see if radius is large enough and whether the electronwill have enough energy.
				{

					tempLoss = energyLosses(tempKEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
					if (tempEnergy-tempLoss > 0.511)
					{
						hit+=tempLoss;
					}
					else
					{
						hit+=(tempEnergy-0.511);
					}
					
					tempEnergy = tempEnergy - tempLoss;
					if (tempEnergy<0.511)
					{
					break;
					}

					double loops(0);
					while(loops==0)
					{
						tempRadius = electronRadius(tempEnergy, tempAngle);

						tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1e-4)) * 100 , landau);
						tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, (heightFromSensor-0.001), heightFromSensor) * 100 , landau);
						tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, 0, thickness) * 100 , landau);
						if ((2 * tempRadius) >  (heightFromSensor + 1.5e-5) && tempEnergy > 0.511)
						{
							
							tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
							if (tempEnergy-tempLoss > 0.511)
							{
								hit+=tempLoss;
							}
							else
							{
								hit+=(tempEnergy-0.511);
							}
							
							tempEnergy = tempEnergy - tempLoss;
							if (tempEnergy<0.511)
							{
							break;
							}
						}
						else
						{
							loops=1;
						}
					}
				}

			}
		}
	}

	else if(material=="air")
	{
		for (int i = 1; i <= tree.GetEntries(); ++i)
		{
			tree.GetEntry(i);
			double tempEnergy = tree.GetLeaf("EventData","electronEnergy")->GetValue(0); //Getting leaf values for electron KE and radius
			double tempRadius = tree.GetLeaf("EventData","electronRadius")->GetValue(0);
			double tempKEnergy = tree.GetLeaf("EventData","electronKEnergy")->GetValue(0);
			double tempAngle = tree.GetLeaf("EventData","electronAngle")->GetValue(0);
			double tempLoss(0);
			//Air - assume no losses through air.
			if(!BU)
			{

				if ((2 * tempRadius) > (heightFromSensor + 1.5e-5) && tempEnergy > 0.511) //Test to see if radius is large enough and whether the electronwill have enough energy.
				{

					tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
							if (tempEnergy-tempLoss > 0.511)
							{
								hit+=tempLoss;
							}
							else
							{
								hit+=(tempEnergy-0.511);
							}
							
					tempEnergy = tempEnergy - tempLoss;
					if (tempEnergy<0.511)
					{
					break;
					}
					double loops(0);
					while(loops==0)
					{
						
						tempRadius = electronRadius(tempEnergy, tempAngle);
						tempEnergy = tempEnergy - energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1e-4)) * 100 , landau);

						if ((2 * tempRadius) >  (heightFromSensor + 1.5e-5) && tempEnergy > 0.511)
						{
							tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
							if (tempEnergy-tempLoss > 0.511)
							{
								hit+=tempLoss;
							}
							else
							{
								hit+=(tempEnergy-0.511);
							}
							
							tempEnergy = tempEnergy - tempLoss;
							if (tempEnergy<0.511)
							{
							break;
							}
						}
						else
							{
								loops=1;
							}
					}
				} //Calculating loss through sensor

			}
			else
			{
				tempEnergy = tempEnergy - energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, (heightFromSensor-0.001), heightFromSensor) * 100 , landau);
				//Calculating loss through sensor

				if ((2*tempRadius) > (heightFromSensor + 1.5e-5) && tempEnergy > 0.511) //Test to see if radius is large enough and whether the electronwill have enough energy.
				{

					tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
					if (tempEnergy-tempLoss > 0.511)
						{
							hit+=tempLoss;
						}
					else
						{
							hit+=(tempEnergy-0.511);
						}
							
				tempEnergy = tempEnergy - tempLoss;
				if (tempEnergy<0.511)
				{
					break;
				}
					double loops(0);
					while(loops==0)
					{
						tempRadius = electronRadius(tempEnergy, tempAngle);
						tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1e-4)) * 100 , landau);
						tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, (heightFromSensor-0.001), heightFromSensor) * 100 , landau);

						if ((2 * tempRadius) >  (heightFromSensor + 1.5e-5) && tempEnergy > 0.511)
						{
							tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
							if (tempEnergy-tempLoss > 0.511)
							{
								hit+=tempLoss;
							}
							else
							{
								hit+=(tempEnergy-0.511);
							}
							
							tempEnergy = tempEnergy - tempLoss;
							if (tempEnergy<0.511)
							{
							break;
							}
						}
						else
						{
							loops=1;
						}
					}
				}
			}


		}
	}

	else if(material=="pe")
	{
		std::cout<<tree.GetEntries()<<std::endl;
		for (int i = 1; i <= tree.GetEntries(); ++i)
		{
			tree.GetEntry(i);

			double tempEnergy = tree.GetLeaf("EventData","electronEnergy")->GetValue(0); //Getting leaf values for electron KE and radius
			double tempRadius = tree.GetLeaf("EventData","electronRadius")->GetValue(0);
			double tempKEnergy = tree.GetLeaf("EventData","electronKEnergy")->GetValue(0);
			double tempAngle = tree.GetLeaf("EventData","electronAngle")->GetValue(0);
			double tempLoss(0);
			double tempHeightFromSensor(0);
			if ((heightFromSensor-0.001)<0)
			{
				tempHeightFromSensor = 0;
			}
			//std::cout<<"Energy "<<tempEnergy<<" Radius "<<tempRadius<<std::endl;
			//Air - assume no losses through air.
			tempEnergy = tempEnergy - energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, tempHeightFromSensor, heightFromSensor) * 100. , landau );

			//Calculating loss through sensor

			if ((2*tempRadius) > ( heightFromSensor + 1.5e-5) && tempEnergy > 0.511) //Test to see if radius is large enough and whether the electronwill have enough energy.
			{

				tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
					if (tempEnergy-tempLoss > 0.511)
					{
						hit+=tempLoss;
					}
					else
					{
						hit+=(tempEnergy-0.511);
					}
							
				tempEnergy = tempEnergy - tempLoss;
				if (tempEnergy<0.511)
				{
					break;
				}

				double loops(0);
				while(loops==0)
				{
					//std::cout<<"temp energy "<<tempEnergy<<" tempa angle "<<tempAngle<<std::endl;
					tempRadius = electronRadius(tempEnergy, tempAngle);
					//std::cout<<tempRadius<<std::endl;
					tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1e-4)) * 100. , landau); 
					tempEnergy = tempEnergy - 2 * energyLosses(tempEnergy, "h2o", segmentLength(tempRadius, (heightFromSensor-0.001), heightFromSensor) * 100. , landau);
					if ((2 * tempRadius) >  (heightFromSensor + 1.5e-5) && tempEnergy > 0.511)
					{
						tempLoss = energyLosses(tempEnergy, "si", segmentLength(tempRadius, heightFromSensor, (heightFromSensor + 1.5e-5)) * 100, landau );
							if (tempEnergy-tempLoss > 0.511)
							{
								hit+=tempLoss;
							}
							else
							{
								hit+=(tempEnergy-0.511);
							}
							
						tempEnergy = tempEnergy - tempLoss;
						if (tempEnergy<0.511)
					{
					break;
					}
					}
					else
					{
						loops=1;
					}
				}
			}
		}
	}


	std::cout<<"Hits: "<<hit<<std::endl;
	cumilativeHits+=hit;
}
void photonGenerator(std::vector<double> &photonData, TF1 *bremFunc, double numOfEvents)
{
	TStopwatch t;
	t.Start();
	for (int i = 0; i < numOfEvents; ++i)
	{

		double photonEnergy = bremFunc->GetRandom(); //Get random photon energy from given bremsstrahlung spectrum.
		photonData.push_back(photonEnergy);  //Add energy to photon vector.

	}
	t.Stop();
	t.Print();

}

void electronGenerator(TRandom *r, TFile &file, std::vector<double> &photonData, Event &tempEventData, TTree &eventTree, std::string material, double materialThickness, TF1 *kleinNishnaFunc)//Material thickness in (cm)
{

	// TF1 kleinNishnaFunc("kleinNishna", kleinNishna, 0, 3.14159, 1); //Sets up TF1 for Klein-Nishna distribution of scattering angles.

	int dataSize = photonData.size();



	int check(0);
	TStopwatch t;

	t.Start();



	for (int i = 0; i < photonData.size(); ++i) //Loop over all photons in the vector
	{


		//std::cout<<"\r"<<i<<"/"<<dataSize<<std::flush;
		double energy = photonData.at(i); //Get photon energy from vector


		std::map<std::string, double> materialMap = { //Declaring map with the various attenuation coefficients.
				{"si",siMassAttenuation(energy)*materialThickness},
				{"al" ,alMassAttenuation(energy)*materialThickness},
				{"pe", peMassAttenuation(energy)*materialThickness},
				{"h2o",h2oMassAttenuation(energy)*materialThickness}
		};





		if ( r->Uniform(1)<materialMap[material] && energy!=-1 ) //Generates random number to decide whether the photon interacts and checks the photon has not previously.
		{
			check++;
			if (check==1)
			{
				//std::cout<<energy<<std::endl;
			}

			kleinNishnaFunc->SetParameter(0,energy);
			kleinNishnaFunc->Update();
			double tempPhotonAngle = kleinNishnaFunc->GetRandom();
			//std::cout<<kleinNishnaFunc.GetParameter(0)<<std::endl;

			tempEventData.electronKEnergy = electronEnergy(energy, tempPhotonAngle);
			tempEventData.photonAngle = tempPhotonAngle;
			tempEventData.electronEnergy = tempEventData.electronKEnergy + 0.511;
			tempEventData.electronAngle = electronAngle(energy, tempPhotonAngle);
			tempEventData.incPhotonEnergy= energy;
			tempEventData.electronRadius = electronRadius(tempEventData.electronEnergy, tempEventData.electronAngle);

			//photonData.at(i)=-1;
			photonData.at(i) = energy - tempEventData.electronKEnergy; //Deleting photon from vector.
			eventTree.Fill(); //Fill TTree with values.


		}


	}
	t.Stop();
	//t.Print();

	//eventTree.Write(); //Write TTree to file.

}

double electronEnergy(double photonEnergy, double photonAngle) //Returns electron energy from photon energy and scattering angle.
{
	double gamma = photonEnergy/0.511; //Energy in (MeV)
	return ( photonEnergy * gamma * ( 1 - std::cos( photonAngle ) ) ) / ( 1 + gamma * ( 1 - std::cos( photonAngle )) );

}

double electronAngle(double photonEnergy, double photonAngle) //Returns the electron scattering angle.
{
	double gamma = photonEnergy/0.511; //Energy in (MeV)
	return std::atan2(1 , ( ( 1 + gamma ) * std::tan( photonAngle/2 ) ) );
}



double energyLosses(double electronEnergy, std::string material, double materialThickness, TF1 * landau)
{

	if (std::isnan(materialThickness))
	{
		//std::cout<<"Nan Material Thickness"<<std::endl;
		return 0;
	}
	BetheValues constants;
	constants.SetParameters(material);//Z, A, rho, I, a, X0, X1, m, C0
	double electronKEnergy = (electronEnergy - 0.511);
	double tau = electronKEnergy/0.511; //Put into units of m_ec^2.
	double gamma = electronEnergy/0.511; //Gamma factor.
	double beta = std::sqrt( 1 - (1/std::pow(gamma,2) ) );
	double v = beta*c;
	double eta = beta*gamma;



	double N = constants.rho * Na / constants.A;
	//Calculates the energy loss due to collisional processes, based on the Bethe-Bloch equation for electrons
	double collisionalLosses = 0.1535 * constants.rho * ( constants.Z / constants.A ) * ( 1 / std::pow(beta,2) ) * ( std::log( ( std::pow(tau,2) * ( tau + 2 ) ) / ( 2 * std::pow( ( constants.I / 0.511e6 ) , 2 )  ) ) + F(tau,beta) - delta(constants,eta) );
	double radiativeLosses = N * electronEnergy * 4 * std::pow(constants.Z,2) * std::pow(re,2) * alpha * (  std::log( ( 2 * electronEnergy  ) / 0.511 ) - (1/3) - f(constants.Z) );
	//Calculates energy loss due to radiative processes e.g brem.
	
	//std::cout<<"sig"<<sig<<std::endl;

	if (material=="si")
	{
		double mpv = (collisionalLosses) ;
		double sig = sigma(electronKEnergy, material,materialThickness);
		//td::cout<<"sig "<<sig<<" mpv "<<mpv<<"electron KE "<<electronKEnergy<<std::endl;
		setLandauPar(mpv, sig, landau);
		return landau->GetRandom()*materialThickness;
	}


	else
	{
		return(collisionalLosses+ radiativeLosses)*materialThickness;
	}


}



double delta(BetheValues constants, double eta) //Calculates the density correction factor for Bethe-Bloch.
{
	double X = std::log10(eta);

	if (X<constants.X0)
	{
		return 0;
	}

	else if(X>constants.X0&&X<constants.X1)
	{
		return 4.6052 * X + constants.C0 + constants.a * std::pow( (constants.X1-X), constants.m );
	}

	else
	{
		return 4.6052 * X + constants.C0;
	}
}


double F(double tau, double beta) //Calculates term in Bethe-Bloch for electrons calculation.
{
	return 1 - std::pow(beta,2) + ( ( std::pow(tau,2) / 8 ) - ( 2*tau+1 ) * std::log(2) ) / std::pow( tau + 1 , 2 );
}

double f(double Z) //Calculates term in radiative loss calculation.
{
	double a = Z/137;
	return std::pow(a,2) * ( std::pow(1+std::pow(a,2), -1 ) + 0.20206 - 0.0369 * std::pow(a,2) + 0.0083 * std::pow(a,4) - 0.002 * std::pow(a,6) );
}

double electronRadius(double electronEnergy, double electronAngle)
{
	double gamma = electronEnergy/0.511; //Gamma factor.
	double beta = std::sqrt( 1 - (1/std::pow(gamma,2) ) );
	double v = beta*c;

	double vPerp = v * std::cos(electronAngle); //Calculating v that is perpendicular to the B field.
	double betaPerp = vPerp/c;
	double gammaPerp = 1 / std::sqrt( 1 - (std::pow(betaPerp,2) ) );
	double momentum = (gammaPerp * vPerp * (0.511/3e8)) * 1e6; //Calculating the momentum in the perpendicular direction. Note (eV not MeV).

	double radius = momentum  / ( 3e8 * 1.5); //Radius in (m).


	return radius;

}


double segmentLength(double radius, double startHeight, double endHeight)
{

	if (startHeight<endHeight && endHeight<radius && startHeight<radius)
	{
		double phi = std::asin( ( radius - endHeight ) / radius );
		double theta = std::asin( ( radius  - startHeight ) / radius );
		double psi = theta - phi;

		return radius * psi;
	}

	else if (startHeight<endHeight && endHeight > radius && startHeight < radius)
	{
		if (endHeight > 2 * radius)
		{
			endHeight = 2*radius;
		}
		double phi = std::asin( (  endHeight - radius) / radius );
		double theta = std::asin( ( radius  - startHeight ) / radius );
		double psi = theta + phi;

		return radius * psi;
	}

	else if (startHeight<endHeight && endHeight > radius && startHeight > radius)
	{
		if (endHeight > 2 * radius)
		{
			endHeight = 2*radius;
		}
		double phi = std::asin( (  endHeight - radius) / radius );
		double theta = std::asin( ( startHeight - radius ) / radius );
		double psi = phi - theta;

		return radius * psi;
	}
	else if (startHeight > 2 * radius && endHeight > 2 * radius )
	{
		return 0;
	}

	else
	{
		std::cout<<"Start height must be less than end height!"<<std::endl;
		std::cout<<startHeight<<" "<<endHeight<<" "<<radius<<std::endl;
		return 0;
	}
}

void materialBlock(TFile &file, std::vector<double> &photonData, std::string material, double thickness, TF1 *kleinNishnaFunc, TRandom3 *r)
{
	for(int i=0; i<photonData.size(); i++)
	{
		if ( r->Uniform(1)<alMassAttenuation(photonData.at(i))*(thickness) && photonData.at(i)!=-1 ) //Generates random number to decide whether the photon interacts and checks the photon has not previously.
		{
			photonData.at(i)=-1;
		}
	}

}

TH1F * createAndFillHist(std::string name, std::vector<double> input, double bins, double min, double max)
{
	const char * name_c = name.c_str();
	TH1F *temp = new TH1F(name_c, name_c, bins, min, max);
	for (std::vector<double>::iterator i = input.begin(); i != input.end(); ++i)
	{
		temp->Fill(*i);
	}

	return temp;
}



const char* toChar(std::string input)
{
	return input.c_str();
}
std::string toChar(double input)
{
	std::stringstream ss;
	ss << std::fixed << std::setprecision(3) << input;

	return ss.str();

}
void setLandauPar(double mpv, double sigma, TF1 *landau)
{
	landau->SetParameter(0, mpv);
	landau->SetParameter(1, sigma);
	landau->Update();
}

double sigma(double electronKEnergy, std::string material, double materialThickness)
{
	BetheValues constants;
	constants.SetParameters(material);
	double electronEnergy = (electronKEnergy + 0.511);
	double gamma = electronEnergy/0.511; //Gamma factor.
	double beta = std::sqrt( 1 - (1/std::pow(gamma,2) ) );
	double val = (4 * 0.1535)  * (constants.Z/constants.A) * ((materialThickness/constants.rho)/(beta*beta));
	if (std::isnan(val))
	{
		std::cout<<"Gamma = "<<gamma<<" Beta = "<<beta<<" Energy = "<<electronKEnergy<<" rho "<<materialThickness<<" Z "<<constants.Z/constants.A<<std::endl;
	}
	return val;
}