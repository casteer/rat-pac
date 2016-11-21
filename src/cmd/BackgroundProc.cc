///////////////////////////////////////////////////////////////////////////////
/// \file BackgroundProc.cc
/// \brief A Background Processor
///
///
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include <TFile.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <RAT/Processor.hh>
#include <RAT/Log.hh>
#include <RAT/DS/Root.hh>

using namespace RAT;

namespace RAT {

// Class declaration
// (if this gets big, you should move it to its own header file)
class BackgroundProc : public Processor {
public:

  BackgroundProc();
  virtual ~BackgroundProc();

  void CalculateBackgroundOutput();

  Double_t CalculateDoublesRate(Double_t n);

  Double_t CalculateSinglesRate(Double_t n);
  Double_t CalculateSinglesStatError(Double_t cut_events);

  Double_t CalculateDoublesPrecision();
  Double_t CalculateSinglesPrecision();

  Processor::Result ProcessSingleEvent(DS::Root *ds);

  virtual Processor::Result DSEvent(DS::Root *ds);
  // Override this method instead if you want to be called for
  // every triggered DAQ event, rather than every physics event
  // Only use one of these methods!
  //virtual Processor::Result Event(DS::Root *ds, DS::EV *ev);

  virtual void SetI(std::string param, int value);
  virtual void SetD(std::string param, double value);
  virtual void SetS(std::string param, std::string value);

  bool CheckPositive(std::string param, double value){
      if (value <= 0){
            throw ParamInvalid(param, "update interval must be > 0");
            return false;// Don't expect it to get here.
      } else {
          return true;
      }
  }


  //====================================
  // Default functions for pmts/scenarios etc
  //====================================
  void SetPMTSet(std::string set){
      if(set==std::string("inner")){
          SetNumberItems(4330);
      }else if(set==std::string("veto")){
          SetNumberItems(482);
      }
      SetActivity();
    }
  void SetPMTBackground(std::string kut, std::string pmt_type){
      SetSource(kut);
      if(pmt_type==std::string("ETL")){
          SetItemWeight(2.1);
          if(kut==std::string("K")){
              SetSourcePPB(60e3);
          }else if(kut==std::string("U")){
              SetSourcePPB(30);
          }else if(kut==std::string("Th")){
              SetSourcePPB(30);
          }
      }else if(pmt_type==std::string("Hamamatsu")){
          SetItemWeight(2.2);
          if(kut==std::string("K")){
              SetSourcePPB(36e3);
          }else if(kut==std::string("U")){
              SetSourcePPB(43);
          }else if(kut==std::string("Th")){
              SetSourcePPB(133);
          }
      }
      SetActivity();
  }

  void SetRockBackground(std::string kut,std::string location){

      SetSource(kut);
	
      // Boulby background from Alner et al, Astroparticle Physics 28 287-302 (2007)
      if(location==std::string("Boulby")){
	  std::cout << "BackgroundProc: Changing location to BOULBY" << std::endl;
	  std::cout << "BackgroundProc: Changing Source to " << kut << std::endl;
          if(kut==std::string("K")){
              SetSourcePPB(1100e3);
          }else if(kut==std::string("U")){
              SetSourcePPB(65);
          }else if(kut==std::string("Th")){
              SetSourcePPB(130);
          }
      }else if(location==std::string("Perry")){
          if(kut==std::string("K")){
              SetSourcePPB(1100e3);
          }else if(kut==std::string("U")){
              SetSourcePPB(65);
          }else if(kut==std::string("Th")){
              SetSourcePPB(130);
          }
      }
      SetActivity();
  }
  //====================================
  // Set functions for private variables
  //====================================

  // Calculate the total activity based on the weight per item, number of items, background ppb etc
  void SetActivity(){
      activity = number_items*item_weight*ppb*ppb2Bq;
      std::cout << "BackgroundProc: Activity = " << activity << " Bq " << std::endl;
  }

  // Set activity based on an input activity
  void SetActivity(Double_t activity_in){
      activity = activity_in;
      std::cout << "BackgroundProc: Activity = " << activity << " Bq " << std::endl;
  }

  void SetItemWeight(Double_t weight_in){
      if(CheckPositive(std::string("item_weight"),weight_in)){
          item_weight = weight_in;
          SetActivity();
          std::cout << "BackgroundProc: Weight per item = " << item_weight << " kg " << std::endl;
      }
  }

  void SetStopAtPrecision(Double_t value){
      if(CheckPositive(std::string("stop_at_precision"),value)){
          stop_at_precision = value;
          std::cout << "BackgroundProc: Stop at doubles rate precision of " << stop_at_precision << std::endl;
      }
  }

  void SetTimeWindow(Double_t window_in){
      if(CheckPositive(std::string("time_window"),window_in)){
        time_window = window_in;
        std::cout << "BackgroundProc: Time window = " << time_window << " s " << std::endl;
      }
  }

  void SetUpdate(Double_t update_in){
      if(CheckPositive(std::string("update"),update_in)){
        updateInterval = (int) update_in;
        std::cout << "BackgroundProc: Update frequency = " << updateInterval << std::endl;
      }
  }

  void SetUpdate(int update_in){
      if(CheckPositive(std::string("update"),(double) update_in)){
        updateInterval = (int) update_in;
        std::cout << "BackgroundProc: Update frequency = " << updateInterval << std::endl;
      }
  }

  void SetNumberItems(Double_t items_in){
      number_items = items_in;
      SetActivity();
      std::cout << "BackgroundProc: Number of items = " << number_items << std::endl;
  }

  void SetPECuts(Double_t min_pe,Double_t max_pe){

      if(CheckPositive(std::string("min_pe_cut"),min_pe)){
        min_pe_cut = min_pe;
      }

      if(CheckPositive(std::string("max_pe_cut"),max_pe)){
        max_pe_cut = max_pe;
      }

      std::cout << "BackgroundProc: Min. PE cut = " << min_pe_cut << std::endl;
      std::cout << "BackgroundProc: Max. PE cut = " << max_pe_cut << std::endl;

  }

  void SetNHitCuts(Double_t min_nhit,Double_t max_nhit){
      if(CheckPositive(std::string("min_nhit_cut"),min_nhit)){
        min_nhit_cut = min_nhit;
      }
      if(CheckPositive(std::string("max_nhit_cut"),max_nhit)){
        max_nhit_cut = max_nhit;
      }
      std::cout << "BackgroundProc: Min. Nhit cut = " << min_nhit_cut << std::endl;
      std::cout << "BackgroundProc: Max. Nhit cut = " << max_nhit_cut << std::endl;
  }

  // Set the strength of the source
  void SetSourcePPB(double ppb_in){
      if(CheckPositive(std::string("ppb"),ppb_in)){
        ppb=ppb_in;
      }
      SetActivity();
      std::cout << "BackgroundProc: Source ppb = " << ppb << std::endl;
  }

  void SetSource(std::string k_u_th){
      source_string=std::string(k_u_th);
      if(k_u_th==std::string("U")) ppb2Bq = 12.35e-3;
      if(k_u_th==std::string("Th")) ppb2Bq = 4.06e-3;
      if(k_u_th==std::string("K")) ppb2Bq = 31e-6;
      SetActivity();
      std::cout << "BackgroundProc: Source conversion factor = " << ppb2Bq << std::endl;
  }

  void Print(){

      std::cout << "BackgroundProc: Source Summary" << std::endl;
      std::cout << " - Source type: " << source_string << std::endl;
      std::cout << " - Source activity: " << activity << std::endl;
      std::cout << " - Source fraction (ppb): " << ppb << std::endl;
      std::cout << " - Coincidence time window: " << time_window << std::endl;

      std::cout << "BackgroundProc: Summary Stats" << std::endl;
      CalculateBackgroundOutput();

  }

  void OpenFile(std::string filename_in){

      std::cout << "BackgroundProc: Opening file " << filename_in << std::endl;

      filename = filename_in;
      // Needs more error checking here... zero-length strings etc.

      // First open the file, destroying the previous file
      std::cout << "BackgroundProc: Creating root file" << std::endl;
      f = new TFile(filename.c_str(),"Recreate");
      Log::Assert(f != NULL, "BackgroundProc: Unable to open ROOT ntuple output file.");

      // Associate the ntuple with the file
      data_out = new TNtuple("data","Ntuple for Watchman Gamma Background Studies","nhit:pe:mean_t:e_true:r_true:z_true:r_centroid:z_centroid:r_pathfit:z_pathfit:r_bonsai:z_bonsai:sub_ev");

  }


protected:
  TFile *f;
  // TH1F *hNumPE;

  Double_t activity;// The source activity
  Double_t time_window;// The time window
  
  Double_t min_pe_cut;// The pe cut minimum
  Double_t max_pe_cut;// The pe cut maximum
 
  Double_t min_nhit_cut;// The nhit cut
  Double_t max_nhit_cut;// The nhit cut

  Double_t rcut;// The fiducial radius cut
  Double_t zcut;// The fiducial z cut

  Double_t min_rcut;// The fiducial radius cut
  Double_t max_rcut;// The fiducial radius cut
  Double_t min_zcut;// The fiducial radius cut
  Double_t max_zcut;// The fiducial radius cut

  // For activity calculations
  Double_t number_items;// The number of items
  Double_t item_weight;// The part's weight in kg
  Double_t ppb;// The background source's ppb
  Double_t ppb2Bq;// Conversion from ppb to Bq/kg

  Int_t nevent;
  Double_t stop_at_precision;

  int dscount;        ///< Number of physics events
  int evcount;        ///< Number of triggered events
  int updateInterval; ///< Number of physics events per update line

  TNtuple* data_out;

  std::string source_string;
  std::string filename;

  const Double_t secs2day = 24.0*60.0*60.0;
};

}


namespace RAT {
     // Override user processor utility function
     Processor *construct_user_proc(std::string ){
    // If you have several different user processors, check the
    // userProcName parameter and allocate the proper one.  Only one
    // kind of user processor in this example, so we ignore userProcName.
    return new BackgroundProc;
  }


// Class definition
BackgroundProc::BackgroundProc() : Processor("user") {
  // Construct all your private objects and do any initialization
  // This gets called when your processor is created and
  // added to the event loop in the command file.

  std::cout << "======= BackgroundProc ======= " << std::endl;

  activity=1.0;// in Bq
  time_window=100e-6;// The time window in seconds
  
  min_pe_cut=5;// The pe cut minimum
  max_pe_cut=25;// The pe cut maximum

  min_nhit_cut=10;// The pe cut minimum
  max_nhit_cut=10;// The pe cut maximum

  number_items=1.0;// The number of items
  item_weight=1.0;// The part's weight in kg
  ppb=1.0;// The background source's ppb
  ppb2Bq=1.0;// Conversion from ppb to Bq/kg

  rcut=5400;
  zcut=5400;

  min_rcut=0.0;// The fiducial radius cut
  max_rcut=5400.0;// The fiducial radius cut

  min_zcut=0.0;// The fiducial z cut
  max_zcut=5400.0;// The fiducial z cut

  nevent=0;

  dscount=0;        ///< Number of physics events
  evcount=0;        ///< Number of triggered events

  updateInterval=100; ///< Number of physics events per update line

  //The simulation stops if it reaches this precision level (in events/day) for the doubles rate
  stop_at_precision=-1;// A -ve value means that it will not stop until the number of events requested is exhausted

  ppb = 1.0;// ppb of background material

  // Initialize the source type
  SetSource(std::string("U"));

  // In order to print out the activity
  // SetActivity(activity);

  filename=std::string("ntuple_output.root");


  // std::cout << "BackgroundProc: Final Set of Parameters... " << std::endl;
  // Print();

}

BackgroundProc::~BackgroundProc(){
  // Delete your private objects and print any summary information.
  // This gets called when the event loop is cleared by the user or at
  // the end of execution.

  std::cout << "Deleting ~BackgroundProc()" << std::endl;

  // Print summary of the set-up and stats
  Print();

  // Write to the file
  f->cd();
  data_out->Write();
  f->Close();

  delete f; // deletes hNumPE too -- ROOT does weird stuff like this
}

void BackgroundProc::SetI(std::string param, int value)
{
    if (param == std::string("update")){ SetUpdate(value); }
    else throw ParamUnknown(param);
}

void BackgroundProc::SetS(std::string param, std::string value)
{
    if (param == std::string("source")){ SetSource(value); }
    else if (param == std::string("filename")){ OpenFile(value);}
    else if (param == std::string("location")){ SetRockBackground(source_string,value); }
    else if (param == std::string("pmt_set")){ SetPMTSet(value); }
    else if (param == std::string("pmt_type")){ SetPMTBackground(source_string,value); }
    else throw ParamUnknown(param);
}

void BackgroundProc::SetD(std::string param, double value)
{

         if (param == std::string("update")){ SetUpdate(value); }
         else if (param == std::string("activity")) { SetActivity(value); }
         else if (param == std::string("min_fiducial_r")) { }
         else if (param == std::string("max_fiducial_r")) { }
         else if (param == std::string("min_fiducial_z")) { }
         else if (param == std::string("max_fiducial_z")) { }
         else if (param == std::string("min_pe")) { SetPECuts(value,max_pe_cut); }
         else if (param == std::string("max_pe")) { SetPECuts(min_pe_cut,value); }
         else if (param == std::string("min_nhit")) { SetNHitCuts(value,max_nhit_cut); }
         else if (param == std::string("max_nhit")) { SetNHitCuts(min_nhit_cut,value); }
         else if (param == std::string("ppb")) {SetSourcePPB(value);}
         else if (param == std::string("ppb2Bq")) { }
         else if (param == std::string("stop_at_precision")) { SetStopAtPrecision(value); }
         else if (param == std::string("number_items")) { SetNumberItems(value);}
         else if (param == std::string("item_weight")) { SetItemWeight(value);}
         else throw ParamUnknown(param);

}


// Calculates the doubles rate precision 
Double_t BackgroundProc::CalculateDoublesPrecision(){
        return 2.0*time_window*activity*activity*secs2day / (((double) dscount)*((double) dscount));
}

// Calculates the doubles rate
Double_t BackgroundProc::CalculateDoublesRate(Double_t cut_events){
  return 2.0*time_window*activity*activity*secs2day*((cut_events*cut_events)/ (((double) dscount)*((double) dscount)));
}

// Calculates the doubles rate
Double_t BackgroundProc::CalculateSinglesRate(Double_t cut_events){
  return ppb*ppb2Bq*activity*secs2day*(cut_events/((double) dscount));
}

// Calculates the doubles rate
Double_t BackgroundProc::CalculateSinglesPrecision(){
  return ppb*ppb2Bq*activity*secs2day*(1.0/((double) dscount));
}



// Calculates the statistical error from the ppb measurements and from the number of cut events
Double_t BackgroundProc::CalculateSinglesStatError(Double_t cut_events){

    if(cut_events==0) return 0.0;

    double err_ppb = 0.5;
    double err_events = 1.0/sqrt(cut_events);
    double err_tot = sqrt(err_ppb*err_ppb + err_events*err_events);

    return err_tot*ppb*ppb2Bq*activity*secs2day*(cut_events/((double) evcount));
}




// Calculates the doubles rate for a range of cuts
void BackgroundProc::CalculateBackgroundOutput(){

  // Write some of the numbers independent of pe and nhit for this update
  std::cout << "Singles rate precision (events/day): " << std::setprecision(6) << CalculateSinglesPrecision() << std::endl;
  std::cout << "Doubles rate precision (events/day): " << std::setprecision(6) << CalculateDoublesPrecision() << std::endl;

  if(!((max_pe_cut>=min_pe_cut)&&(max_nhit_cut>=min_nhit_cut))) return;

  // Find and write out the doubles rate
  std::cout << std::scientific ;
  std::cout << " - PE " ;
  std::cout << "\t \t Nhit " ;
  std::cout << "\t \t Ncut ";
  std::cout << "\t \t Singles Rate (events/day) " ;
  std::cout << "\t \t Doubles Rate (events/day) " << std::endl;
  // Loop over the pe and nhit cuts
  for(Double_t pe=min_pe_cut;pe<=max_pe_cut;pe+=2){
        for(Double_t nhit=min_nhit_cut;nhit<=max_nhit_cut;nhit++){

		// Create the cuts string char buffer
                char buff[200];
                // sprintf(buff, "(r_bonsai<%f)&&(abs(z_bonsai)<%f)&&(pe>%f)&&(nhit>%f)",rcut,zcut,pe,nhit);
                sprintf(buff, "(pe>%f)&&(nhit>%f)",pe,nhit);

		// Find how many pass the cut
                int nselected = data_out->Draw("pe",buff,"goff");
		
		// Find and write out the doubles rate 		
                std::cout << (int) pe << "\t \t";
                std::cout << (int) nhit << "\t \t";
                std::cout << nselected << std::scientific << "\t \t";
                std::cout << std::setprecision(6) << CalculateSinglesRate((Double_t) nselected) << "\t \t \t \t";
                std::cout << std::setprecision(6) << CalculateDoublesRate((Double_t) nselected) << std::endl;
     }
  }

  if(CalculateDoublesPrecision()<stop_at_precision){
      // Write the root file then die
      // Write to the file
      f->cd();
      data_out->Write();
      f->Close();
      Log::Die("BackgroundProc: Run stopping as it has reached the target precision");
  }

}

Processor::Result BackgroundProc::DSEvent(DS::Root *ds){

/*
 *     if (!f) { // Assume no file specified, and so we must open the default
      info << "BackgroundProc: No output file specified, opening " << filename << newline;
      OpenFile(filename);
    }
*/

  dscount++;
  evcount += ds->GetEVCount();
 
  // Continually update the ntuple for each event
  Processor::Result single_result = ProcessSingleEvent(ds);

  //... but only periodically produce stats on the doubles background rate etc... 
  if (dscount % updateInterval == 0){
	info << dformat("BackgroundProc: Event %d (%d triggered events)\n",dscount, evcount);
        CalculateBackgroundOutput();
  }
  
  return single_result;
}

Processor::Result BackgroundProc::ProcessSingleEvent(DS::Root *ds){
  // Process one event. This is called once per DETECTOR event, and not once
  // per PHYSICS event, i.e. once for every RAT::EV in the list.
  // hNumPE->Fill(ds->GetMC()->GetNumPE());

  // To save writing this out alot 
  RAT::DS::MC *mc = ds->GetMC();

  // Update the internal event counter
  nevent = mc->GetMCParticleCount();

  // Get the primary 
  RAT::DS::MCParticle *prim = mc->GetMCParticle(0);

  // - Find the (r,z) coordinates of the primary
  Double_t r_pFitFprim = sqrt(pow(prim->GetPosition().X(),2)+ pow(prim->GetPosition().Y(),2));
  Double_t z_pFitFprim = prim->GetPosition().Z();
  Double_t e_pFitFprim = prim->GetKE();

  Double_t r_pFitFB = -1.0;
  Double_t z_pFitFB = 0.0;

  Double_t r_pFitFC = -1.0;
  Double_t z_pFitFC = 0.0;

  Double_t r_pFitFP = -1.0;
  Double_t z_pFitFP = 0.0;
  
  Double_t totPE=0.0;

  Double_t pmtCount=0.0;

  // get the sub events
  Int_t subevents = ds->GetEVCount();
  for(int k = 0; k<subevents; k++) {
            
        // Get the event pointer 
  	RAT::DS::EV *ev = ds->GetEV(k);

        // Need all fitters to be present
        // bool fitters_all=(ev->ExistCentroid())&&(ev->ExistPathFit())&&(ev->ExistBonsaiFit());
        // Log::Assert(fitters_all==true, "BackgroundProc: Need all fitters first.");

        // The pmt count for this event
        pmtCount = Double_t(ev->GetPMTCount());

	// Centroid fit 
            RAT::DS::Centroid *fc = ev->GetCentroid();
            TVector3 pFitFC = fc->GetPosition();
            // - Find the (r,z) coordinates
            r_pFitFC = sqrt(pow(pFitFC.X(),2)+ pow(pFitFC.Y(),2));
            z_pFitFC = pFitFC.Z();

  	// Pathfit 
            RAT::DS::PathFit *pf = ev->GetPathFit();
            TVector3 pFitFP = pf->GetPosition();
            // - Find the (r,z) coordinates
            r_pFitFP = sqrt(pow(pFitFP.X(),2)+ pow(pFitFP.Y(),2));
            z_pFitFP = pFitFP.Z();

	// - Bonsai fit                     
            RAT::DS::BonsaiFit *pb = ev->GetBonsaiFit();
            TVector3 pFitFB = pb->GetPosition();
            // - Find the (r,z) coordinates
            r_pFitFB = sqrt(pow(pFitFB.X(),2)+ pow(pFitFB.Y(),2));
            z_pFitFB = pFitFB.Z();

	Double_t mean_t=0.0;
        for(int j = 0; j<ev->GetPMTCount();j++) {
		RAT::DS::PMT *pmt = ev->GetPMT(j);
                totPE+=pmt->GetCharge();// Charge in pC
		mean_t += pmt->GetTime()/((double) ev->GetPMTCount());

	}// Loop over the pmt count 

  // Fill the ntuple
  data_out->Fill(pmtCount,totPE,mean_t,e_pFitFprim,r_pFitFprim,z_pFitFprim,r_pFitFC,z_pFitFC,r_pFitFP,z_pFitFP,r_pFitFB,z_pFitFB,Double_t(k));

  }// End of the loop over sub events

  // Tell the event loop handler the result of processing.
  return Processor::OK;

}

}// Namespace RAT
