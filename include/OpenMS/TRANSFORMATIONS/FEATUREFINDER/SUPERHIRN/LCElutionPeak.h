// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#ifndef LC_elution_peak_H
#define LC_elution_peak_H

#include <vector>
#include <map>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ConsensusIsotopePattern.h>

namespace OpenMS
{

	typedef std::multimap<int, MSPeak> elution_peak;
	typedef std::vector<elution_peak> MZ_series;
	typedef std::vector<elution_peak>::iterator MZ_series_ITERATOR;
	typedef std::multimap<int, MSPeak>::iterator SIGNAL_iterator;

	class OPENMS_DLLAPI LCElutionPeak
	{

			////////////////////////////////////////////////
			// declaration of the private members:

		private:

			// isotopic pattern:
			ConsensusIsotopePattern* isotopePattern;
			int fNrIsotopes;
			double f_observed_Mass;
			double fIsotopMass;

		protected:

			double fMonoMass;
			double fVolume;
			int fCharge;
			int fScanNumberStart;
			int fScanNumberApex;
			int fScanNumberEnd;
			double fapex_intensity;
			double fRT;
			double fStartTR;
			double fEndTR;
			double fpeak_area;
			double fSignalToNoise;
			double fSNIntensityThreshold;
			MSPeak* APEX;

			std::string elutionPeakExtraInfo;

			// the raw signals assigned to this peak
			std::multimap<int, MSPeak> intens_signals;
			//multimap<int, MSPeak> raw_intens_signals;
			std::multimap<int, int> CHRG_MAP;

			////////////////////////////////////////////////
			// declaration of the public members:

		public:

			// parameters to debug a ceratain mass range
			static double DEBUG_MASS_START;
			static double DEBUG_MASS_END;

			// cut off, where everything small than this precentile of the
			// apex is discarded
			// static float intensity_apex_percentil_cutoff;

			// resolution of the retention time, for peak area copmuting:
			//  static float TR_RESOLUTION;

			// class destructor
			~LCElutionPeak();

			// class constructor
			LCElutionPeak();
			// class constructor
			LCElutionPeak(MZ_series_ITERATOR, double);
			// class copy constructor
			LCElutionPeak(const LCElutionPeak&);
			// constructor for the object feature:
			LCElutionPeak(const LCElutionPeak*);

			//////////////////////////////////////////////////////
			// Analyze the LC elution peak
			void analyzeLCElutionPeak()
			{

				if (get_nb_ms_peaks() > 1)
				{

					CHRG_MAP.clear();

					// determine the intensity background baseline based on S/N
					// value:
					setSNIntensityThreshold();

					// Compute a varietiy of parameters for the LC elution peak
					computeLCElutionPeakParameters();

					// define parameters such as chrg, score
					compute_CHRG();

					// create the consensus pattern:
					createConsensIsotopPattern();
				}
				else
				{
					defineLCElutionPeakParametersFromMSPeak();
				}
			}

			//////////////////////////////////////////////////////
			// determine the intensity background baseline based on S/N
			// value:
			void setSNIntensityThreshold();

			//////////////////////////////////////////////////////
			// Compute a varietiy of parameters for the LC elution peak
			void computeLCElutionPeakParameters();

			// removes background peaks and computes the total peak area:
			// void compute_LC_peak_area();
			// computes the area of between 2 peaks:
			double compute_delta_area(double, double, double, double);
			// define the apex into the elution profile::
			// void define_apex();
			// removes peaks which have lower intensity than x percentile
			// of the apex:
			void remove_background_peak();
			// compute the charge state of the LC peak
			void compute_CHRG();
			// compute the score of the LC peak:
			//void compute_SCORE_and_SN();
			// define all required peak parameters from a single MS peak:
			void defineLCElutionPeakParametersFromMSPeak();
			//////////////////////////////////////////////////////

			//////////////////////////////////////////////////////
			// print all monositopic peak cluster along the LC profile:
			void createConsensIsotopPattern();

			// print the elution profile from a peak:
			void print_profile(std::ofstream*);
			// find the closest existing mz peak in the elution profile:
			MSPeak* find_true_peak(float);
			// print the elution profile from a peak:
			void show_info();

			//////////////////////////////////////////////////
			// overload operators:
			LCElutionPeak& operator=(const LCElutionPeak&);
			LCElutionPeak& operator<=(const LCElutionPeak&);
			LCElutionPeak& operator>=(const LCElutionPeak&);
			LCElutionPeak& operator<(const LCElutionPeak&);
			LCElutionPeak& operator>(const LCElutionPeak&);

			////////////////////////////////////////////////////////////////////////////////
			// print all monositopic peak cluster along the LC profile:
			//void printIsotopClusters();
			// print the consensus isotope pattern:
			//void printConsensIsotopPattern();

			void setElutionPeakExtraInfo(std::string in)
			{
				elutionPeakExtraInfo = in;
			}

			std::string getElutionPeakExtraInfo()
			{
				return elutionPeakExtraInfo;
			}


			///////////////////////////////
			// start here all the get / set
			// function to access the
			// variables of the class

			////////////////////////////
			// access signal_intens map:
			SIGNAL_iterator get_signal_list_start()
			{
				return intens_signals.begin();
			}

			SIGNAL_iterator get_signal_list_end()
			{
				return intens_signals.end();
			}


			////////////////////////////
			// access the raw signal intens map:
			//SIGNAL_iterator get_raw_signal_list_start(){return raw_intens_signals.begin();};
			//SIGNAL_iterator get_raw_signal_list_end(){return raw_intens_signals.end();};

			// update the retention time by the current tmp_scan_apex:
			void set_apex_retention_time(double IN)
			{
				fRT = IN;
			}


			// to update the list of score and charge state:
			void update_CHRGMAP(MSPeak* IN)
			{
				std::multimap<int, int>::iterator T = CHRG_MAP.find(IN->get_charge_state());
				if (T == CHRG_MAP.end())
				{
					CHRG_MAP.insert(std::make_pair(IN->get_charge_state(), 1));
				}
				else
				{
					(*T).second++;
				}
			}

			///////////////////
			// get scan apex:
			int get_scan_apex()
			{
				return fScanNumberApex;
			}

			double get_apex_intensity()
			{
				return fapex_intensity;
			}

			double get_apex_retention_time()
			{
				return fRT;
			}

			double get_apex_MZ()
			{
				return get_MZ(get_scan_apex());
			}


			//////////////////
			// get an intensity of a ms_peak
			float get_intensity(int IN)
			{
				return (*(intens_signals.find(IN))).second.get_intensity();
			}

			// get the original M/Z of a ms_peak
			double get_MZ(int);

			/////////////////
			// get the total peak area:
			double get_total_peak_area()
			{
				return fpeak_area;
			}

			////////////////
			// get start / end scan:
			int get_start_scan()
			{
				return fScanNumberStart;
			}

			int get_end_scan()
			{
				return fScanNumberEnd;
			}

			void set_start_retention_time(double IN)
			{
				fStartTR = IN;
			}

			double get_start_retention_time()
			{
				return fStartTR;
			}

			void set_end_retention_time(double IN)
			{
				fEndTR = IN;
			}

			double get_end_retention_time()
			{
				return fEndTR;
			}


			////////////////
			// get number of peaks in the elution profile:
			int get_nb_ms_peaks()
			{
				return (int) intens_signals.size();
			}


			//////////////
			// access teh charge state of the LC elutino peak:
			int get_charge_state()
			{
				return fCharge;
			}


			////////////
			// get signal to noise ratio:
			double getSignalToNoise()
			{
				return fSignalToNoise;
			}

			double getSignalToNoiseBackground()
			{
				return fSNIntensityThreshold;
			}

	};

} // ns

#endif
