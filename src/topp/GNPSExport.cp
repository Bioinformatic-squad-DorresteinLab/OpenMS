// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Dorrestein Lab - University of California San Diego - https://dorresteinlab.ucsd.edu/$
// $Authors: Abinesh Sarvepalli and Louis Felix Nothias$
// $Contributors: Fabian Aicheler and Oliver Alka from Oliver Kohlbacher's group at Tubingen University$
// --------------------------------------------------------------------------

//----------------------------------------------------------
// Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_GNPSExport GNPSExport
  @brief Process and export MS/MS data (.MGF format) from a consensusXML file.

This tool was developed for the Feature Based Molecular Networking (FBMN) workflow on GNPS (https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash2.jsp)
See the FBMN workflow documentation here (https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)

In brief, after running an OpenMS "metabolomics" pipeline, the consensusXML file is used to export files needed for FBMN on GNPS.
These two files are:
	- The MS/MS spectral data file (.MGF format) which is generated  with the GNPSExport util.
	- The feature quantification table (.TXT format) which is generated with the TextExport util.
	
Requirements: 
	- The IDMapper has to be ran on the featureXML files, in order to associate MS/MS scan(s) (peptide identification) 
	with each features. These peptide identifications are used by the GNPSExport.
	- The FileFilter has to be ran on the consensusXML file, prior to the GNPSExport, in order to remove consensusElements 
	without MS/MS scans (peptide identification).

Parameters: 
	- Cosine Score Treshold @Abi please describe what is is doing EXACTLY 
	- Binning @Abi please describe what is is doing EXACTLY
	
	@Abi: are we using a max number of peptide annotations ? If yes how is it define currently ?

Options for the GNPSExport spectral processing are:
	- Most intense: the GNPSExport will output the MS/MS scan with the highest precursor ion intensity 
	as a representative MS/MS scan per consensusElement in the .MGF file.
	- Merge: the GNPSExport will first merge all the MS/MS scans for a consensusElement, 
	using the user-specified parameters (cosine score threshold, binning width), and output 
	the merged MS/MS scan as as a representative MS/MS scan per consensusElement in the .MGF file.
	- All MS/MS: the GNPSExport will output all the MS/MS scan(s) for consensusElements in the .MGF file.

A representative OpenMS-GNPS workflow has the following steps:
  1. Input mzML files
  2. Run the FeatureFinderMetabo tool on the mzML files.
  3. Run the IDMapper tool on the featureXML and mzML files.
  4. Run the MapAlignerPoseClustering tool on the featureXML files.
  5. Run the MetaboliteAdductDecharger on the featureXML files.
  6. Run the FeatureLinkerUnlabeledKD tool or FeatureLinkerUnlabeledQT, on the featureXML files and output a consensusXML file.
  8. Run the FileFilter on the consensusXML file to keep only consensusElements with at least MS/MS scan (peptide identification).  
  9. Run the GNPSExport on the "filtered consensusXML file" to export an .MGF file.
  10. Run the TextExport on the "filtered consensusXML file" to export an .TXT file.
  11. Upload your files to GNPS and run the Feature-Based Molecular Networking workflow. Instructions are here:
https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/

The GitHub for that ProteoSAFe workflow and an OpenMS python wrappers is available here:
https://github.com/Bioinformatic-squad-DorresteinLab/openms-gnps-workflow

An online version of the OpenMS-GNPS pipeline for FBMN running on CCMS server (http://proteomics.ucsd.edu/) is available on GNPS:
https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-OpenMS

GNPS (Global Natural Products Social Molecular Networking, https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash2.jsp)
is an open-access knowledge base for community-wide organisation and sharing of raw, processed
or identified tandem mass (MS/MS) spectrometry data. 
The GNPS web-platform makes possible to perform spectral library search against public MS/MS spectral libraries, 
as well as to perform various data analysis such as MS/MS molecular networking, Network Annotation Propagation 
Network Annotation Propagation (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006089)
and the DEREPLICATOR (https://www.nature.com/articles/nchembio.2219)
The GNPS paper is available here (https://www.nature.com/articles/nbt.3597)

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iostream>
#include <fstream>

using namespace OpenMS;
using namespace std;

class TOPPGNPSExport : public TOPPBase
{
public:
  TOPPGNPSExport() :
  TOPPBase("GNPSExport", "Tool to export consensus features into MGF format", false) {}

private:
  double DEF_COSINE_SIMILARITY = 0.8;
  double DEF_MERGE_BIN_SIZE = BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES;
  // double DEF_PRECURSOR_MZ_TOLERANCE = 0.0001;
  // double DEF_PRECURSOR_RT_TOLERANCE = 5;

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in_cm", "<file>", "", "input file containing consensus elements with \'peptide\' annotations");
    setValidFormats_("in_cm", ListUtils::create<String>("consensusXML"));

    registerInputFileList_("in_mzml", "<files>", ListUtils::create<String>(""), "original mzml files containing ms/ms spectrum information");
    setValidFormats_("in_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Output MGF file");
    setValidFormats_("out", ListUtils::create<String>("mgf"));

    registerStringOption_("output_type", "<choice>", "full_spectra", "specificity of mgf output information", false);
    setValidStrings_("output_type", ListUtils::create<String>("full_spectra,merged_spectra,most_intense"));

    registerTOPPSubsection_("merged_spectra", "Options for exporting mgf file with merged spectra per feature");
    registerDoubleOption_("merged_spectra:cos_similarity", "<num>", DEF_COSINE_SIMILARITY, "Cosine similarity threshold for merged_spectra output", false);
    registerDoubleOption_("merged_spectra:ms2_bin_size", "<num>", DEF_MERGE_BIN_SIZE, "Bin size (Da) when merging ms2 scans", false);
  }

  // the main function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    ProgressLogger progress_logger;
    progress_logger.setLogType(log_type_);

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String consensus_file_path(getStringOption_("in_cm"));
    StringList mzml_file_paths(getStringList_("in_mzml"));
    String out(getStringOption_("out"));
    String output_type(getStringOption_("output_type"));
    double cos_sim_threshold(getDoubleOption_("merged_spectra:cos_similarity"));

    ofstream output_file(out);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    // ConsensusMap
    ConsensusXMLFile consensus_file;
    consensus_file.setLogType(log_type_);
    ConsensusMap consensus_map;
    consensus_file.load(consensus_file_path, consensus_map);

    // MSExperiment
    vector<MSExperiment> ms_maps;
    for (auto mzml_file_path : mzml_file_paths)
    {
      MzMLFile mzml_file;
      MSExperiment map;
      mzml_file.setLogType(log_type_);
      mzml_file.load(mzml_file_path, map);

      ms_maps.push_back(map);
    }


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    progress_logger.startProgress(0, consensus_map.size(), "parsing features and ms2 identifications...");
    // std::stringstream output_stream;
    Size feature_count = 1;
    for (Size i = 0; i < consensus_map.size(); ++i)
    {
      progress_logger.setProgress(i);
      // current feature
      const ConsensusFeature feature = consensus_map[i];

      // store "mz rt" information from each scan
      stringstream scans_output;
      scans_output << setprecision(2) << setfill('0') << fixed;

      // determining charge and most intense feature for header
      BaseFeature::ChargeType charge = feature.getCharge();
      for (ConsensusFeature::HandleSetType::const_iterator feature_iter = feature.begin();
        feature_iter != feature.end(); ++feature_iter)
      {
        if (feature_iter->getCharge() > charge)
        {
          charge = feature_iter->getCharge();
        }
      }

      // print spectra information (PeptideIdentification tags)
      vector<PeptideIdentification> peptide_identifications = feature.getPeptideIdentifications();

      cout << "\nfor feature " << i << " e_" << feature.getUniqueId() << endl;

      // clean peptide identifications outside mz rt tol ranges

      // vector of <<map index, spectrum index>, most intense ms2 scan>
      vector<pair<pair<double,PeptideIdentification>,pair<int,int>>> peptides;

      // determine if current feature has peptide annotations
      bool should_skip_feature;
      if (!(should_skip_feature = peptide_identifications.empty()))
      {
        for (Size peptide_index = 0; peptide_index < peptide_identifications.size(); peptide_index++)
        {
          auto peptide_identification = peptide_identifications[peptide_index];

          // append spectra information to scans_output
          int map_index = -1, spectrum_index = -1;
          if (peptide_identification.metaValueExists("spectrum_index"))
          {
            spectrum_index = peptide_identification.getMetaValue("spectrum_index");
          }
          if (peptide_identification.metaValueExists("map_index"))
          {
            map_index = peptide_identification.getMetaValue("map_index");
          }

          if (map_index != -1 && spectrum_index != -1)
          {
            // TEMP: log debug map index and spectrum index values once they are found
            cout << "\tmap index\t" << map_index << "\tspectrum index\t" << spectrum_index << endl;

            // retrieve spectrum for current peptide annotation
            auto ms2_scan = ms_maps[map_index][spectrum_index];
            ms2_scan.sortByIntensity(true);

            if (ms2_scan.getMSLevel() == 2 && !ms2_scan.empty())
            {
              double similarity_index = 0;
              // double similarity_index = abs(feature.getMZ() - peptide_identification.getMZ()) + abs(feature.getRT() - peptide_identification.getRT());
              for(auto ms2_iter = ms2_scan.begin(); ms2_iter != ms2_scan.end(); ms2_iter++)
              {
                similarity_index += ms2_iter->getIntensity();
              }
              auto first_pair = pair<double,PeptideIdentification>(similarity_index, peptide_identification);
              auto second_pair = pair<int,int>(map_index, spectrum_index);

              peptides.push_back(pair<pair<double,PeptideIdentification>,pair<int,int>>(first_pair,second_pair));
            }
            // for debug purposes
            else
            {
              // should_skip_feature = true;

              if (ms2_scan.empty())
              {
                cout << "\t\t-ms2 scan is empty\t" << peptide_identification.getMetaValue("spectrum_reference") << endl;
              }
            }
          }
        }
      }
      else
      {
        cout << "empty peptide identification list" << endl;
      }

      // peptides list of < <similarity_index, PeptideIdentification>, <map_index, feature_index> >
      // with the remaining peptides left within mz/rt tol of most intense
      if (!should_skip_feature && !peptides.empty())
      {
        // prepare peptides for output with smallest difference in mz
        sort (peptides.begin(), peptides.end(), [](const pair<pair<double,PeptideIdentification>,pair<int,int>> &a,
          const pair<pair<double,PeptideIdentification>,pair<int,int>> &b)
        {
          return a.first.first > b.first.first;
        });

        // tmp stream for current feature
        output_file << setprecision(4) << fixed;

        // full spectra
        if (output_type == "full_spectra")
        {
          for (auto peptide : peptides)
          {
            output_file << "BEGIN IONS" << endl;

            output_file << "SCANS=" << to_string(feature_count) << endl;

            string filename = mzml_file_paths[peptide.second.first];
            Size parse_index = filename.rfind("/") + 1;
            filename = filename.substr(parse_index);
            output_file << "FEATURE_ID=e_" << feature.getUniqueId() << endl;

            output_file << "MSLEVEL=2" << endl;
            output_file << "CHARGE=" << to_string(charge == 0 ? 1 : charge) << "+" << endl;
            output_file << "PEPMASS=" << feature.getMZ() << endl;
            output_file << "FILE_INDEX=" << peptide.second.second << endl;
            output_file << "RTINSECONDS=" << peptide.first.second.getRT() << endl;

            auto ms2_scan = ms_maps[peptide.second.first][peptide.second.second];
            // sort spectra
            sort (ms2_scan.begin(), ms2_scan.end(), [](const Peak1D& a, const Peak1D& b)
            {
              return a.getMZ() < b.getMZ();
            });

            for (Size l = 0; l < ms2_scan.size(); l++)
            {
              output_file << ms2_scan[l].getMZ() << "\t" << (int) ms2_scan[l].getIntensity() << endl;
            }

            output_file << "END IONS" << endl << endl;
          }
          feature_count++;
        }
        // merged spectra
        else if(output_type == "merged_spectra")
        {
          // map mz to intensity
          map<double,double> ms2_block;

          // // MapType exp;

          const BinnedSpectrum binned_highest_int(ms_maps[peptides[0].second.first][peptides[0].second.second], BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

          // MERGE: merge all peptide annotation spectra (mz and intensity)
          vector<pair<double,double>> mz_intensity_all;
          Size peptide_count = 0;
          for(auto peptide : peptides)
          {
              int map_index = peptide.second.first;
              int spectra_index = peptide.second.second;

              auto spectrum = ms_maps[map_index][spectra_index];
              const BinnedSpectrum binned_spectrum(spectrum, BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

              BinnedSpectralContrastAngle bsca;
              double cosine_sim = bsca(binned_highest_int, binned_spectrum);
              cout << "\tsimilarity of peptide " << peptide_count++ << " = " << cosine_sim << endl;
              if(cosine_sim < cos_sim_threshold)
              {
                cout << "\t" << "merged out\t" << map_index << "\t" << spectra_index << endl;
                break;
              }

              for(auto spectrum_iter = spectrum.begin(); spectrum_iter != spectrum.end(); spectrum_iter++)
              {
                mz_intensity_all.push_back(pair<double,double>(spectrum_iter->getMZ(), spectrum_iter->getIntensity()));
              }
          }
          sort(mz_intensity_all.begin(), mz_intensity_all.end(),[](pair<double,double>a, pair<double,double>b)
          {
            return a.first > b.first;
          });

          // generate new spectrum
          vector<double> mz_merged;
          vector<double> intensity_merged;
          double last_mz = numeric_limits<double>::min();
          double delta_mz(getDoubleOption_("merged_spectra:ms2_bin_size"));
          double sum_mz = 0;
          double sum_intensity = 0;
          Size count(0);
          for(auto it_mz = mz_intensity_all.begin(); it_mz != mz_intensity_all.end(); it_mz++)
          {
            if(abs(it_mz->first-last_mz)>delta_mz && count>0)
            {
              mz_merged.push_back(sum_mz/count);
              intensity_merged.push_back(sum_intensity);

              sum_mz = 0;
              sum_intensity = 0;
              last_mz = it_mz->first;
              count = 0;
            }
            sum_mz += it_mz -> first;
            sum_intensity += it_mz->second;
            count++;
          }
          if(count > 0)
          {
            mz_merged.push_back(sum_mz/count);
            intensity_merged.push_back(sum_intensity);
          }

          if(mz_merged.size() < mz_intensity_all.size())
            cout << "\tmz_merged.size() = " << mz_merged.size() << endl;

          // zip mz and intensity

          for(Size ms2_block_index = 0; ms2_block_index < mz_merged.size(); ms2_block_index++)
          {
            ms2_block[mz_merged[ms2_block_index]] = intensity_merged[ms2_block_index];
          }

          // print
          output_file << "BEGIN IONS" << endl;

          output_file << "SCANS=" << feature_count++ << endl;
          output_file << "FEATURE_ID=e_" << feature.getUniqueId() << endl;

          output_file << "MSLEVEL=2" << endl;
          output_file << "CHARGE=" << std::to_string(charge == 0 ? 1 : charge) << "+" << endl;
          output_file << "PEPMASS=" << feature.getMZ() << endl;
          output_file << "FILE_INDEX=" << peptides[0].second.second << endl;
          output_file << "RTINSECONDS=" << peptides[0].first.second.getRT() << endl;

          for (auto ms2_iter = ms2_block.begin(); ms2_iter != ms2_block.end(); ++ms2_iter)
          {
            output_file << ms2_iter->first << "\t" << (int) ms2_iter->second << endl;
          }
          output_file << "END IONS" << endl << endl;
        }
        // most intense ms2 block
        else if(output_type == "most_intense")
        {
          // print
          output_file << "BEGIN IONS" << endl;

          output_file << "SCANS=" << feature_count++ << endl;
          output_file << "FEATURE_ID=e_" << feature.getUniqueId() << endl;

          output_file << "MSLEVEL=2" << endl;
          output_file << "CHARGE=" << std::to_string(charge == 0 ? 1 : charge) << "+" << endl;
          output_file << "PEPMASS=" << feature.getMZ() << endl;
          output_file << "FILE_INDEX=" << peptides[0].second.second << endl;
          output_file << "RTINSECONDS=" << peptides[0].first.second.getRT() << endl;

          auto ms2_scan = ms_maps[peptides[0].second.first][peptides[0].second.second];
          for (auto ms2_iter = ms2_scan.begin(); ms2_iter != ms2_scan.end(); ++ms2_iter)
          {
            output_file << ms2_iter->getMZ() << "\t" << (int) ms2_iter->getIntensity() << endl;
          }
          output_file << "END IONS" << endl << endl;
        }

        // output feature information to general outputStream
        output_file << endl;
      }
      // should skip printing feature ms2 spectra block
      else
      {
        // print empty block
        output_file << "BEGIN IONS" << endl;

        output_file << "SCANS=" << feature_count++ << endl;
        output_file << "FEATURE_ID=e_" << feature.getUniqueId() << endl;
        output_file << "MSLEVEL=2" << endl;
        output_file << "CHARGE=" << std::to_string(charge == 0 ? 1 : charge) << "+" << endl;
        output_file << "PEPMASS=" << feature.getMZ() << endl;
        output_file << "RTINSECONDS=" << feature.getRT() << endl;

        output_file << "END IONS" << endl << endl << endl;
      }
    }
    progress_logger.endProgress();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    output_file.close();

    return EXECUTION_OK;
  }
};

// the actual main functioned needed to create an executable
int main (int argc, const char** argv)
{
  TOPPGNPSExport tool;
  return tool.main(argc, argv);
}
/// @endcond
