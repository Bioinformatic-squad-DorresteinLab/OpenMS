// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
Please cite our preprint: Nothias, L.F. et al, Feature-based Molecular Networking in the GNPS Analysis Environment
bioRxiv 812404 (2019) (https://www.biorxiv.org/content/10.1101/812404v1)
See the FBMN workflow documentation here (https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)
In brief, after running an OpenMS "metabolomics" pipeline, the GNPSExport TOPP tool can be used
on the consensusXML file and corresponding mzML files to generate the files needed for FBMN on GNPS.
These two files are:
	- The MS/MS spectral data file (.MGF format) which is generated  with the GNPSExport util.
	- The feature quantification table (.TXT format) which is generated with the TextExport util.
For each consensusElement in the consensusXML file, the GNPSExport produces one representative consensus
MS/MS spectrum (named peptide annotation in OpenMS jargon) outputed in the MS/MS spectral file (.MGF file).
Several mode for the generation of the consensus MS/MS spectrum are available and described below.
Note that these parameters are defined in the GNPSExport parameters file (.INI file, available with that link.
Representative command:
GNPSExport -ini iniFile-GNPSExport.ini -in_cm filefilter.consensusXML -in_mzml inputFile0.mzML inputFile1.mzML -out GNPSExport_output.mgf
The GNPSExport TOPP tool can be ran on a consensusXML file and the corresponding mzML files to generate a MS/MS spectral file (MGF format)
and corresponding feature quantification table (.TXT format) that contains the LC-MS peak area intensity.
Requirements:
	- The IDMapper has to be ran on the featureXML files, in order to associate MS2 scan(s) (peptide annotation) with each
	features. These peptide annotations are used by the GNPSExport.
	- The FileFilter has to be ran on the consensusXML file, prior to the GNPSExport, in order to remove consensusElements
	without MS2 scans (peptide annotation).
Parameters:
	- Cosine Score Treshold - defines Cosine Similarity Treshold for the pairwise cosine similarity between the MS/MS scan with the highest precursor intensity and the other MS/MS scans.
	- Binning - defines the Binning width of fragment ions during the merging of eligible MS/MS spectra.
Options for the GNPSExport spectral processing are:
	- [RECOMMENDED]: merged_spectra - For each consensusElement, the GNPSExport will merge all the
	eligible MS/MS scans into one representative consensus MS/MS spectrum. Eligible MS/MS scans have a
	pairwise cosine similarity with the MS/MS scan of highest precursor intensity above the Cosine Similarity Treshold.
	The fragment ions of merged MS/MS scans are binned in m/z (or Da) range defined by the Binning width parameter.
		- Cosine Similarity Treshold: merged_spectra:cos_similarity (float, default: 0.9) - Parameter that defines
		Cosine Similarity Treshold for the pairwise cosine similarity between the MS/MS scan with the highest
	precursor intensity and the other MS/MS scans.
		- Binning width: merged_spectra:ms2_binned_size (float, default: 0.02 Daltons) - Parameter that defines the
		Binning width of fragment ions during the merging of eligible MS/MS spectra.
	- Most intense: most_intense - For each consensusElement, the GNPSExport will output the most
	intense MS/MS scan (with the highest precursor ion intensity) as consensus MS/MS spectrum.
	- All MS/MS: full_spectra - For each consensusElement, the GNPSExport will output All MS/MS scans.
Note that mass accuracy and the retention time window for the pairing between MS/MS scans and a LC-MS feature
orconsensusElement is defined at the IDMapper tool step.
A representative OpenMS-GNPS workflow would sequencially use these OpenMS TOPP tools:
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

#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
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
  TOPPBase("GNPSExport", "Tool to export representative consensus MS/MS scan per consensusElement into a .MGF file format. \
  See the documentation on https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking_with_openms", false) {}

  map<double,double> generateMSMSSpectrumBins(vector<pair<double,double>> mz_int_pairs, double delta_mz)
  {
    //
    // sort input spectrum
    //
    sort(mz_int_pairs.begin(), mz_int_pairs.end(), [](const pair<double,double> &a, const pair<double,double> &b)
    {
      return a.first < b.first;
    });

    //
    // generate new spectrum
    //
    vector<double> mz_merged = new vector<double>();
    vector<double> intensity_merged = new vector<double>();
    double last_mz = numeric_limits<double>::min();
    double sum_mz = 0;
    double sum_intensity = 0;
    Size count = 0;
    for (auto it_mz=mz_int_pairs.begin(); it_mz!=mz_int_pairs.end(); ++it_mz)
    {
      if (abs(it_mz->first-last_mz) > delta_mz && count > 0)
      {
        if (sum_intensity > 0)
        {
          mz_merged.push_back(1.*sum_mz/count);
          intensity_merged.push_back(sum_intensity);
        }

        last_mz = it_mz->first;
        sum_mz = 0;
        count = 0;
        sum_intensity = 0;
      }

      sum_mz += it_mz -> first;
      sum_intensity += it_mz->second;
      count++;
    }
    //remaining scans in last bucket
    if(count > 0)
    {
      mz_merged.push_back(sum_mz/count);
      intensity_merged.push_back(sum_intensity);
    }

    // map mz and intensity
    map<double,double> ms2_block;
    for(Size ms2_block_index = 0; ms2_block_index < mz_merged.size(); ms2_block_index++)
    {
      ms2_block[mz_merged[ms2_block_index]] = intensity_merged[ms2_block_index];
    }

    return ms2_block;
  }

  vector<int> sortElementMapsByIntensity(ConsensusFeature feature)
  {
    // convert element maps to vector of pair<int,double>(map, intensity)
    vector<pair<int,double>> element_ints = new vector<pair<int,double>>();
    for (ConsensusFeature::HandleSetType::const_iterator feature_iter = feature.begin();\
          feature_iter != feature.end(); ++feature_iter)
    {
      element_ints.push_back(pair<int,double>(feature_iter->getMapIndex(), feature_iter->getIntensity()));
    }

    // sort elements by intensity
    sort(element_ints.begin(), element_ints.end(), [](const pair<int,double> &a, const pair<int,double> &b)
    {
      return a.second > b.second;
    });

    // flatten vector of pairs<int,double>(map,intensity) to just vector of element maps
    vector<int> element_maps = new vector<int>();
    for(auto element_int : element_ints)
    {
      element_maps.push_back(element_int.first);
    }
    return element_maps;
  }

  vector<PeptideIdentification> getElementPeptideIdentificationsByElementIntensity(ConsensusFeature feature, vector<int> sorted_element_maps)
  {
    vector<PeptideIdentification> pepts = vector<PeptideIdentification>;
    for (int element_map : sorted_element_maps)
    {
      for (auto pept_id : feature.getPeptideIdentifications())
      {
        if (pept_id.metaValueExists("spectrum_index") && pept_id.metaValueExists("map_index") \
            && (int)pept_id.getMetaValue("map_index") == element_map)
        {
          pepts.push_back(pept_id);
          break;
        }
      }
    }

    return pepts;
  }

private:
  constexpr double static DEF_COSINE_SIMILARITY = 0.9;
  constexpr double static DEF_MERGE_BIN_SIZE = BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES;

  constexpr double static DEF_PREC_MASS_TOL = 0.5;
  constexpr bool static DEF_PREC_MASS_TOL_ISPPM = false;

  constexpr double static DEF_PEPT_CUTOFF = 5;
  constexpr double static DEF_MSMAP_CACHE = 50;
  // double DEF_PRECURSOR_MZ_TOLERANCE = 0.0001;
  // double DEF_PRECURSOR_RT_TOLERANCE = 5;

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in_cm", "<file>", "", "input consensusXML file containing only consensusElements with \"peptide\" annotations.");
    setValidFormats_("in_cm", ListUtils::create<String>("consensusXML"));

    registerInputFileList_("in_mzml", "<files>", ListUtils::create<String>(""), "original mzml files containing the ms2 spectra (aka peptide annotation)");
    setValidFormats_("in_mzml", ListUtils::create<String>("mzML"));

    registerDoubleOption_("peptide_cutoff", "<num>", DEF_PEPT_CUTOFF, "Number of most intense peptides to consider when merging", false);
    registerIntOption_("msmap_cache", "<num>", DEF_MSMAP_CACHE, "Number of msmaps that can be cached during export for optimized performance ()", false);

    registerOutputFile_("out", "<file>", "", "Output MGF file");
    setValidFormats_("out", ListUtils::create<String>("mgf"));

    registerStringOption_("output_type", "<choice>", "merged_spectra", "specificity of mgf output information", false);
    setValidStrings_("output_type", ListUtils::create<String>("full_spectra,merged_spectra,most_intense"));

    registerTOPPSubsection_("merged_spectra", "Options for exporting mgf file with merged spectra per consensusElement");
    registerDoubleOption_("merged_spectra:precursor_mass_tolerance", "<num>", DEF_PREC_MASS_TOL, "Precursor mass tolerance (Da) for ms annotations", false);
    registerDoubleOption_("merged_spectra:cos_similarity", "<num>", DEF_COSINE_SIMILARITY, "Cosine similarity threshold for merged_spectra output", false);
    registerDoubleOption_("merged_spectra:ms2_bin_size", "<num>", DEF_MERGE_BIN_SIZE, "Bin size (Da) for fragment ions when merging ms2 scans", false);
  }

  // the main function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    ProgressLogger progress_logger;
    progress_logger.setLogType(log_type_);

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    double pept_cutoff(getDoubleOption_("peptide_cutoff"));
    double cos_sim_threshold(getDoubleOption_("merged_spectra:cos_similarity"));
    double bin_width(getDoubleOption_("merged_spectra:ms2_bin_size"));

    int max_msmap_cache(getIntOption_("msmap_cache"));

    String consensus_file_path(getStringOption_("in_cm"));
    StringList mzml_file_paths(getStringList_("in_mzml"));
    String out(getStringOption_("out"));
    String output_type(getStringOption_("output_type"));

    ofstream output_file(out);


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    // ConsensusMap
    ConsensusXMLFile consensus_file;
    ConsensusMap consensus_map;
    consensus_file.load(consensus_file_path, consensus_map);

    cout << "output_type..." << output_type << endl; // LOG:

    // vector<String> element_ids = new vector<String>();
    // for (Size cons_i = 0; cons_i < consensus_map.size(); ++cons_i)
    // {

    // }


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    progress_logger.startProgress(0, consensus_map.size(), "parsing features and ms2 identifications...");

    map<int,MSExperiment> ms_maps;

    for (Size cons_i = 0; cons_i < consensus_map.size(); ++cons_i)
    {
      ConsensusFeature feature = consensus_map[cons_i];

      //
      // determine feature's charge
      //
      BaseFeature::ChargeType charge = feature.getCharge();
      for (ConsensusFeature::HandleSetType::const_iterator feature_iter = feature.begin();\
            feature_iter != feature.end(); ++feature_iter)
      {
        // determine feature's charge
        if (feature_iter->getCharge() > charge)
        {
          charge = feature_iter->getCharge();
        }
      }

      //
      // compute most intense peptide identifications (based on precursor intensity)
      //
      vector<int> element_maps = sortElementMapsByIntensity(feature);
      vector<PeptideIdentification> pepts = \
        getElementPeptideIdentificationsByElementIntensity(feature, element_maps);
      if (output_type == "most_intense")
      {
        GNPSExportBlock block(cons_i, feature.getUniqueId(), charge, feature.getMZ(), \
        {pair<int,int>((int)pepts[0].getMetaValue("map_index"),(int)pepts[0].getMetaValue("spectrum_index"))});
        export_blocks.push_back(block);
      }
      else if (output_type == "merged_spectra")
      {
        if (pepts.size() > pept_cutoff)
        {
          pepts.erase(pepts.begin()+pept_cutoff, pepts.end());
        }

        // clean up unused maps from cache
        if (ms_maps.size() > max_msmap_cache)
        {
          for (auto ms_maps_iter=ms_maps.begin();ms_maps_iter!=ms_maps.end();++ms_maps_iter)
          {
            int map_index = ms_maps_iter->first;
            bool should_delete = true;

            for (PeptideIdentification pept : pepts)
            {
              if (map_index == (int) pept.getMetaValue("map_index")) { should_delete = false; }
            }

            if (should_delete)
            {
              cout << "\tdeleting..." << ms_maps_iter->first << endl;
              ms_maps.erase(map_index);
              ms_maps_iter=ms_maps.begin();
            }
          }
        }

        vector<pair<int,int>> map_refs new vector<pair<int,int>>;

        // make sure all maps are added
        for (PeptideIdentification pept : pepts)
        {
          if (ms_maps.find(pept.getMetaValue("map_index")) == ms_maps.end())
          {
            cout << "\tloading " << pept.getMetaValue("map_index") << " map...";
            MzMLFile mzml_file;
            MSExperiment spec_map;
            mzml_file.load(mzml_file_paths[pept.getMetaValue("map_index")], spec_map);

            ms_maps.insert({pept.getMetaValue("map_index"), spec_map});

            cout << "ms_maps.size " << ms_maps.size() << endl;
          }
        }


        BinnedSpectrum binned_highest_int(ms_maps[pepts[0].getMetaValue("map_index")][pepts[0].getMetaValue("spectrum_index")],\
          BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);
        map_refs.push_back(pair<int,int>(pepts[0].getMetaValue("map_index"),pepts[0].getMetaValue("spectrum_index")));
        for (PeptideIdentification pept : pepts)
        {
          auto spec = ms_maps[pept.getMetaValue("map_index")][pept.getMetaValue("spectrum_index")];

          if (pept == pepts[0]) { continue; }

          BinnedSpectrum binned_spectrum(spec,BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES,false,\
            1,BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

          BinnedSpectralContrastAngle bsca;
          if (bsca(binned_highest_int, binned_spectrum) < cos_sim_threshold) { continue; }

          map_refs.push_back(pair<int,int>(pept.getMetaValue("map_index"),pept.getMetaValue("spectrum_index")));
        }


        GNPSExportBlock block(cons_i, feature.getUniqueId(), charge, feature.getMZ(), map_refs);
        export_blocks.push_back(block);

        cout << "pushed back..." << endl;
      }

      progress_logger.setProgress(cons_i);
    }
    progress_logger.endProgress();
    cout << "export_blocks.size()..." << export_blocks.size() << endl; // LOG:


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // compute map reference frequencies for effective "caching"
    map<int,double> map_freq;
    unsigned int total_map_refs = 0;
    for (GNPSExportBlock block : export_blocks)
    {
      total_map_refs += block.map_refs.size();
      for (pair<int,int> block_map : block.map_refs)
      {
        if (map_freq.find(block_map.first) == map_freq.end())
        {
          map_freq.insert({block_map.first, 1});
        }
        else
        {
          map_freq.insert({block_map.first, map_freq[block_map.first]+1});
        }
      }
    }
    for (auto map_freq_iter=map_freq.begin(); map_freq_iter!=map_freq.end(); ++map_freq_iter)
    {
      // map_freq[map_freq_iter->first] = map_freq_iter->second/total_map_refs*1.;
      cout << "\tmap[" << map_freq_iter->first << "]..." << map_freq[map_freq_iter->first] << endl;
    }


    progress_logger.startProgress(0, export_blocks.size(), "exporting " + output_type + " of ms2 identifications...");
    for (Size block_iter=0; block_iter<export_blocks.size(); ++block_iter)
    {
      GNPSExportBlock block = export_blocks[block_iter];
      ConsensusFeature feature = consensus_map[block.scan_id];

      if (output_type == "most_intense")
      {
        //
        // read most_intense scan info
        //
        int map_index = block.map_refs[0].first;
        int spec_index = block.map_refs[0].second;

        if (ms_maps.find(map_index) == ms_maps.end())
        {
          MzMLFile mzml_file;
          MSExperiment spec_map;
          mzml_file.load(mzml_file_paths[map_index], spec_map);

          ms_maps.insert({map_index, spec_map});
        }

        vector<pair<double,double>> mz_int_pairs = new vector<pair<double,double>>();

        for (auto spec_iter=ms_maps[map_index][spec_index].begin(); spec_iter!=ms_maps[map_index][spec_index].end(); ++spec_iter)
        {
          mz_int_pairs.push_back(pair<double,double>(spec_iter->getMZ(),spec_iter->getIntensity()));
        }
        map<double,double> ms2_block(generateMSMSSpectrumBins(mz_int_pairs, -1));

        output_file << "BEGIN IONS" << endl;
        output_file << "OUTPUT=" << output_type << endl;

        output_file << "SCANS=" << (block_iter+1) << endl;
        output_file << "FEATURE_ID=e_" << block.feature_id << endl;

        output_file << "MSLEVEL=2" << endl;
        output_file << "CHARGE=" << to_string(block.charge == 0 ? 1 : block.charge) << "+" << endl;
        output_file << "PEPMASS=" << block.pepmass << endl;
        output_file << "FILE_INDEX=" << spec_index << endl;
        output_file << "RTINSECONDS=" << ms_maps[map_index][spec_index].getRT() << endl;

        output_file << fixed << setprecision(4);
        for (auto ms2_iter = ms2_block.begin(); ms2_iter != ms2_block.end(); ++ms2_iter)
        {
          if ((int) ms2_iter->second > 0)
          {
            output_file << ms2_iter->first << "\t" << (int) ms2_iter->second << endl;
          }
        }

        output_file << "END IONS" << endl << endl;
      }
      else if (output_type == "merged_spectra")
      {
        cout << "outputting maps..." << endl;
        //
        // assert all maps are loaded
        //
        vector<pair<double,double>> mz_int_pairs = new vector<pair<double,double>>();

        for (auto maps_iter=block.map_refs.begin(); maps_iter!=block.map_refs.end(); maps_iter++)
        {
          Size map_index = maps_iter->first;
          Size spec_index = maps_iter->second;

          if (ms_maps.find(map_index) == ms_maps.end())
          {
            MzMLFile mzml_file;
            MSExperiment spec_map;
            mzml_file.load(mzml_file_paths[map_index], spec_map);

            ms_maps.insert({map_index,spec_map});
          }

          for (auto spec_iter=ms_maps[map_index][spec_index].begin();spec_iter!=ms_maps[map_index][spec_index].end(); ++spec_iter)
          {
            mz_int_pairs.push_back(pair<double,double>(spec_iter->getMZ(), spec_iter->getIntensity()));
          }
        }

        map<double,double> ms2_block = generateMSMSSpectrumBins(mz_int_pairs, bin_width);

        output_file << "BEGIN IONS" << endl;
        output_file << "OUTPUT=" << output_type << endl;

        output_file << "SCANS=" << (block_iter+1) << endl;
        output_file << "FEATURE_ID=e_" << block.feature_id << endl;

        output_file << "MSLEVEL=2" << endl;
        output_file << "CHARGE=" << to_string(block.charge == 0 ? 1 : block.charge) << "+" << endl;
        output_file << "PEPMASS=" << block.pepmass << endl;
        output_file << "FILE_INDEX=" << block.map_refs[0].second << endl;
        output_file << "RTINSECONDS=" << ms_maps[block.map_refs[0].first][block.map_refs[0].second].getRT() << endl;

        output_file << fixed << setprecision(4);
        for (auto ms2_iter = ms2_block.begin(); ms2_iter != ms2_block.end(); ++ms2_iter)
        {
          if ((int) ms2_iter->second > 0)
          {
            output_file << ms2_iter->first << "\t" << (int) ms2_iter->second << endl;
          }
        }

        output_file << "END IONS" << endl << endl;
      }

      progress_logger.setProgress(block_iter);
    }
    progress_logger.endProgress();

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
