// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <string>
#include <iostream>

#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/util/XMLString.hpp>

namespace OpenMS
{

  namespace StringUtils
  {
    std::streampos stringToStreampos(std::string s)
    {
      // Try to cast the string to a type that can hold the integer value
      // stored in the std::streampos type (which can vary from 32 to 128 bits
      // depending on the system).
      //
      // long long has a minimal size of 64 bits and can address a range up to
      // 16 Exbibit (or 2 Exbibyte), we can hopefully expect our files to be
      // smaller than an Exabyte and should be safe.
      std::streampos res;
      try
      {

        res = boost::lexical_cast< unsigned long long >(s);
        // TESTING CAST: res = (int) res; // use this only when emulating 32 bit systems

      }
      catch (boost::bad_lexical_cast&)
      {
        std::cerr << "Trying to convert corrupted / unreadable value to std::streampos : " << s << std::endl;
        std::cerr << "This can also happen if the value exceeds 63 bits, please check your input." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            String("Could not convert string '") + s + "' to a 64 bit integer.");
      }

      // Check if the value can fit into std::streampos
      if ( fabs( boost::lexical_cast< long double >(s) - res) > 0.1)
      {
        std::cerr << "Your system may not support addressing a file of this size,"
          << " only addresses that fit into a " << sizeof(std::streamsize)*8 <<
          " bit integer are supported on your system." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            String("Could not convert string '") + s + "' to an integer on your system.");
      }

      return res;
    }
  }

  int IndexedMzMLDecoder::parseOffsets(String filename, std::streampos indexoffset, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {
    //-------------------------------------------------------------
    // Open file, jump to end and read last indexoffset bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(filename.c_str());
    // get length of file:
    f.seekg(0, f.end);
    std::streampos length = f.tellg();

    if (indexoffset < 0 || indexoffset > length)
    {
      std::cerr << "IndexedMzMLDecoder::parseOffsets Error: Offset was " <<
        indexoffset << " (not between 0 and " << length << ")." << std::endl;
      return -1;
    }

    //-------------------------------------------------------------
    // Read full end of file to parse offsets for spectra and chroms
    //-------------------------------------------------------------
    // read data as a block:
    // allocate memory:
    std::streampos readl = length - indexoffset;
    char* buffer = new char[ readl + std::streampos(1)];
    f.seekg(-readl, f.end);
    f.read(buffer, readl);
    buffer[readl] = '\0';

    //-------------------------------------------------------------
    // Add a sane start element and then give it to a DOM parser
    //-------------------------------------------------------------
    // http://stackoverflow.com/questions/4691039/making-xerces-parse-a-string-insted-of-a-file
    std::string tmp_fixed_xml = "<indexedmzML>" +  String(buffer) + "\n";
    int res = domParseIndexedEnd_(tmp_fixed_xml, spectra_offsets, chromatograms_offsets);

    delete[] buffer;

    return res;
  }

  std::streampos IndexedMzMLDecoder::findIndexListOffset(String filename, int buffersize)
  {
    // return value
    std::streampos indexoffset = -1;

    //-------------------------------------------------------------
    // Open file, jump to end and read last n bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(filename.c_str());

    if (!f.is_open())
    {
      return indexoffset;
    }

    // Read the last few bytes and hope our offset is there to be found
    char* buffer = new char[buffersize + 1];
    f.seekg(-buffersize, f.end);
    f.read(buffer, buffersize);
    buffer[buffersize] = '\0';

#ifdef DEBUG_READER
    std::cout << " reading file " << filename  << " with size " << buffersize << std::endl;
    std::cout << buffer << std::endl;
#endif

    //-------------------------------------------------------------
    // Since we could be anywhere in the XML structure, use regex to find
    // indexListOffset and read its content.
    //-------------------------------------------------------------
    boost::regex listoffset_rx("<[^>/]*indexListOffset\\s*>\\s*(\\d*)");
    boost::cmatch matches;
    boost::regex_search(buffer, matches, listoffset_rx);
    String thismatch(matches[1].first, matches[1].second);
    if (thismatch.size() > 0)
    {
      try
      {
        indexoffset = OpenMS::StringUtils::stringToStreampos(thismatch);
      }
      catch (Exception::ConversionError& e)
      {
        std::cerr << "Corrupted / unreadable value in <indexListOffset> : " << thismatch << std::endl;
        // free resources and re-throw
        delete[] buffer;
        f.close();
        throw;
      }
    }
    else
    {
      std::cerr << "IndexedMzMLDecoder::findIndexListOffset Error: Could not find element indexListOffset in the last " <<
        buffersize << " bytes. Maybe this is not a indexedmzML." << std::endl;
      std::cerr << buffer << std::endl;
    }

    f.close();
    delete[] buffer;

    return indexoffset;
  }

  int IndexedMzMLDecoder::domParseIndexedEnd_(std::string in, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {
    //-------------------------------------------------------------
    // Create parser from input string using MemBufInputSource
    //-------------------------------------------------------------
    xercesc::MemBufInputSource myxml_buf(reinterpret_cast<const unsigned char*>(in.c_str()), in.length(), "myxml (in memory)");
    xercesc::XercesDOMParser* parser = new xercesc::XercesDOMParser();
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD(false);
    parser->parse(myxml_buf);

    //-------------------------------------------------------------
    // Start parsing
    // see http://www.yolinux.com/TUTORIALS/XML-Xerces-C.html
    //-------------------------------------------------------------

    // no need to free this pointer - owned by the parent parser object
    xercesc::DOMDocument* doc =  parser->getDocument();
    // Get the top-level element ("indexedmzML")
    xercesc::DOMElement* elementRoot = doc->getDocumentElement();
    if (!elementRoot)
    {
      std::cerr << "IndexedMzMLDecoder::domParseIndexedEnd Error: No root element found:" << std::endl << std::endl << in << std::endl;
      delete parser;
      return -1;
    }

    // Extract the indexList tag (there should only be one)
    XMLCh* tag = xercesc::XMLString::transcode("indexList");
    xercesc::DOMNodeList* li = elementRoot->getElementsByTagName(tag);
    xercesc::XMLString::release(&tag);
    if (li->getLength() != 1)
    {
      std::cerr << "IndexedMzMLDecoder::domParseIndexedEnd Error: no indexList element found:" << std::endl << std::endl << in << std::endl;
      delete parser;
      return -1;
    }
    xercesc::DOMNode* indexListNode = li->item(0);

    XMLCh* idref_tag = xercesc::XMLString::transcode("idRef");
    XMLCh* name_tag = xercesc::XMLString::transcode("name");

    // Iterate through indexList elements
    xercesc::DOMNodeList* index_elems = indexListNode->getChildNodes();
    const  XMLSize_t nodeCount_ = index_elems->getLength();
    for (XMLSize_t j = 0; j < nodeCount_; ++j)
    {
      xercesc::DOMNode* currentNode = index_elems->item(j);
      if (currentNode->getNodeType() && // true is not NULL
          currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE) // is element
      {
        std::vector<std::pair<std::string, std::streampos> > result;
        xercesc::DOMNodeList* offset_elems = currentNode->getChildNodes();
        for (XMLSize_t k = 0; k < offset_elems->getLength(); ++k)
        {
          xercesc::DOMNode* currentONode = offset_elems->item(k);
          if (currentONode->getNodeType() && // true is not NULL
              currentONode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE) // is element
          {
            xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentONode);

            char* name = xercesc::XMLString::transcode(currentElement->getAttribute(idref_tag));
            char* offset = xercesc::XMLString::transcode(currentONode->getTextContent());

            std::streampos thisOffset = OpenMS::StringUtils::stringToStreampos( String(offset) );
            result.push_back(std::make_pair(std::string(name), thisOffset));

            xercesc::XMLString::release(&name);
            xercesc::XMLString::release(&offset);
          }
        }

        // should be either spectrum or chromatogram ...
        xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentNode);
        char* tmp_name = xercesc::XMLString::transcode(currentElement->getAttribute(name_tag));
        std::string name(tmp_name);
        xercesc::XMLString::release(&tmp_name);

        if (name == "spectrum")
        {
          spectra_offsets = result;
        }
        else if (name == "chromatogram")
        {
          chromatograms_offsets = result;
        }
        else
        {
          std::cerr << "IndexedMzMLDecoder::domParseIndexedEnd Error: expected only " <<
            "'spectrum' or 'chromatogram' below indexList but found instead '" << name << "'." << std::endl;
          xercesc::XMLString::release(&idref_tag);
          xercesc::XMLString::release(&name_tag);
          delete parser;
          return -1;
        }
      }
    }
    xercesc::XMLString::release(&idref_tag);
    xercesc::XMLString::release(&name_tag);


    delete parser;
    return 0;
  }

}